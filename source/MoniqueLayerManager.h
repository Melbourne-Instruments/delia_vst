#pragma once
#include "MoniqueCalibrator.h"
#include "MoniqueCommon.h"
#include "MoniqueEffectsEngine.h"
#include "MoniqueLayer.h"
#include "common.h"
#include <pthread.h>

namespace Monique {

static void *_run_fx_synth(void *data);

class LayerManager {
  public:
    LayerManager() {
        constexpr int FX_RT_THREAD_PRIORITY = 75;
        constexpr int FX_RT_THREAD_CPU_AFINITY = 2;

        printf("\n init the layer manager");

        // setup layer 0 with max voices and layer 1 with zero
        _layer_0.updateParameter(gen_param_id(PresetCommonParameters::LAYER_1_NUM_VOICES), 1);
        _layer_1.updateParameter(gen_param_id(PresetCommonParameters::LAYER_2_NUM_VOICES), 0);

        // Create the FX RT thread condition variable and mutex
        __cobalt_pthread_cond_init(&_fx_synth_cond, nullptr);
        __cobalt_pthread_mutex_init(&_fx_synth_mutex, nullptr);

        // Initialise the FX RT thread attributes
        pthread_attr_t task_attributes;
        __cobalt_pthread_attr_init(&task_attributes);
        pthread_attr_setdetachstate(&task_attributes, PTHREAD_CREATE_JOINABLE);
        pthread_attr_setinheritsched(&task_attributes, PTHREAD_EXPLICIT_SCHED);
        pthread_attr_setschedpolicy(&task_attributes, SCHED_FIFO);

        // Set the FX RT thread affinity to CPU core 2
        cpu_set_t cpus;
        CPU_ZERO(&cpus);
        CPU_SET(FX_RT_THREAD_CPU_AFINITY, &cpus);
        pthread_attr_setaffinity_np(&task_attributes, sizeof(cpu_set_t), &cpus);

        // Create and spawn the FX RT thread
        _exit_fx_rt_thread = false;
        _fx_run = false;
        _fx_done = false;
        struct sched_param rt_params = {.sched_priority = FX_RT_THREAD_PRIORITY};
        pthread_attr_setschedparam(&task_attributes, &rt_params);
        __cobalt_pthread_create(&_fx_rt_thread, &task_attributes, _run_fx_synth, this);
    }

    ~LayerManager(){};

    void handleNote(MidiNote note_on) {
        _layer_0.handleNote(note_on);
        _layer_1.handleNote(note_on);
        // printf("\n note %d", note_on.pitch);
    }

    void polyPressure(PolyPressure pp) {

        _layer_0.polyPressure(pp);
        _layer_1.polyPressure(pp);
    }

    void updateParameter(const uint param_id, const float value) {
        _in_change_ids.push_back(param_id);
        _in_change_values.push_back(value);
    }

    std::vector<ParamChange> &paramRefresh() {
        float denormalised;
        for (int i = 0; i < _in_change_ids.size(); i++) {
            uint param_id = _in_change_ids.at(i);
            float value = _in_change_values.at(i);
            switch (param_id) {
            case (int)GlobalParams::RUN_FILTER_CAL:
                if (value < 0.25) {
                    _run_calibration = false;

                } else {
                    _calibrator.startFilterCal(true);
                    _run_calibration = true;
                }
                break;
            case (int)GlobalParams::RUN_VCA_CAL:
                if (value < 0.25) {
                    _run_calibration = false;
                } else {
                    // old cal
                    //_calibrator.startMainVcaCal(true);

                    // new cal
                    _calibrator.startVcaCal(true);
                    _run_calibration = true;
                }
                break;

            case (int)GlobalParams::RUN_DC_CAL:
                if (value < 0.25) {
                    _run_calibration = false;
                } else {
                    _calibrator.startDcCal(true);
                    _run_calibration = true;
                }
                break;
            case (int)GlobalParams::RUN_DC_OD_CAL:
                if (value < 0.25) {
                    _run_calibration = false;
                } else {
                    _calibrator.startDcCalOverdrive(true);
                    _run_calibration = true;
                }
                break;
            case (int)GlobalParams::RELOAD_CAL:
                _layer_0.reloadCal();

                _layer_1.reloadCal();
                _calibrator.reloadCal();
                break;
            case (int)GlobalParams::ALL_NOTES_OFF:
                _layer_0.allNotesOff();
                _layer_1.allNotesOff();
                break;
            case (int)GlobalParams::CV_1_MODE: {
                denormalised = from_normalised_float(GlobalParams::CV_1_MODE, value);
                _parameters.cv_1_mode = (CvInputMode)(int)denormalised;
                float offset = _parameters.cv_1_offset_setting;
                float gain = _parameters.cv_1_gain_setting;
                cvSwitchHelper(_parameters.cv_1_mode, offset, gain);
                _cv_1_gain = gain;

                _cv_1_offset = offset;
            } break;
            case (int)GlobalParams::CV_2_MODE: {
                denormalised = from_normalised_float(GlobalParams::CV_1_MODE, value);

                _parameters.cv_2_mode = (CvInputMode)(int)denormalised;
                float offset = _parameters.cv_2_offset_setting;
                float gain = _parameters.cv_2_gain_setting;
                cvSwitchHelper(_parameters.cv_2_mode, offset, gain);
                _cv_2_gain = gain;
                _cv_2_offset = offset;
            } break;
            case (int)GlobalParams::CV_1_GAIN: {
                denormalised = from_normalised_float(GlobalParams::CV_1_GAIN, value);
                _parameters.cv_1_gain_setting = denormalised;
                float offset = _parameters.cv_1_offset_setting;
                float gain = _parameters.cv_1_gain_setting;
                cvSwitchHelper(_parameters.cv_1_mode, offset, gain);
                _cv_1_gain = gain;
                _cv_1_offset = offset;

            } break;
            case (int)GlobalParams::CV_2_GAIN: {
                denormalised = from_normalised_float(GlobalParams::CV_2_GAIN, value);
                _parameters.cv_2_gain_setting = denormalised;
                float offset = _parameters.cv_2_offset_setting;
                float gain = _parameters.cv_2_gain_setting;
                cvSwitchHelper(_parameters.cv_2_mode, offset, gain);
                _cv_2_gain = gain;
                _cv_2_offset = offset;

            } break;
            case (int)GlobalParams::CV_1_OFFSET: {
                denormalised = from_normalised_float(GlobalParams::CV_1_OFFSET, value);
                _parameters.cv_1_offset_setting = denormalised;
                float offset = _parameters.cv_1_offset_setting;
                float gain = _parameters.cv_2_gain_setting;
                cvSwitchHelper(_parameters.cv_1_mode, offset, gain);
                _cv_1_gain = gain;
                _cv_1_offset = offset;

            } break;
            case (int)GlobalParams::CV_2_OFFSET: {
                denormalised = from_normalised_float(GlobalParams::CV_2_OFFSET, value);
                _parameters.cv_2_offset_setting = denormalised;

                float offset = _parameters.cv_2_offset_setting;
                float gain = _parameters.cv_2_gain_setting;
                cvSwitchHelper(_parameters.cv_2_mode, offset, gain);
                _cv_2_gain = gain;
                _cv_2_offset = offset;

            } break;
            default:
                const ParamDecoded param = ParamDecoded(Parameter(param_id));
                _layer_0.updateParameter(param, value);
                _layer_1.updateParameter(param, value);
                _effects.processParamChange(param, param_id, value);
                break;
            }
        }
        _in_change_ids.clear();
        _in_change_values.clear();

        _layer_0.paramRefresh();
        _layer_1.paramRefresh();
        return _morph_changes;
    }

    std::vector<ParamChange> &getChanges() {
        return _morph_changes;
    }

    void setTempo(float tempo) {
        _global_tempo = (tempo) / 60.0;

        _layer_0.setTempo(_global_tempo);
        _layer_1.setTempo(_global_tempo);
        _effects.setTempo(_global_tempo);
    };

    void cvSwitchHelper(const CvInputMode &mode, float &offset_in, float &gain_in) {
        constexpr float offset_max = 0.12;
        constexpr float gain_max = 0.12;
        constexpr float _0_10_gain = -1.0;
        constexpr float _10_10_gain = _0_10_gain / 2;

        constexpr float _10_10_offset = 0;
        constexpr float _0_10_offset = +0.25;
        constexpr float _5_5_gain = _0_10_gain;
        constexpr float _5_5_offset = 0;
        constexpr float _0_5_gain = _0_10_gain * 2;
        constexpr float _0_5_offset = _0_10_offset / 2;
        const float offset = offset_in;
        const float gain = gain_in;
        switch (mode) {
        case CvInputMode::NEG_10_TO_10:
            offset_in = (offset * offset_max);
            gain_in = _10_10_gain * (1 + gain * gain_max);
            break;
        case CvInputMode::ZERO_TO_10:
            offset_in = _0_10_offset + (offset * offset_max);
            gain_in = _0_10_gain * (1 + gain * gain_max);
            break;
        case CvInputMode::ZERO_TO_5:

            offset_in = _0_5_offset + (offset * offset_max);
            gain_in = _0_5_gain * (1 + gain * gain_max);
            break;
        case CvInputMode::NEG_5_TO_5:
            offset_in = (offset * offset_max);
            gain_in = _5_5_gain * (1 + gain * gain_max);
            break;

        default:
            break;
        }
        return;
    }

    void runSynth(ProcessAudioData &data) {
        if (_run_calibration) {
            _calibrator.run(data);
            return;
        }

        // Signal the FX run condition - to process the FX synth data
        __cobalt_pthread_mutex_lock(&_fx_synth_mutex);
        _audio_data = &data;
        _fx_run = true;
        __cobalt_pthread_cond_signal(&_fx_synth_cond);
        __cobalt_pthread_mutex_unlock(&_fx_synth_mutex);

        // recalculate CV inputs
        constexpr int cv_step = SAMPLE_RATE / CV_SAMPLE_RATE;
        float s1, s2;
        for (uint i = 0; i < CV_BUFFER_SIZE; i++) {
            s1 = ((*data.input_1)[i * cv_step]);
            s2 = ((*data.input_1)[i * cv_step]);
            _parameters._cv_1[i] = _cv_1_gain * ((*data.input_1)[i * cv_step] + _cv_1_offset);
            _parameters._cv_2[i] = _cv_2_gain * ((*data.input_2)[i * cv_step] + _cv_2_offset);
        }
        // printf("\n cv %f %f %f %f", _parameters._cv_1[0], s1, _cv_1_gain, _cv_1_offset);

        // Process the layers
        _layer_0.runSynth(data);
        _layer_1.runSynth(data);

        // Wait for the FX done condition to be signalled
        __cobalt_pthread_mutex_lock(&_fx_synth_mutex);
        while (!_fx_done) {
            __cobalt_pthread_cond_wait(&_fx_synth_cond, &_fx_synth_mutex);
        }
        _fx_done = false;
        __cobalt_pthread_mutex_unlock(&_fx_synth_mutex);

        // Set the FX modulation destination values
        float fx_snd_mod, fx_macro_mod;
        _layer_0.getFxMod(fx_snd_mod, fx_macro_mod);
        _effects.setModDests(fx_snd_mod, fx_macro_mod);
    }

    void runFxSynth() {
        // Run forever until exited
        while (!_exit_fx_rt_thread) {
            // Wait for the FX run condition to be signalled
            __cobalt_pthread_mutex_lock(&_fx_synth_mutex);
            while (!_fx_run) {
                __cobalt_pthread_cond_wait(&_fx_synth_cond, &_fx_synth_mutex);
            }
            _fx_run = false;
            __cobalt_pthread_mutex_unlock(&_fx_synth_mutex);

            // Process the effects data
            if (_audio_data) {
                _effects.run(*_audio_data->input_l, *_audio_data->input_r, *_audio_data->output_l, *_audio_data->output_r);
            }

            // Signal the FX done condition - FX synth data has been processed
            __cobalt_pthread_mutex_lock(&_fx_synth_mutex);
            _fx_done = true;
            __cobalt_pthread_cond_signal(&_fx_synth_cond);
            __cobalt_pthread_mutex_unlock(&_fx_synth_mutex);
        }

        // Thread exited, destroy the condition variable and mutex
        __cobalt_pthread_cond_destroy(&_fx_synth_cond);
        __cobalt_pthread_mutex_destroy(&_fx_synth_mutex);
    }

    void stopSynth() {
        // Stop the FX RT thread
        __cobalt_pthread_mutex_lock(&_fx_synth_mutex);
        _exit_fx_rt_thread = true;
        _audio_data = nullptr;
        _fx_run = true;
        __cobalt_pthread_cond_signal(&_fx_synth_cond);
        __cobalt_pthread_mutex_unlock(&_fx_synth_mutex);
        __cobalt_pthread_join(_fx_rt_thread, nullptr);
        _fx_rt_thread = 0;
    }

    GlobalSynthParameters _parameters;
    float _cv_1_gain = 1;
    float _cv_1_offset = 0;
    float _cv_2_gain = 1;
    float _cv_2_offset = 0;
    std::vector<ParamChange> _morph_changes;
    std::vector<float> _in_change_values;
    std::vector<uint> _in_change_ids;
    std::array<int, NUM_VOICES> _num_voices_array = {3, 3};
    bool _run_calibration = false;
    MoniqueCalSequencer _calibrator;
    MoniqueFx _effects;
    MoniqueLayer _layer_0 = MoniqueLayer(0, _morph_changes, _num_voices_array, _parameters);
    MoniqueLayer _layer_1 = MoniqueLayer(1, _morph_changes, _num_voices_array, _parameters);
    pthread_t _fx_rt_thread;
    bool _exit_fx_rt_thread;
    pthread_cond_t _fx_synth_cond;
    pthread_mutex_t _fx_synth_mutex;
    bool _fx_run;
    bool _fx_done;
    ProcessAudioData *_audio_data = nullptr;
    float _global_tempo = 120.f / 60.f;
};

//----------------------------------------------------------------------------
// _run_fx_synth
//----------------------------------------------------------------------------
static void *_run_fx_synth(void *data) {
    auto layer_manager = static_cast<LayerManager *>(data);
    layer_manager->runFxSynth();

    // To suppress warnings
    return nullptr;
}

} // namespace Monique