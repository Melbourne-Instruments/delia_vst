#pragma once

#include "MoniqueAnalogModel.h"
#include "MoniqueCommon.h"
#include "SynthMath.h"
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace Monique {
static std::string CAL_FILE_DONE_PATH = "/udata/delia/tuning/cal_complete";

static std::string ssystem(const char *command) {
    char tmpname[L_tmpnam];
    auto res = std::tmpnam(tmpname);
    std::string scommand = command;
    std::string cmd = scommand + " >> " + tmpname;
    auto res_2 = std::system(cmd.c_str());
    std::ifstream file(tmpname, std::ios::in | std::ios::binary);
    std::string result;
    if (file) {
        while (!file.eof())
            result.push_back(file.get());
        file.close();
    }

    // print errors
    else {
        std::cout << res << res_2;
    }
    remove(tmpname);
    return result;
}

enum class VCA_State {
    BEGIN,
    CLICKWAIT,
    RUNNING,
    CLICKWAIT_2,
    RUNNING_2,
    FINISH,
};

enum class DC_State {
    BEGIN,
    CLICKWAIT,
    RUNNING,
    CLICKWAIT_2,
    RUNNING_2,
    FINISH,
    TEST
};

class BiquadFilter2 {
  public:
    BiquadFilter2(float a1, float a2, float b0, float b1, float b2) :
        _a1(a1), _a2(a2), _b0(b0), _b1(b1), _b2(b2) {
        // Reset the filer
        reset();
    };

    ~BiquadFilter2() {
    }

    inline void reset() {
        // Reset the filter
        _m1 = 0;
        _m2 = 0;
        _dn = 1e-20f;
    };

    inline float process(float input) {
        float w = input - (_a1 * _m1) - (_a2 * _m2) + _dn;
        float out = (_b1 * _m1) + (_b2 * _m2) + (_b0 * w);
        _m2 = _m1;
        _m1 = w;
        _dn = -_dn;
        return out;
    }

  private:
    const float _a1, _a2, _b0, _b1, _b2;
    float _m1, _m2;
    float _dn;
};

template <uint FC>
class AnalysisLowpassFilter {
  public:
    AnalysisLowpassFilter() :
        _biquad_1((_a1 / _a0), (_a2 / _a0), (_b0 / _a0), (_b1 / _a0), (_b2 / _a0)),
        _biquad_2((_a1 / _a0), (_a2 / _a0), (_b0 / _a0), (_b1 / _a0), (_b2 / _a0)),
        _biquad_3((_a1 / _a0), (_a2 / _a0), (_b0 / _a0), (_b1 / _a0), (_b2 / _a0)),
        _biquad_4((_a1 / _a0), (_a2 / _a0), (_b0 / _a0), (_b1 / _a0), (_b2 / _a0)){};

    inline void reset() {
        // Reset the biquad filters
        _biquad_1.reset();
        _biquad_2.reset();
        _biquad_3.reset();
        _biquad_4.reset();
    }

    ~AnalysisLowpassFilter() {
        printf("\n");
    }

    inline float process(float input) {
        // Process the biquad filters in sequence
        return _biquad_4.process(_biquad_3.process(_biquad_2.process(_biquad_1.process(input))));
    }

  private:
    static constexpr float _q = 0.70f;
    static constexpr float _w0 = (2 * M_PIf32 * FC) / SAMPLE_RATE;
    static constexpr float _alpha = std::sin(_w0) / (2 * _q);
    static constexpr float _b0 = (1 - std::cos(_w0)) / 2;
    static constexpr float _b1 = 1 - std::cos(_w0);
    static constexpr float _b2 = (1 - std::cos(_w0)) / 2;
    static constexpr float _a0 = 1 + _alpha;
    static constexpr float _a1 = -2 * std::cos(_w0);
    static constexpr float _a2 = 1 - _alpha;
    BiquadFilter2 _biquad_1;
    BiquadFilter2 _biquad_2;
    BiquadFilter2 _biquad_3;
    BiquadFilter2 _biquad_4;
};

class biquadBandpass {
  public:
    biquadBandpass(){};
    ~biquadBandpass(){};

    void reset() {

        _x1 = 0;
        _y1 = 0;
        _x2 = 0;
        _y2 = 0;
    };

    void setup(float freq, float q) {
        float dbGain = 0;
        float A = std::pow(10, dbGain / 40);
        float omega = 2 * M_PI * freq / _sample_rate;
        float sn = std::sin(omega);
        float cs = std::cos(omega);
        float alpha = sn / (2 * q);
        float beta = std::sqrt(A + A);
        _b0 = alpha;
        _b1 = 0;
        _b2 = -alpha;
        _a0 = 1 + alpha;
        _a1 = -2 * cs;
        _a2 = 1 - alpha;
        _b0 /= _a0;
        _b1 /= _a0;
        _b2 /= _a0;
        _a1 /= _a0;
        _a2 /= _a0;
    }

    float run(float x) {
        float y = _b0 * x + _b1 * _x1 + _b2 * _x2 - _a1 * _y1 - _a2 * _y2;
        _x2 = _x1;
        _x1 = x;
        _y2 = _y1;
        _y1 = y;
        return y;
    }

  private:
    const float _sample_rate = SAMPLE_RATE;
    float _a0, _a1, _a2, _b0, _b1, _b2;
    float _x1 = 0;
    float _x2 = 0;
    float _y1 = 0;
    float _y2 = 0;
};

class VoiceCalController {
  public:
    MoniqueAnalogModel &_voice;
    VoiceData &voice_input;
    float _pitch_1 = 6. / NOTE_GAIN;
    float _pitch_2 = 6. / NOTE_GAIN;
    float _filter = 0;
    float _res = 0;
    float _mix_tri1 = 0;
    float _mix_sq1 = 0;
    float _mix_tri2 = 0;
    float _mix_sqr2 = 0;
    float _mix_xor = 0;
    float _mix_l = 0;
    float _mix_r = 0;
    bool _mute = true;
    float _shape_1 = 0;
    float _shape_2 = 0;
    bool _overdrive = false;

    VoiceCalController(MoniqueAnalogModel &voice, VoiceData &v_in) :
        _voice(voice),
        voice_input(v_in) {
    }

    ~VoiceCalController(){};

    void run() {
        voice_input.filt_2_cut.fill(_filter);
        voice_input.filt_2_res.fill(_res);
        voice_input.tri_1_lev.fill(_mix_tri1);
        voice_input.sqr_1_lev.fill(_mix_sq1);
        voice_input.tri_2_lev.fill(_mix_tri2);
        voice_input.sqr_2_lev.fill(_mix_sqr2);
        voice_input.xor_lev.fill(_mix_xor);
        voice_input.vca_l.fill(_mix_l);
        voice_input.vca_r.fill(_mix_r);
        voice_input.mute_1 = (_mute);
        voice_input.mute_2 = (_mute);

        voice_input.osc_1_pitch.fill(_pitch_1);
        voice_input.osc_2_pitch.fill(_pitch_2);
        voice_input.osc_1_shape.fill(_shape_1);
        voice_input.osc_2_shape.fill(_shape_2);
        voice_input.sub_osc = (false);
        voice_input.hard_sync = (false);
        voice_input.overdrive = _overdrive;
        _voice.setMute(_mute);
    };

    void setOff() {
        _filter = .8;
        _res = 1;
        _mix_l = 1;
        _mix_r = 1;
        _mix_tri1 = 0;
        _mix_tri2 = 0;
        _mix_sq1 = 0;
        _mix_sqr2 = 0;
        _mix_xor = 0;
        _mute = true;
        _pitch_1 = 6. / NOTE_GAIN;
        _pitch_2 = 6. / NOTE_GAIN;
        _shape_1 = 0;
        _shape_2 = 0;
    }

    void setFilterTuneTest(float cut) {
        _filter = cut;
        _res = 1;
        _mix_l = 1;
        _mix_r = 1;
        _mix_tri1 = 0;
        _mix_tri2 = 0;
        _mix_sq1 = 0;
        _mix_sqr2 = 0;
        _mix_xor = 0;
        _mute = false;
        _pitch_1 = 6. / NOTE_GAIN;
        _pitch_2 = 6. / NOTE_GAIN;
        _shape_1 = 0;
        _shape_2 = 0;
    }

    void setDCtest(bool overdrive) {
        _filter = .9;
        _res = 0;
        _mix_l = 0;
        _mix_r = 0;
        _mix_tri1 = 0;
        _mix_tri2 = 0;
        _mix_sq1 = 0;
        _mix_sqr2 = 0;
        _mix_xor = 0;
        _mute = false;
        _pitch_1 = 6. / NOTE_GAIN;
        _pitch_2 = 6. / NOTE_GAIN;
        _shape_1 = 0;
        _shape_2 = 0;
        _overdrive = overdrive;
    }

    void setMainOutVca2(float left, float right) {
        _filter = .8;
        _res = 0;
        _mix_l = left;
        _mix_r = right;
        _mix_tri1 = 0;
        _mix_tri2 = 0;
        _mix_sq1 = 0;
        _mix_sqr2 = 0;
        _mix_xor = 0.0;
        _mute = false;
        _pitch_1 = 6. / NOTE_GAIN;
        _pitch_2 = 6. / NOTE_GAIN;
        _shape_1 = 0;
        _shape_2 = 0;
    }

    void setMainOutVca(float left, float right) {
        _filter = .9;
        _res = 1;
        _mix_l = left;
        _mix_r = right;
        _mix_tri1 = 0;
        _mix_tri2 = 0;
        _mix_sq1 = 0;
        _mix_sqr2 = 0;
        _mix_xor = 0.0;
        _mute = false;
        _pitch_1 = 6. / NOTE_GAIN;
        _pitch_2 = 6. / NOTE_GAIN;
        _shape_1 = 0;
        _shape_2 = 0;
    }
};

class MoniqueCalSequencer {
  public:
    std::array<VoiceData, NUM_VOICES> _voice_data;
    std::array<MoniqueAnalogModel, NUM_VOICES> _analog_voices = {
        MoniqueAnalogModel(0, _voice_data.at(0)),
        MoniqueAnalogModel(1, _voice_data.at(1)),
        MoniqueAnalogModel(2, _voice_data.at(2)),
        MoniqueAnalogModel(3, _voice_data.at(3)),
        MoniqueAnalogModel(4, _voice_data.at(4)),
        MoniqueAnalogModel(5, _voice_data.at(5))};

    std::array<VoiceCalController, NUM_VOICES> _controller = {
        VoiceCalController(_analog_voices[0], _voice_data[0]),
        VoiceCalController(_analog_voices[1], _voice_data[1]),
        VoiceCalController(_analog_voices[2], _voice_data[2]),
        VoiceCalController(_analog_voices[3], _voice_data[3]),
        VoiceCalController(_analog_voices[4], _voice_data[4]),
        VoiceCalController(_analog_voices[5], _voice_data[5])};

    std::ofstream _out;
    DC_State dc_state = DC_State::BEGIN;
    VCA_State vca_state = VCA_State::BEGIN;
    std::ofstream _out_2;
    std::vector<float> _audio_data_out;
    std::vector<float> _stimulus_out;
    std::vector<float> _rms_1;
    std::vector<float> _rms_2;
    std::vector<float> _stim_l;
    std::vector<float> _stim_r;
    std::vector<float> _dc_cal_status;
    uint _call_counter = 0;
    AnalysisLowpassFilter<2100> _lowpass_l;
    AnalysisLowpassFilter<1900> _highpass_l;
    AnalysisLowpassFilter<2100> _lowpass_r;
    AnalysisLowpassFilter<2000> _highpass_r;
    AnalysisLowpassFilter<800> _highpass_mix;
    biquadBandpass _bandpass_filter;
    std::vector<float> _audio_dft_buffer;
    uint dft_size = 4906;
    uint N = 4096;
    uint dft_buf_ptr = 0;
    float _shape_1 = 0;
    float _shape_2 = 0;
    int _test_startup_counter = 0;

    void setupMainVcaCal(uint voicen) {
        std::stringstream file;
        file << "/udata/delia/tuning/"
             << "voice_" << voicen << "_main_vca.dat";
        _out.open(file.str(), std::ios::out | std::ios::trunc | std::ios::binary);
        for (uint i = 0; i < NUM_VOICES; i++) {
            _controller.at(i).setOff();
        }
        _highpass_l.reset();
        _lowpass_l.reset();
        _highpass_r.reset();
        _lowpass_r.reset();
    }

    void setupFilter(uint voicen) {
        std::stringstream file;
        file << "/udata/delia/tuning/"
             << "voice_" << voicen << "_filter.dat";
        std::cout << file.str();
        _out.open(file.str(), std::ios::out | std::ios::trunc | std::ios::binary);
        file.str("");

        file << "/udata/delia/tuning/"
             << "voice_" << voicen << "_signal.dat";
        std::cout << file.str();
        _out_2.open(file.str(), std::ios::out | std::ios::trunc | std::ios::binary);
        for (uint i = 0; i < NUM_VOICES; i++) {
            _controller.at(i).setOff();

            // TODO add some kind of defaults here
            // auto &model = _analog_voices.at(i).getModelData();
        }
    }

    void setVoiceShape(float shape_1, float shape_2) {
        _shape_1 = shape_1;
        _shape_2 = shape_2;
    };

    void runMainVca(ProcessAudioData data) {
        if (test_counter == NUM_VOICES) {
            if (!_test_complete) {
                _test_complete = true;
                printf("\n save all cal");
                for (auto &voice : _analog_voices) {
                    voice.saveCalibration();
                }
                std::ofstream outfile;
                outfile.open(CAL_FILE_DONE_PATH, std::ios_base::app); // append instead of overwrite
                outfile << " ";
            }
            return;
        }
        for (auto &a : data.voice_high_res) {
            a->fill(0.f);
        }

        if (buffer_count == 0) {
            _analog_voices.at(test_counter).setMainLVcaOffset(0);
            _analog_voices.at(test_counter).setMainRVcaOffset(0);
            min_l = 1000000;
            min_r = 1000000;
            _lowpass_l.reset();
            _highpass_l.reset();
            _lowpass_r.reset();
            _highpass_r.reset();
            setupMainVcaCal(test_counter);
            for (auto &val : queue_l) {
                val = 1;
            }
            for (auto &val : queue_r) {
                val = 1;
            }
            queue_i = 0;
            old_level_l = 100;
            old_level_r = 100;
            vca_state_l = 0;
            vca_state_r = 0;
        }

        for (int i = 0; i < BUFFER_SIZE; i++) {
            output_sine_phase += 2 * M_PI * 2000. / 96000.;

            if (output_sine_phase > 2 * M_PI) {
                output_sine_phase -= 2 * M_PI;
            }

            data.voice_high_res.at(test_counter)->at(i) = 0.9 * std::sin(output_sine_phase);
            float lin = data.input_l->at(i);
            float rin = data.input_r->at(i);
            float rin_lp = _lowpass_r.process(rin);
            rin = rin_lp - _highpass_r.process(rin);
            float lin_lp = _lowpass_l.process(lin);
            lin = lin_lp - _highpass_l.process(lin);
            max = lin > max ? lin : max;
            min = lin < min ? lin : min;
            float l = lin;
            float r = rin;
            level_l += l * l;
            level_r += r * r;
        }

        if (test_counter == NUM_VOICES) {
            return;
        }

        _controller.at(test_counter).setMainOutVca(vca_state_l, vca_state_r);

        for (auto &controller : _controller) {
            controller.run();
        }

        if (buffer_count % 400 == 0 && buffer_count > 0) {
            test_counter_2++;
            if (test_counter_2 > 0) {
                level_l = std::sqrt(level_l);
                level_r = std::sqrt(level_r);
                _audio_data_out.push_back((vca_state_l));
                _audio_data_out.push_back((level_l));
                _audio_data_out.push_back((level_r));
                old_delta_l = -old_level_l + level_l;
                old_delta_r = -old_level_r + level_r;
                printf("\n v %d levels: %f %f %f %f %f %f ", test_counter, vca_state_l, level_l, old_delta_l, vca_state_r, level_r, old_delta_r);
                if (buffer_count % 400 == 0) {

                    if (min_l > level_l) {
                        min_l = level_l;
                        min_val_l = vca_state_l;
                        printf("\n new min %f %f ", min_l, min_val_l);
                    }
                    queue_l.at(queue_i) = vca_state_l;
                    queue_r.at(queue_i) = vca_state_r;
                    queue_i += 1;
                    queue_i = queue_i == 20 ? 0 : queue_i;
                    float oldest_val_l = queue_l.at(queue_i);
                    float oldest_val_r = queue_r.at(queue_i);
                    printf("\ncheck %f %f", oldest_val_l, vca_state_l);
                    if (abs(oldest_val_l - vca_state_l) < abs(2 * step_l + 0.0000001)) {
                        float sum = 0;
                        for (auto val : queue_l) {
                            sum += val;
                            printf("cal:  %f ", val);
                        }
                        sum = sum / 20.;
                        vca_output_l = sum;
                        printf("\n min: %f %f ", min_l, min_val_l);
                        printf("\nfinal val L: %f %d", sum, test_counter);
                        done_l = true;
                    }
                    if (abs(oldest_val_r - vca_state_r) < abs(2 * step_r + 0.0000001)) {
                        _first_pass_r = true;
                        if (!_first_pass_r) {
                            printf("\n first pass");
                        } else {
                            float sum = 0;
                            for (auto val : queue_r) {
                                sum += val;
                            }
                            sum = sum / 20.;
                            vca_output_r = sum;
                            printf("\nfinal val R: %f %d", sum, test_counter);
                            done_r = true;
                        }
                    }

                    old_level_l = level_l;
                    old_level_r = level_r;
                }
            }
            level_l = 0;
            level_r = 0;
        }

        if (buffer_count % 400 == 0) {
            test_counter_2 = 0;
            if (old_delta_l > 0) {
                step_l = -step_l;
            } else {
            }
            if (old_delta_r > 0) {
                step_r = -step_r;
            } else {
            }
            float oldest_val_r = queue_r.at(queue_i);
            if (abs(oldest_val_r - vca_state_r) < abs(3 * step_l)) {
                if (level_l < 0.0009) {
                }
            }
            vca_state_l += step_l;
            vca_state_r += step_r;
        }
        if (done_l && done_r) {
            _analog_voices.at(test_counter).setMainRVcaOffset(vca_output_r);
            _analog_voices.at(test_counter).setMainLVcaOffset(vca_output_l);
            done_l = false;
            done_r = false;
            buffer_count = -1;
            printf("\nfinished voice");
            _out.write(reinterpret_cast<const char *>(_audio_data_out.data()), sizeof(float) * (_audio_data_out.size()));
            _out.close();
            test_counter++;
            printf("\ndone");
        }
        buffer_count++;
    }

    void writeOutVector(std::vector<float> &data, std::string file_name) {

        std::ofstream out;
        out.open(file_name, std::ios::out | std::ios::trunc | std::ios::binary);
        out.write(reinterpret_cast<const char *>(data.data()), sizeof(float) * (data.size()));
        out.close();
        printf("\n write file %s %d\n", file_name.c_str(), (int)data.size());
    }

    void startVcaCal(bool start) {
        printf("\n setup VCA cal ");
        _test_complete = false;
        _voice_counter = 0;
        _testing = testRunning::MainVca;
        _audio_data_out.clear();
        _stimulus_out.clear();
    }

    void runVCA2(ProcessAudioData &data) {
        constexpr float INIT_TEST_VAL = -0.02;
        constexpr float FINISH_TEST_VAL = 0.02;
        constexpr float SEARCH_WIDTH_2 = 0.008;
        constexpr float TEST_RUNTIME = 3;
        constexpr float CLICK_WAIT_TIME = .2;
        constexpr float DC_BLOCK = 0.001;
        if (_test_complete) {
            return;
        }

        switch (vca_state) {
        case VCA_State::BEGIN: {
            printf("\n start vca test %d ", _voice_counter);
            test_counter = 0;
            _stimulus_out.clear();
            _rms_1.clear();
            _rms_2.clear();
            vca_state = VCA_State::CLICKWAIT;

            // set all voices to off so everything execpt the target voice will be muted
            for (uint i = 0; i < NUM_VOICES; i++) {
                _controller.at(i).setOff();
            }
            _analog_voices.at(_voice_counter).setMainLVcaOffset(0);
            _analog_voices.at(_voice_counter).setMainRVcaOffset(0);
        } break;
        case VCA_State::CLICKWAIT: {
            bool finish = (test_counter > (int)(CLICK_WAIT_TIME * BUFFER_RATE));
            if (finish) {
                vca_state = VCA_State::RUNNING;
                test_counter = 0;
            }
            //
            float val = 0;
            _controller.at(_voice_counter).setMainOutVca2(INIT_TEST_VAL, INIT_TEST_VAL);
            for (auto &controller : _controller) {
                controller.run();
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }
            float old_dc_cal_val = dc_cal_value;

            // set the test output wave
            float offset = _analog_voices.at(_voice_counter).getModelData().voice_offset;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = (((float)s) / ((float)BUFFER_SIZE) - 0.5) * 0.5 + offset;
            }

            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_r->at(s);
                float s_r = data.input_l->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
            }
            test_counter++;
        } break;
        case VCA_State::RUNNING: {

            bool finish = (test_counter > (int)(TEST_RUNTIME * BUFFER_RATE));
            if (finish) {
                vca_state = VCA_State::CLICKWAIT_2;
                test_counter = 0;

                float min_rms_l = 1;
                float min_rms_r = 1;
                int index_l, index_r;
                for (int i = 0; i < _rms_1.size(); i++) {
                    if (_rms_1.at(i) < min_rms_l) {
                        min_rms_l = _rms_1.at(i);
                        index_l = i;
                    }
                    if (_rms_2.at(i) < min_rms_r) {
                        min_rms_r = _rms_2.at(i);
                        index_r = i;
                    }
                }
                min_l = _stimulus_out.at(index_l);
                min_r = _stimulus_out.at(index_r);

                printf("\n finish stage 1 %d %f %f %f %f", _voice_counter, min_l, _stimulus_out.at(index_l), min_r, _stimulus_out.at(index_r));
                // write out analysis
                std::stringstream file;
                file << "/udata/delia/tuning/"
                     << "voice_" << _voice_counter << "_vca_coarse_rms_l.dat";
                std::cout << file.str();
                writeOutVector(_rms_1, file.str());
                // clear file
                file.str("");
                file << "/udata/delia/tuning/"
                     << "voice_" << _voice_counter << "_vca_coarse_rms_r.dat";
                std::cout << file.str();
                writeOutVector(_rms_2, file.str());
                file.str("");
                file << "/udata/delia/tuning/"
                     << "voice_" << _voice_counter << "_vca_coarse_stim_l.dat";
                std::cout << file.str();
                writeOutVector(_stimulus_out, file.str());
                file.str("");
                file << "/udata/delia/tuning/"
                     << "voice_" << _voice_counter << "_vca_coarse_stim_r.dat";
                std::cout << file.str();
                writeOutVector(_stimulus_out, file.str());

                _rms_1.clear();
                _rms_2.clear();
                _stim_l.clear();
                _stim_r.clear();
                break;
            }
            float cal_stim = INIT_TEST_VAL + (FINISH_TEST_VAL - INIT_TEST_VAL) * ((float)test_counter / (float)(TEST_RUNTIME * BUFFER_RATE));
            _controller.at(_voice_counter).setMainOutVca2(cal_stim, cal_stim);
            for (auto &controller : _controller) {
                controller.run();
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }

            // set the test output wave
            float offset = _analog_voices.at(_voice_counter).getModelData().voice_offset;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = (((float)s) / ((float)BUFFER_SIZE) - 0.5) * 0.5 + offset;
            }

            // get RMS level for current buffer
            float rms_l = 0;
            float rms_r = 0;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_r->at(s);
                float s_r = data.input_l->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
                rms_l += s_l * s_l;
                rms_r += s_r * s_r;
            }
            rms_l = std::sqrt(rms_l / BUFFER_SIZE);
            rms_r = std::sqrt(rms_r / BUFFER_SIZE);

            _stimulus_out.push_back(cal_stim);
            _rms_1.push_back(rms_l);
            _rms_2.push_back(rms_r);
            test_counter++;
        } break;
        case VCA_State::CLICKWAIT_2: {
            bool finish = (test_counter > (int)(CLICK_WAIT_TIME * BUFFER_RATE));
            if (finish) {
                vca_state = VCA_State::RUNNING_2;
                test_counter = 0;
                break;
            }
            float cal_stim_l = min_l - (SEARCH_WIDTH_2 / 2.f);
            float cal_stim_r = min_r - (SEARCH_WIDTH_2 / 2.f);
            _controller.at(_voice_counter).setMainOutVca2(cal_stim_l, cal_stim_r);
            for (auto &controller : _controller) {
                controller.run();
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }
            float old_dc_cal_val = dc_cal_value;

            // set the test output wave
            float offset = _analog_voices.at(_voice_counter).getModelData().voice_offset;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = (((float)s) / ((float)BUFFER_SIZE) - 0.5) * 0.5 + offset;
            }

            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_r->at(s);
                float s_r = data.input_l->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
            }
            test_counter++;
        } break;
        case VCA_State::RUNNING_2: {

            bool finish = (test_counter > (int)(TEST_RUNTIME * BUFFER_RATE));
            if (finish) {
                vca_state = VCA_State::FINISH;
                test_counter = 0;
                break;
            }
            float cal_stim_l = (min_l - SEARCH_WIDTH_2 / 2.f) + SEARCH_WIDTH_2 * ((float)test_counter / (float)(TEST_RUNTIME * BUFFER_RATE));
            float cal_stim_r = (min_r - SEARCH_WIDTH_2 / 2.f) + SEARCH_WIDTH_2 * ((float)test_counter / (float)(TEST_RUNTIME * BUFFER_RATE));
            _controller.at(_voice_counter).setMainOutVca2(cal_stim_l, cal_stim_r);
            for (auto &controller : _controller) {
                controller.run();
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }

            // set the test output wave
            float offset = _analog_voices.at(_voice_counter).getModelData().voice_offset;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = (((float)s) / ((float)BUFFER_SIZE) - 0.5) * 0.5 + offset;
            }

            // get RMS level for current buffer
            float rms_l = 0;
            float rms_r = 0;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_r->at(s);
                float s_r = data.input_l->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
                rms_l += s_l * s_l;
                rms_r += s_r * s_r;
            }
            rms_l = std::sqrt(rms_l / BUFFER_SIZE);
            rms_r = std::sqrt(rms_r / BUFFER_SIZE);

            _stim_l.push_back(cal_stim_l);
            _stim_r.push_back(cal_stim_r);
            _rms_1.push_back(rms_l);
            _rms_2.push_back(rms_r);
            test_counter++;
        } break;
        case VCA_State::FINISH: {

            constexpr float filter_coeff = .05;
            // filter audio data
            {
                auto &filter_data = _rms_1;
                float filter = filter_data.at(0);
                for (int i = 0; i < filter_data.size(); i++) {
                    filter += (filter_data.at(i) - filter) * filter_coeff;
                    filter_data.at(i) = filter;
                }
                // filter backwards to remove phase offset
                filter = filter_data.at(filter_data.size() - 1);
                for (int i = filter_data.size() - 1; i > -1; i--) {
                    filter += (filter_data.at(i) - filter) * filter_coeff;
                    filter_data.at(i) = filter;
                }
            }
            {
                auto &filter_data = _rms_2;
                float filter = filter_data.at(0);
                for (int i = 0; i < filter_data.size(); i++) {
                    filter += (filter_data.at(i) - filter) * filter_coeff;
                    filter_data.at(i) = filter;
                }
                // filter backwards to remove phase offset
                filter = filter_data.at(filter_data.size() - 1);
                for (int i = filter_data.size() - 1; i > -1; i--) {
                    filter += (filter_data.at(i) - filter) * filter_coeff;
                    filter_data.at(i) = filter;
                }
            }

            // write out analysis
            std::stringstream file;
            file << "/udata/delia/tuning/"
                 << "voice_" << _voice_counter << "_vca_rms_l.dat";
            std::cout << file.str();
            writeOutVector(_rms_1, file.str());
            // clear file
            file.str("");
            file << "/udata/delia/tuning/"
                 << "voice_" << _voice_counter << "_vca_rms_r.dat";
            std::cout << file.str();
            writeOutVector(_rms_2, file.str());
            file.str("");
            file << "/udata/delia/tuning/"
                 << "voice_" << _voice_counter << "_vca_stim_l.dat";
            std::cout << file.str();
            writeOutVector(_stim_l, file.str());
            file.str("");
            file << "/udata/delia/tuning/"
                 << "voice_" << _voice_counter << "_vca_stim_r.dat";
            std::cout << file.str();
            writeOutVector(_stim_r, file.str());

            min_l = 1;
            min_r = 1;
            int index_l, index_r;
            for (int i = 0; i < _rms_1.size(); i++) {
                if (_rms_1.at(i) < min_l) {
                    min_l = _rms_1.at(i);
                    index_l = i;
                }
                if (_rms_2.at(i) < min_r) {
                    min_r = _rms_2.at(i);
                    index_r = i;
                }
            }
            // reset vars
            printf("\n finish vca test %d min: %f %f  %f %f", _voice_counter, min_l, _stim_l.at(index_l), min_r, _stim_r.at(index_r));
            _analog_voices.at(_voice_counter).setMainLVcaOffset(_stim_l.at(index_l));
            _analog_voices.at(_voice_counter).setMainRVcaOffset(_stim_r.at(index_r));
            vca_state = VCA_State::BEGIN;
            _voice_counter++;
            if (_voice_counter >= NUM_VOICES) {
                for (auto &item : _analog_voices) {
                    item.saveCalibration();
                }
                _test_complete = true;
                std::ofstream outfile;
                outfile.open(CAL_FILE_DONE_PATH, std::ios_base::app); // append instead of overwrite
                outfile << " ";
                printf("\nfinished vca cal");
            }
        }

        default:
            break;
        }
    }

    void
    startDcCal(bool start) {

        printf("\n setup DC cal ");
        _test_complete = false;
        _voice_counter = 0;
        _testing = DC;
        _audio_data_out.clear();
        _stimulus_out.clear();
        dc_overdrive = false;
    }

    void startDcCalOverdrive(bool start) {
        printf("\n setup DC cal ");
        _test_complete = false;
        _voice_counter = 0;
        _testing = DC;
        _audio_data_out.clear();
        _stimulus_out.clear();
        dc_overdrive = true;
    }

    int sine_phase = 0;

    void runDC(ProcessAudioData &data) {
        constexpr float DC_INIT_TEST_VAL = -0.2;
        constexpr float DC_FINISH_TEST_VAL = 0.2;
        constexpr float DC_2ND_STAGE_TEST_WIDTH = 0.005;
        constexpr float DC_TEST_RUNTIME = 2;
        constexpr float CLICK_WAIT_TIME = .8;
        constexpr float DC_BLOCK = 0.0003;
        if (_test_complete) {
            return;
        }

        switch (dc_state) {
        case DC_State::BEGIN: {
            printf("\n start dc test %d ", _voice_counter);
            test_counter = 0;
            _audio_data_out.clear();
            _stimulus_out.clear();
            dc_cal_value = DC_INIT_TEST_VAL;
            dc_state = DC_State::CLICKWAIT;
            sine_phase = 0;

            // set all voices to off so everything execpt the target voice will be muted
            for (uint i = 0; i < NUM_VOICES; i++) {
                _controller.at(i).setOff();
            }
        } break;
        case DC_State::CLICKWAIT: {

            bool finish = (test_counter > (int)(CLICK_WAIT_TIME * BUFFER_RATE));
            if (finish) {
                dc_state = DC_State::RUNNING;
                test_counter = 0;
            }

            // override the VCA output levels to the target voice
            float val = 0;
            _controller.at(_voice_counter).setDCtest(dc_overdrive);
            _controller.at(_voice_counter)._mix_l = val;
            _controller.at(_voice_counter)._mix_r = val;
            for (auto &controller : _controller) {
                controller.run();
            }

            for (int s = 0; s < CV_BUFFER_SIZE; s++) {
                val = 0.f;
                if (s >= (CV_BUFFER_SIZE / 2))
                    val = 1.f;
                _controller.at(_voice_counter).voice_input.vca_l.at(s) = val;
                _controller.at(_voice_counter).voice_input.vca_r.at(s) = val;
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }
            float old_dc_cal_val = dc_cal_value;
            dc_cal_value = DC_INIT_TEST_VAL;

            // printf("\n click wait %f %d ", dc_cal_value, test_counter);
            //   set the DC level to test
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = dc_cal_value;
            }

            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_l->at(s);
                float s_r = data.input_r->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
            }
            test_counter++;
        } break;
        case DC_State::RUNNING: {
            bool finish = (test_counter > (int)(DC_TEST_RUNTIME * BUFFER_RATE));
            if (finish) {
                dc_state = DC_State::CLICKWAIT_2;
                test_counter = 0;

                // filter audio data
                float filter = _audio_data_out.at(0);
                constexpr float filter_coeff = 0.05;
                for (int i = 0; i < _audio_data_out.size(); i++) {
                    filter += (_audio_data_out.at(i) - filter) * filter_coeff;
                    _audio_data_out.at(i) = filter;
                }
                // filter backwards to remove phase offset
                filter = _audio_data_out.at(_audio_data_out.size() - 1);
                for (int i = _audio_data_out.size() - 1; i > -1; i--) {
                    filter += (_audio_data_out.at(i) - filter) * filter_coeff;
                    _audio_data_out.at(i) = filter;
                }
                std::stringstream file;
                file << "/udata/delia/tuning/"
                     << "voice_" << _voice_counter << "_dc_coarse_rms.dat";
                std::cout << file.str();
                writeOutVector(_audio_data_out, file.str());
                // clear file
                file.str("");
                file << "/udata/delia/tuning/"
                     << "voice_" << _voice_counter << "_dc_coarse_stim.dat";
                std::cout << file.str();
                writeOutVector(_stimulus_out, file.str());
                // clear file
                float min = 1;
                int index;
                for (int i = 0; i < _audio_data_out.size(); i++) {
                    if (_audio_data_out.at(i) < min) {
                        min = _audio_data_out.at(i);
                        index = i;
                    }
                }
                dc_first_stage_min = _stimulus_out.at(index);
                _stimulus_out.clear();
                _audio_data_out.clear();
            }

            // override the VCA output levels to the target voice
            float val = 0;
            _controller.at(_voice_counter).setDCtest(dc_overdrive);
            _controller.at(_voice_counter)._mix_l = val;
            _controller.at(_voice_counter)._mix_r = val;
            for (auto &controller : _controller) {
                controller.run();
            }

            for (int s = 0; s < CV_BUFFER_SIZE; s++) {
                val = 0.f;
                if (s >= (CV_BUFFER_SIZE / 2))
                    val = 1.f;
                _controller.at(_voice_counter).voice_input.vca_l.at(s) = val;
                _controller.at(_voice_counter).voice_input.vca_r.at(s) = val;
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }

            float old_dc_cal_val = dc_cal_value;
            dc_cal_value = DC_INIT_TEST_VAL + (DC_FINISH_TEST_VAL - DC_INIT_TEST_VAL) * ((float)test_counter / (float)(DC_TEST_RUNTIME * BUFFER_RATE));
            // printf("\n val %f %d ", dc_cal_value, test_counter);
            //   set the DC level to test
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = dc_cal_value;
            }

            // get RMS level for current buffer
            float rms = 0;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_l->at(s);
                float s_r = data.input_r->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
                rms += s_l * s_l + s_r * s_r;
            }
            rms = std::sqrt(rms / BUFFER_SIZE * 2);

            // printf(" rms %f", rms);
            _stimulus_out.push_back(old_dc_cal_val);
            _audio_data_out.push_back(rms);
            test_counter++;
        } break;
        case DC_State::CLICKWAIT_2: {

            bool finish = (test_counter > (int)(CLICK_WAIT_TIME * BUFFER_RATE));
            if (finish) {
                dc_state = DC_State::RUNNING_2;
                _stimulus_out.clear();
                _audio_data_out.clear();
                test_counter = 0;
                dc_cal_value = dc_first_stage_min + -DC_2ND_STAGE_TEST_WIDTH / 2 + (DC_2ND_STAGE_TEST_WIDTH) * ((float)0 / (float)(DC_TEST_RUNTIME * BUFFER_RATE));
                break;
            }

            // override the VCA output levels to the target voice
            float val = 0;
            _controller.at(_voice_counter).setDCtest(dc_overdrive);
            _controller.at(_voice_counter)._mix_l = val;
            _controller.at(_voice_counter)._mix_r = val;
            for (auto &controller : _controller) {
                controller.run();
            }

            for (int s = 0; s < CV_BUFFER_SIZE; s++) {
                val = 0.f;
                if (s > (CV_BUFFER_SIZE / 2))
                    val = 1.f;
                _controller.at(_voice_counter).voice_input.vca_l.at(s) = val;
                _controller.at(_voice_counter).voice_input.vca_r.at(s) = val;
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }
            float old_dc_cal_val = dc_cal_value;
            dc_cal_value = dc_first_stage_min + -DC_2ND_STAGE_TEST_WIDTH / 2 + (DC_2ND_STAGE_TEST_WIDTH) * ((float)0 / (float)(DC_TEST_RUNTIME * BUFFER_RATE));

            // printf("\n click wait 2 %f %d ", dc_cal_value, test_counter);
            //    set the DC level to test
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = dc_cal_value;
            }

            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_l->at(s);
                float s_r = data.input_r->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
            }
            test_counter++;
        } break;
        case DC_State::RUNNING_2: {
            bool finish = (test_counter > (int)(DC_TEST_RUNTIME * BUFFER_RATE));
            if (finish) {
                dc_state = DC_State::FINISH;
            }

            // override the VCA output levels to the target voice
            float val = 0;
            _controller.at(_voice_counter).setDCtest(dc_overdrive);
            _controller.at(_voice_counter)._mix_l = val;
            _controller.at(_voice_counter)._mix_r = val;
            for (auto &controller : _controller) {
                controller.run();
            }

            for (int s = 0; s < CV_BUFFER_SIZE; s++) {
                val = 0.f;
                if (s > CV_BUFFER_SIZE / 2)
                    val = 1.f;
                _controller.at(_voice_counter).voice_input.vca_l.at(s) = val;
                _controller.at(_voice_counter).voice_input.vca_r.at(s) = val;
            }

            // init all voice output data to 0
            for (int v = 0; v < NUM_VOICES; v++) {
                for (int s = 0; s < BUFFER_SIZE; s++) {
                    data.voice_high_res.at(v)->at(s) = 0;
                }
            }
            float old_dc_cal_val = dc_cal_value;
            dc_cal_value = dc_first_stage_min + -DC_2ND_STAGE_TEST_WIDTH / 2 + (DC_2ND_STAGE_TEST_WIDTH) * ((float)test_counter / (float)(DC_TEST_RUNTIME * BUFFER_RATE));

            // printf("\n val %f %d ", dc_cal_value, test_counter);
            //    set the DC level to test
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(_voice_counter)->at(s) = dc_cal_value;
            }

            // get RMS level for current buffer
            float rms = 0;
            for (int s = 0; s < BUFFER_SIZE; s++) {
                float s_l = data.input_l->at(s);
                float s_r = data.input_r->at(s);
                s_l -= dc_l;
                dc_l += s_l * DC_BLOCK;
                s_r -= dc_r;
                dc_r += s_r * DC_BLOCK;
                rms += s_l * s_l + s_r * s_r;
            }
            rms = std::sqrt(rms / BUFFER_SIZE * 2);
            _stimulus_out.push_back(old_dc_cal_val);
            _audio_data_out.push_back(rms);
            test_counter++;
        } break;
        case DC_State::FINISH: {
            // filter audio data
            float filter = _audio_data_out.at(0);
            constexpr float filter_coeff = 1;
            for (int i = 0; i < _audio_data_out.size(); i++) {
                filter += (_audio_data_out.at(i) - filter) * filter_coeff;
                _audio_data_out.at(i) = filter;
            }
            // filter backwards to remove phase offset
            filter = _audio_data_out.at(_audio_data_out.size() - 1);
            for (int i = _audio_data_out.size() - 1; i > -1; i--) {
                filter += (_audio_data_out.at(i) - filter) * filter_coeff;
                _audio_data_out.at(i) = filter;
            }
            // write out analysis

            std::stringstream file;
            file << "/udata/delia/tuning/"
                 << "voice_" << _voice_counter << "_dc_rms.dat";
            std::cout << file.str();
            writeOutVector(_audio_data_out, file.str());
            // clear file
            file.str("");
            file << "/udata/delia/tuning/"
                 << "voice_" << _voice_counter << "_dc_stim.dat";
            std::cout << file.str();
            writeOutVector(_stimulus_out, file.str());
            // clear file
            float min = 1;
            int index;

            for (int i = 0; i < _audio_data_out.size(); i++) {
                if (_audio_data_out.at(i) < min) {
                    min = _audio_data_out.at(i);
                    index = i;
                }
            }
            // reset vars
            printf("\n finish dc test %d min: %f %f", _voice_counter, min, _stimulus_out.at(index));
            if (dc_overdrive) {
                _analog_voices.at(_voice_counter).setVoiceOffsetOd(_stimulus_out.at(index));
            } else {
                _analog_voices.at(_voice_counter).setVoiceOffset(_stimulus_out.at(index));
            }

            dc_state = DC_State::BEGIN;
            _voice_counter++;
            if (_voice_counter >= NUM_VOICES) {
                for (auto &item : _analog_voices) {
                    item.saveCalibration();
                }
                _test_complete = true;
                std::ofstream outfile;
                outfile.open(CAL_FILE_DONE_PATH, std::ios_base::app); // append instead of overwrite
                outfile << " ";
                printf("\nfinished dc cal");
            }
        }
        case DC_State::TEST: {

            bool finish = (test_counter > (int)(0.5 * BUFFER_RATE));
        }

        default:
            break;
        }
    }

    static constexpr int test_len = 36;
    static constexpr float stim_values[test_len] =
        {0.000, 0.018, 0.036, 0.055, 0.073, 0.091, 0.109, 0.127, 0.145, 0.164, 0.182, 0.200, 0.218, 0.236, 0.255, 0.273, 0.291, 0.309, 0.327, 0.345, 0.364, 0.382, 0.400, 0.418, 0.436, 0.455, 0.473, 0.491, 0.509, 0.527, 0.545, 0.564, 0.582, 0.600, 0.700, 0.8};
    std::array<float, test_len> stim_times =
        {500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 500.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.000, 250.0};

    void runFilter(ProcessAudioData &data) {
        if (_test_complete) {
            return;
        }
        if (counter >= test_len) {

            printf("\n Test Complete %d\n", test_counter);
            FILE *fp = popen("vcgencmd pmicrd 1c", "r");
            char buff[512];
            auto res = fgets(buff, sizeof buff, fp);
            std::string temp_str = std::string(res);
            temp_str = temp_str.substr(7, temp_str.length());
            float temp = (float)(int)std::stoi(temp_str.c_str(), 0, 16);
            printf("\n temp: %f", temp);
            _stimulus_out.push_back(temp);
            printf("\nsize:  %d", (int)_audio_data_out.size());
            _out_2.write(reinterpret_cast<const char *>(_audio_data_out.data()), sizeof(float) * (_audio_data_out.size() - 1));
            printf("\nsize:  %d", (int)_stimulus_out.size());
            _out.write(reinterpret_cast<const char *>(_stimulus_out.data()), sizeof(float) * (_stimulus_out.size()));
            _out.close();
            _out_2.close();
            _call_counter = -1;

            counter = 1;
            _audio_data_out.clear();
            _stimulus_out.clear();
            test_counter++;
            if (test_counter == NUM_VOICES) {
                if (_test_complete) {
                    return;
                }

                _test_complete = true;
                std::ofstream outfile;
                outfile.open(CAL_FILE_DONE_PATH, std::ios_base::app); // append instead of overwrite
                outfile << " ";
                return;
            }
            counter = 0;
            setupFilter(test_counter);
            printf("\n filter setup voice %d", test_counter);
            stim = stim_values[counter];
            _call_counter = 0;
        }

        if (_call_counter >= stim_times.at(counter)) {
            counter++;
            stim = stim_values[counter];
            _call_counter = 0;
            printf("\n set filter test %f %d %d\n", stim, _call_counter, counter);
        }

        _controller.at(test_counter).setFilterTuneTest(stim);
        for (auto &controller : _controller) {
            controller.run();
        }
        for (int v = 0; v < NUM_VOICES; v++) {
            for (int s = 0; s < BUFFER_SIZE; s++) {
                data.voice_high_res.at(v)->at(s) = 0;
            }
        }
        if (_call_counter % 100 == 0) {
        }

        if (_call_counter > 0 && _call_counter < stim_times[counter]) {
            for (int i = 0; i < BUFFER_SIZE; i++) {
                _audio_data_out.push_back(data.input_r->at(i));
                data.output_l->at(i) = 0;
                data.output_r->at(i) = 0;
            }
            _stimulus_out.push_back(_analog_voices.at(test_counter).getFilterOut().at(0));
            _stimulus_out.push_back(_analog_voices.at(test_counter).getFilterOut().at(1));
            _stimulus_out.push_back(stim);
        }

        _call_counter++;
    }

  public:
    MoniqueCalSequencer() {

        _audio_data_out.reserve(SAMPLE_RATE * 2);
        _stimulus_out.reserve(SAMPLE_RATE * 2);
    };

    bool isTestRunning() {
        return !_test_complete;
    };

    ~MoniqueCalSequencer(){

    };

    /**
     * @brief Init the Main LR vca cal. This works in a similar way to the mix VCA cal. A tone is run though the voice digital output into the filer. the LR vca's are adjusted to minimise the RMS level at the output. Each voice is set one after the other
     *
     * @param run
     */
    void startMainVcaCal(bool run) {
        if (run) {
            printf("\n setup main vca cal");
            _testing = MainVca;
            _test_complete = false;
            vca_start = -0.001;
            vca_state_l = vca_start;
            vca_state_r = vca_start;
            level_l = 0;
            level_r = 0;
            dc_l = 0;
            dc_r = 0;
            wo = false;
            step = .0005;
            step_l = step;
            step_r = step;
            old_level_l = 100;
            old_level_r = 100;
            old_delta_l = 0;
            old_delta_r = 0;
            buffer_count = 0;
            min = 1000;
            max = -1000;
            sample_counter = 0;
            output_sine_phase = 0;
            queue_i = 0;
            done_l = false;
        }
    }

    /**
     * @brief The filter testing works by turning the filter resonance up, and saving audio samples of a series of test points. this data is analysed in a python script which calculates the res freq of each testpoint, calculating the filter cal values that best match. the Osc 'temp' reading is also saved as part of this data, so it can be used to set the relationship between the osc 'temp' and the filter cal. (this is because we open loop correct the temp of the filter, by using the temp of the osc as a stand in value)
     *
     * @param run
     */
    void startFilterCal(bool run) {
        if (run) {
            printf("\n setup filter testing");

            _test_complete = false;
            _testing = Filter;
            _call_counter = 0;
            done_l = false;
            done_r = false;
            buffer_count = 0;
            stim = 0.95;
            _audio_data_out.clear();
            _stimulus_out.clear();
            test_counter = -1;
            counter = 100;
        }
    };

    void runCal(uint voice_num);

    float stim = 0.95;
    bool _first_pass_l = false;
    bool _first_pass_r = false;
    float min_l = 1;
    float min_r = 1;
    float min_val_l, min_val_r;
    int l = 20;
    int counter = 1;
    int _test_counter = 0;
    int test_counter = 0;
    int test_counter_2 = 0;
    uint _voice_counter = 0;
    bool _test_complete = false;

    enum testRunning {
        NoTest = 0,
        MainVca,
        Filter,
        DC
    };

    testRunning _testing = NoTest;
    float vca_start = -0.001;
    float vca_state_l = vca_start;
    float vca_state_r = vca_start;
    float vca_output_l;
    float vca_output_r;
    float level_l = 0;
    float level_r = 0;
    float dc_l = 0;
    float dc_first_stage_min = 0;
    float dc_r = 0;
    float dc_cal_value = 0;
    bool dc_overdrive = false;
    bool wo = false;
    float step = .005;
    float step_l = step;
    float step_r = step;
    int buf_ave_cnt = 1;
    int buf_ave_cnt_store = 1;
    float old_level_l = 100;
    float old_level_r = 100;
    float old_delta_l = 0;
    float old_delta_r = 0;
    uint buffer_count = 0;
    float min = 1000;
    float max = -1000;
    int sample_counter = 0;
    float output_sine_phase = 0;
    uint queue_i = 0;
    std::array<float, 20> queue_l;
    std::array<float, 20> queue_r;
    bool done_l = false;
    bool done_r = false;
    float _tri_0, _tri_1, _sqr_0, _sqr_1, _xor;
    uint _current_vca = 0;

    void reloadCal() {
        for (auto &voice : _analog_voices) {
            voice.reloadCal();
        }
    }

    void run(ProcessAudioData &buffer_data) {
        switch (_testing) {
        case NoTest:
            break;
        case MainVca:
            // old (nina) vca cal
            // runMainVca(buffer_data);

            // new VCA cal routine
            runVCA2(buffer_data);
            break;

        case Filter:
            runFilter(buffer_data);
            break;
        case DC:
            runDC(buffer_data);
        default:
            break;
        }
        for (int v = 0; v < NUM_VOICES; v++)

        {
            _analog_voices.at(v).run(*buffer_data.voice_cv.at(v));
        }
    };
};

} // namespace Monique