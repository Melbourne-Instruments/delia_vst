
#pragma once

#include "MoniqueSlaveEnvelope.h"
#include "SynthMath.h"
#include "common.h"

namespace Monique {
constexpr float UNISON_GAIN_ADJ = exp10f(-2.5f / 20);

class NinaOutPanner {
  public:
    NinaOutPanner(float *left, float *right, float &overdrive_comp, float &max_mix_out_level, VoiceData &voice_data, float &unison_vel_1, float &unison_vel_2, GenAdsrEnvelope &vca_env, SlaveEnvelope &slave_vca, MoniqueSynthParameters &params) :
        _max_level(max_mix_out_level), _out_left(left), _out_right(right), _overdrive_comp(overdrive_comp),

        _voice_data(voice_data),
        _unison_vel_1(unison_vel_1),
        _unison_vel_2(unison_vel_2),
        _vca_env(vca_env),
        _slave_vca(slave_vca),
        _layer_params(params)

            {};

    ~NinaOutPanner(){};

    void run() {
        _spin_pan += _spin_inc;

        if (_voice_data.unison_osc_enable) {
            float vca = _vca_in;
            constexpr float QUIET_THRESH = 0.00001;
            const bool alloc_1 = _vca_env._output_pregain < QUIET_THRESH;
            const bool alloc_2 = _slave_vca._output_pregain < QUIET_THRESH;
            const float slave_lvl = _slave_vca._output_pregain;
            const float main_lvl = _vca_env._output_pregain;
            const bool main_vca_max = _vca_env._output > _slave_vca._output;
            float &unison_1 = _voice_data.unison_1_tmp;
            float &unison_2 = _voice_data.unison_2_tmp;
            if (main_vca_max) {
                unison_1 = _vca_env._output_gain * (1 - (_slave_vca._output * (1 - UNISON_GAIN_ADJ)));
                unison_2 = _slave_vca._output * UNISON_GAIN_ADJ;
                _vca_in = main_lvl * UNISON_GAIN_ADJ + (_vca_env._output * (1 - UNISON_GAIN_ADJ));
            } else {
                unison_1 = _vca_env._output * UNISON_GAIN_ADJ;
                unison_2 = _slave_vca._output_gain * (1 - (_vca_env._output * (1 - UNISON_GAIN_ADJ)));
                _vca_in = slave_lvl * UNISON_GAIN_ADJ + (_vca_env._output * (1 - UNISON_GAIN_ADJ));
            }

            if (!(alloc_1 && alloc_2)) {
                _vca_in = _vca_in * UNISON_GAIN_ADJ;
            }
        } else {
            _voice_data.unison_1_tmp = 1.f;
            _voice_data.unison_2_tmp = 0.f;
        }

        float vol = cv_clip((_vca_in)*_overdrive_comp);
        if (_layer_params.output_mode == OutputRouting::ONE_TWO) {
            _left_pan = _sin_pan * _filter_b + _left_pan * _filter_a;
            _right_pan = _cos_pan * _filter_b + _right_pan * _filter_a;
        } else {
            _left_pan = 0.707;
            _right_pan = 0.707;
        }
        *_out_left = _left_pan * vol;
        *_out_right = _right_pan * vol;

        // find the max output level for this buffer
        _max_level = vol > _max_level ? vol : _max_level;
        if (print) {
            // printf("\noutput: %f %d %d %f %f %d %d", vol, _vca_env._voice_gate_on, _slave_vca._voice_gate_on, _voice_data.unison_1_tmp, _voice_data.unison_2_tmp, main_vca_max, alloc_2);
        }
    }

    void reCalculate() {
        float spin_nor = fabsf32(_spin_rate);
        _spin_inc = (spin_nor * _spin_rate) * (max_spin_rate / (float)CV_SAMPLE_RATE);
        float spin_cut = -6 * (8 + 40 * spin_nor) / (float)BUFFER_RATE + 1.f;
        _filter_a = std::exp(-2 * M_PI * (4 + 40 * spin_nor) / CV_SAMPLE_RATE);
        _filter_b = 1 - _filter_a;

        // if the spin value is near 0, we decay the spin position to move the panning back to center smoothly
        if (spin_nor < SPIN_THRESH) {
            _spin_inc = 0;
            _spin_pan -= _spin_pan * spin_pan_decay;
        }
        _spin_pan = _spin_pan > M_PI ? _spin_pan - 2 * M_PI : _spin_pan;
        _spin_pan = _spin_pan < -M_PI ? _spin_pan + 2 * M_PI : _spin_pan;
        float scaled_pan_pos = (M_PI / 4) * _pan_position;
        const float pan = (_spin_pan - scaled_pan_pos) + M_PIf32 / (4);
        if (dump) {
        }
        _sin_pan = (std::sin(pan));
        _cos_pan = (std::cos(pan));

        // reset the max level value
        _max_level = 0;
    }

    float *getSpinIn() {
        return &_spin_rate;
    }

    void resetSpin() {
        _spin_pan = 0;
    }

    float *getPanIn() {
        return &_pan_position;
    }

    float *getVcaIn() {
        return &_vca_in;
    }

    bool dump = false;
    bool print = false;

    float &getMaxVcaVolume() {
        return _max_level;
    }

    static constexpr float PAN_SMOOTHING = std::exp(-2 * M_PI * 45.0 / BUFFER_RATE);
    static constexpr float SPIN_THRESH = 0.03;
    static constexpr float spin_pan_decay = 2.f / BUFFER_RATE;
    static constexpr float max_spin_rate = 4.5 * 2 * M_PI;
    /**
     * @brief phase increment calculation = 2^(log(2,2*pi/fs) + freq_exponential)
     * so we precalculate the log(2,2*pi/fs factor) TODO: [NINA-268] investigate making the phase factor a constexpr
     *
     */
    const float phase_factor = std::exp2f(2.0f * M_PIf32 / CV_SAMPLE_RATE);
    float _scaled_pan_pos = 0;
    float _sin_pan = 0;
    float _cos_pan = 0;
    float _left_pan = 0;
    float _right_pan = 0;
    float _spin_cut = 0;
    float _filter_a = 0;
    float _filter_b = 0;
    float &_max_level;
    float &_unison_vel_1;
    float &_unison_vel_2;
    VoiceData &_voice_data;
    GenAdsrEnvelope &_vca_env;
    SlaveEnvelope &_slave_vca;
    MoniqueSynthParameters &_layer_params;

    float _spin_inc = 0;
    float _spin_pan = 0;
    float _spin_rate = 0;
    float _pan_position = 0;
    float _volume = 1.0;
    float _vca_in = 0.5;

    float *const _out_left;
    float *const _out_right;
    float &_overdrive_comp;
};

} // namespace Monique