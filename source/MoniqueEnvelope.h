#pragma once

#include "MoniqueCommon.h"
#include "SynthMath.h"
#include "common.h"

namespace Monique {

enum class AdsrState { ATT,
    DEC,
    REL
};

constexpr float SLOW_LINEAR_THRESH = 0.15;
constexpr float SLOW_MULTIPLIER = 0.01;

class GenAdsrEnvelope {
  public:
    GenAdsrEnvelope(float &velocity, float &vel_sense, bool &reset, bool &drone, bool &voice_gate_on, bool &voice_trigger) :
        _velocity(velocity),
        _velocity_sense(vel_sense),
        _reset(reset),
        _drone(drone),
        _voice_gate_on(voice_gate_on),
        _voice_trigger(voice_trigger) {}

    ~GenAdsrEnvelope() {}

    /**
     * @brief set if the envelope should be looping or not TODO: [NINA-267] implement envelope looping
     *
     * @param should_loop enable for loop on
     */
    void setloop(bool should_loop);

    /**
     * @brief force the envelope to zero
     *
     */
    void forceReset() {
        _output_pregain = 0;
        _env_signal = 0;
    }

    void setSlowMode(bool slow) {
        _slow_mode = slow;
    }

    /**
     * @brief generate 1 sample of the env output. run at the CV rate
     *
     */
    void run() {

        // each state has an exp decay or release. only att->rel transition is triggered here so its more efficient
        float old = _env_signal;
        switch (_env_state) {
        case AdsrState::ATT: {
            _env_signal += _att_coeff * (att_asymtote - _env_signal);
            if (_env_signal > att_trigger) {
                _env_signal = att_trigger;
                _env_state = AdsrState::DEC;
            }
            _output_pregain = _env_signal;
            break;
        }
        case AdsrState::DEC: {
            _env_signal += _dec_coeff * (-_env_signal);
            _output_pregain = _env_signal + (1 - _env_signal) * _sus_level;
            break;
        }
        case AdsrState::REL: {
            _env_signal += _rel_coeff * (idle_level - _env_signal);
            _output_pregain = _env_signal;
            break;
        }
        default:
            break;
        }
        if (_drone) {
            _output_pregain = _sus_level;
        }
        _output = _output_pregain * _output_gain;
        if (!floatIsValid(_env_signal)) {
            printf("\n adsr  %f %f %f %f %f ", _output, _env_signal, _rel_input, _output_gain, old);
            // assert(false);
        }
    }

    void reCalculate() {
        // we remove these state transitions out of the run function so that they arn't called for every sample since it seems to have little impact on sound
        switch (_env_state) {
        case AdsrState::ATT: {
            if (!_voice_gate_on) {
                _env_state = AdsrState::REL;
            }
            break;
        }
        case AdsrState::DEC: {
            if (!_voice_gate_on) {
                _env_state = AdsrState::REL;
                _env_signal = _env_signal + (1 - _env_signal) * _sus_level;
            }
            if (_voice_trigger) {
                _env_state = AdsrState::ATT;
                if (_reset) {
                    _env_signal = 0;
                    // set the env signal to current output level to avoid clicks
                } else {
                    _env_signal = _env_signal + (1 - _env_signal) * _sus_level;
                }
            }
            break;
        }
        case AdsrState::REL: {
            if (_voice_gate_on) {
                _env_state = AdsrState::ATT;
            }
            if (_voice_trigger) {
                _env_state = AdsrState::ATT;
                if (_reset) {
                    _env_signal = 0;
                }
            }
            break;
        }
        }
        // clip each input signal so we dont have any numerical issues
        _att_input = unipolarClip(_att_input);
        _dec_input = unipolarClip(_dec_input);
        _rel_input = unipolarClip(_rel_input);

        // based off NINA time coeff scaling, math is refactored based off the following equations
        // NINA knob to coeff scaling: y=1-2^-1.4/(CV_Rate*5x^3)
        // try following in MS Mathmatics
        //  solve({y=1-2^-1.4/(u*5x^3),0.5=(1-y)^t,s=t/u}, {y,t,s})
        if (!_slow_mode) {
            _att_coeff = (1 - fastpow2(-(1.4 / ((float)CV_SAMPLE_RATE * 5. * _att_input * _att_input * _att_input))));
            _dec_coeff = (1 - fastpow2(-(1.4 / ((float)CV_SAMPLE_RATE * 5. * _dec_input * _dec_input * _dec_input))));
            _rel_coeff = (1 - fastpow2(-(1.4 / ((float)CV_SAMPLE_RATE * 5. * _rel_input * _rel_input * _rel_input))));
        } else {
            // slow mode scaling is different, equation is scaled to be faster, also has linear amount added towards 0 to allow instantanious

            const float lin_offset_att = _att_input < SLOW_LINEAR_THRESH ? (1.f - SLOW_MULTIPLIER - 0.001) * (SLOW_LINEAR_THRESH - _att_input) / SLOW_LINEAR_THRESH : 0.f;
            const float att_5 = _att_input * _att_input * _att_input * _att_input * _att_input;
            _att_coeff = (1 - fastpow2(-(1.4 / ((float)CV_SAMPLE_RATE * 5. * att_5)))) * SLOW_MULTIPLIER + lin_offset_att;

            const float lin_offset_dec = _dec_input < SLOW_LINEAR_THRESH ? (1.f - SLOW_MULTIPLIER - 0.001) * (SLOW_LINEAR_THRESH - _dec_input) / SLOW_LINEAR_THRESH : 0.f;
            const float dec_5 = _dec_input * _dec_input * _dec_input * _dec_input * _dec_input;
            _dec_coeff = (1 - fastpow2(-(1.4 / ((float)CV_SAMPLE_RATE * 5. * dec_5)))) * SLOW_MULTIPLIER + lin_offset_dec;

            const float lin_offset_rel = _rel_input < SLOW_LINEAR_THRESH ? (1.f - SLOW_MULTIPLIER - 0.001) * (SLOW_LINEAR_THRESH - _rel_input) / SLOW_LINEAR_THRESH : 0.f;
            const float rel_5 = _rel_input * _rel_input * _rel_input * _rel_input * _rel_input;
            _rel_coeff = (1 - fastpow2(-(1.4 / ((float)CV_SAMPLE_RATE * 5. * rel_5)))) * SLOW_MULTIPLIER + lin_offset_rel;
        }

        if (print) {
            // printf("\n adsr  %d %f %f %f %f ", _slow_mode, _output, _att_input, _att_coeff, _output_gain);
        }
        constexpr float b = 6.7;
        constexpr float a = 10;
        constexpr float c = 10;

        _sus_level = unipolarClip(_sus_input * _sus_input);

        // calculate velocity gain which is a blend of a static gain and velocity
        float vel_gain;
        if (_velocity_sense > 0.f) {
            vel_gain = 1.0f * (1 - _velocity_sense) + _velocity_sense * _velocity;
        } else {
            vel_gain = 1.f * (1 + _velocity_sense) + _velocity_sense * (_velocity - 1);
        }

        if (dump || print) {
            // printf("\n  vel %f %f %f ", vel_gain, _velocity_sense, _velocity);
        }

        // final gain

        _output_gain = cv_clip(_gain * vel_gain);

        if (!floatIsValid(_env_signal)) {
            //    printf("\n adsr  %f %f %f %f %f ", _output, _env_signal, _att_coeff, _output_gain, 0.f);
            // assert(false);
        }
    };

    bool inReleaseState() {
        return _env_state == AdsrState::REL;
    }

    float att_store = 0;

    float *getOutput() { return &_output; }

    float *getAttIn() { return &_att_input; }

    float *getDecIn() { return &_dec_input; }

    float *getSusIn() { return &_sus_input; }

    float *getRelIn() { return &_rel_input; }

    float *getGainIn() { return &_gain; }

    bool dump = false;
    bool print = false;

    static constexpr float env_control_clip_min = 0.f;
    static constexpr float env_control_clip_max = 1.f;
    static constexpr float att_trigger = 1.0f;
    static constexpr float att_asymtote = att_trigger / (1 - 1 / M_Ef32);
    static constexpr float idle_level = 0.f;
    static constexpr int time_const_translate = 19;
    const float time_const_offset = (exp2f32(time_const_mult * 0 - time_const_translate));
    float time_const_mult = 14;

    /**
     * @brief Coeff for env 6db rolloff filter.
     *
     */
    static constexpr float ENV_SMOOTH = 650 / (float)CV_SAMPLE_RATE;

    AdsrState _env_state = AdsrState::REL;
    float _smoothing = 1;
    float _env_signal = 0;
    float _att_coeff = 1;
    float _dec_coeff = 1;
    float _sus_level = 0.5;
    float _rel_coeff = 1;
    float _output = 0;
    float _gain = 1.0;
    float _output_gain = 0;
    float &_velocity;
    float &_velocity_sense;
    bool &_reset;
    bool &_drone;
    const bool &_voice_gate_on;
    const bool &_voice_trigger;

    float _att_input = 0;
    float _dec_input = 0;
    float _sus_input = 0;
    float _rel_input = 0;
    bool _slow_mode = false;
    float _output_pregain = 0;
};

} // namespace Monique
