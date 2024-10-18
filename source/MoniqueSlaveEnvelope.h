#pragma once

#include "MoniqueCommon.h"
#include "MoniqueEnvelope.h"
#include "SynthMath.h"
#include "common.h"

namespace Monique {

class SlaveEnvelope {
  public:
    SlaveEnvelope(GenAdsrEnvelope &master_env, float &vel, bool &gate, bool &trig) :
        _velocity_sense(master_env._velocity_sense),
        _reset(master_env._reset),
        _drone(master_env._drone),
        _att_coeff(master_env._att_coeff),
        _rel_coeff(master_env._rel_coeff),
        _sus_level(master_env._sus_level),
        _dec_coeff(master_env._dec_coeff),
        _gain(master_env._gain),
        _voice_trigger(trig),
        _voice_gate_on(gate),
        _velocity(vel) {}

    ~SlaveEnvelope() {}

    void forceReset() {
        _output_pregain = 0;
        _env_signal = 0;
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
            printf("\n adsr  %f %f %f %f ", _output, _env_signal, _output_gain, old);
            // assert(false);
        }
    }

    void reCalculate() {

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
    };

    bool inReleaseState() {
        return _env_state == AdsrState::REL;
    }

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
    float _env_signal = 0;
    float &_att_coeff;
    float &_dec_coeff;
    float &_sus_level;
    float &_rel_coeff;
    float _output = 0;
    float _output_pregain = 0;
    float &_gain;
    float _output_gain = 0;
    float &_velocity;
    float &_velocity_sense;
    bool &_reset;
    bool &_drone;
    bool &_voice_gate_on;
    bool &_voice_trigger;
};

} // namespace Monique
