/**
 * @file NinaLfo.h
 * @brief Declaration of a simple LFO for nina
 * @date 2022-07-17
 *
 * Copyright (c) 2023 Melbourne Instruments
 *
 */
#pragma once

#include "MoniqueCommon.h"
#include "NoiseOscillator.h"
#include "SynthMath.h"
#include "common.h"

namespace Monique {

class MoniqueLfo {
    int counter = 0;

  public:
    MoniqueLfo(LfoOscShape &shape_a, LfoOscShape &shape_b, float &morph, float &slew, bool &lfo_retrigger, float &tempo, bool &unipolar, bool &sync, bool &global, float &global_out, bool set_global_out, bool &en_trigger) :
        _shape_a(shape_a), _shape_b(shape_b), _morph(morph), _slew_setting(slew), _lfo_retrigger(lfo_retrigger), _tempo(tempo), _lfo_unipolar_mode(unipolar), _tempo_sync(sync), _global(global), _global_out(global_out), _set_global_out(set_global_out), _enable_retrigger(en_trigger) {
        _noise.setVolume(1.0);
    };

    ~MoniqueLfo(){};

    void run() {
        if (_freq > 1.2) {
            //_freq = 1.2;
        }
        float phase_inc;
        float local_phase;

        if (_tempo_sync) {

            int time_multiplier_set = std::round(_freq * (float)TempoSyncMultipliers::NUM_SYNC_TEMPOS);
            float time_multiplier = _tempo_calc.getTempoSyncMultiplier(time_multiplier_set);
            float time = (_tempo * time_multiplier);
            float lfo_rate = time;
            phase_inc = 2.0f * M_PI * lfo_rate / (float)CV_SAMPLE_RATE;
            // printf("\n tmpo %f %f", _tempo, time);
        } else {
            if (_freq < 0.0001) {
                phase_inc = 0;
            } else {
                phase_inc = (2 * M_PI * (0.03 + 30. * _freq * _freq)) / (float)CV_SAMPLE_RATE;
            }
        }
        if (dump) {
            // printf("\n %d  %f", _tempo_sync, _freq);
        }

        _phase += phase_inc;
        _phase = _phase > 2 * M_PIf32 ? _phase - 2.f * M_PIf32 : _phase;
        local_phase = _phase;

        _sine = sinf32(local_phase);
        float square = _sine > 0 ? 1.f : -1.f;
        bool rand_step = (square > 0.f) != (_square > 0.f);
        _square = square;
        _saw_u = (local_phase / M_PIf32) - 1.f;
        _saw_d = 1.0 - (local_phase / M_PIf32);
        _tri = _saw_u > 0 ? 2.f * (.5f - _saw_u) : 2.f * (_saw_u + .5f);
        _rand = (rand_step || _retrigger_lfo) ? (_noise.getSample() - .5) * 2.0f : _rand;
        // Calculate the LFO value at the current phase
        _selectLFO(_shape_a, _lfo_a);
        _selectLFO(_shape_b, _lfo_b);

        float lfo;
        if (_lfo_unipolar_mode) {

            lfo = 0.5 * (_gain_a * _lfo_a + _gain_b * _lfo_b) + 0.5;
        } else {
            lfo = (_gain_a * _lfo_a + _gain_b * _lfo_b);
        }
        _lfo_mix += std::clamp(lfo - _lfo_mix, -_slew, _slew);
        _lfo_out = cv_clip(_lfo_mix * _gain / 2.f);

        if (_set_global_out) {
            _global_out = _lfo_out;
        }
        if (_global) {
            _lfo_out = _global_out;
        }
        _phase_inc = phase_inc;
        _retrigger_lfo = false;
    }

    void reCalculate() {
        // we have to stop the phase value overflowing. but it can go over 1.0. so we only do this here in recalculate, rather than in the run func
        _gain_a = 1 - _morph;
        _gain_b = _morph;
        if (_lfo_retrigger && _enable_retrigger) {
            _retrigger_lfo = true;
            _phase = 0;
        }

        // calculate slew value, this is scaled to the current LFO rate, and the current slew setting
        float phase_inc = (2 * M_PI * (0.03 + 30. * _freq * _freq)) / (float)CV_SAMPLE_RATE;
        float slew_linear_factor = 0;
        if (_slew_setting > 0.75) {
            slew_linear_factor = (_slew_setting - 0.75) * 4;
            slew_linear_factor = slew_linear_factor * slew_linear_factor;
        }
        _slew = _slew_setting * _slew_setting * phase_inc * 25 + slew_linear_factor;
    }

    float *getOutput() {
        return &_lfo_out;
    }

    float *getPitchIn() {
        return &_freq;
    }

    float *getGainIn() { return &_gain; }

    bool dump = false;

  private:
    /**
     * @brief lfo slew rate control is defined by Fc = 2^(15x - 2) where x = (0,1)
     *
     */
    static constexpr float lfo_shape_calc = (float)LfoOscShape::NUM_SHAPES;

    LfoOscShape &_shape_a, &_shape_b;
    const float &_morph;
    float &_slew_setting;
    float _slew = 0;
    float _phase_inc = 0;
    const bool &_lfo_retrigger;
    const bool &_enable_retrigger;
    const float &_tempo;
    const bool &_tempo_sync;
    const bool &_global;
    float &_global_out;
    const bool _set_global_out;
    float _gain = 1.0;
    float _phase = 0.0f;
    float _freq = 0.0;
    NoiseOscillator _noise;
    static TempoSyncMultipliers _tempo_calc;
    bool _phase_wrapped = false;
    bool _retrigger_lfo = false;
    float _sine;
    float _square;
    float _tri;
    float _saw_d;
    float _saw_u;
    float _rand;
    float _lfo_out = 0.0f;
    float _gain_a = 1.0f;
    float _gain_b = 0.0f;
    float _lfo_a = 0.0f;
    float _lfo_b = 0.0f;
    float _lfo_slew = 1.0f;
    float _lfo_mix = 0;
    bool &_lfo_unipolar_mode;

  public:
    const float phase_factor = exp2f(2.0f * M_PIf32 / CV_SAMPLE_RATE);

    /**
     * @brief calculate the output signal for the selected LFO shape and phase
     *
     * @param lfo_shape selected shape
     * @param phase current lfo phase
     * @param wrap if the phase has just wrapped
     * @param output current output. output signal is written to this variable
     * @return float
     */
    void _selectLFO(LfoOscShape lfo_shape, float &output) {
        switch (lfo_shape) {
        case LfoOscShape::SINE:
            output = _sine;
            break;

        case LfoOscShape::TRIANGLE:
            output = _tri;
            break;

        case LfoOscShape::SAWTOOTH_UP:
            output = _saw_u;
            break;
        case LfoOscShape::SAWTOOTH_DOWN:
            output = _saw_d;
            break;

        case LfoOscShape::SQUARE:
            output = _square;
            break;

        case LfoOscShape::RANDOM:
            output = _rand;
            break;

        default: {
        }
        }
    }
};

} // namespace Monique