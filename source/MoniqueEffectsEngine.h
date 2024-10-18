#pragma once
#include "ChorusEngine.h"
#include "MoniqueCommon.h"
#include "MoniqueDelay.h"
#include "MoniqueExpander.h"
#include "MoniqueReverb.h"
#include <array>

namespace Monique {
constexpr float metrenome_freq = 1500.f;
constexpr float metrenome_freq_2 = 2500.f;
constexpr float metrenome_pulse_cycles = 6.f;
constexpr int met_samples = (int)(((float)SAMPLE_RATE / metrenome_freq) * metrenome_pulse_cycles);
constexpr float met_volume = .095;
static const auto macro_params = fx_macro_params();

class FxProcessor {
  public:
    FxProcessor(ChorusEngine &ce, BBDDelay &delay_l, BBDDelay &delay_r, MoniqueReverb &reverb) :
        _chorus(ce),
        _delay_l(delay_l),
        _delay_r(delay_r),
        _reverb(reverb) {
        _reverb.setwetDry(1.0);
        _delay_l.setLFOGain(.2);
        _delay_r.setLFOGain(.2);
        _delay_l.setLFORate(.5);
        _delay_r.setLFORate(.6);
    }

    inline void run(audioBuffer &in_l, audioBuffer &in_r, audioBuffer &out_l, audioBuffer &out_r) {
        switch (fx_type) {
        case FxType::CHORUS: {
            for (int i = 0; i < BUFFER_SIZE; ++i) {
                _chorus.process(in_l.data() + i, in_r.data() + i, out_l.data() + i, out_r.data() + i);
            }
            break;
        }

        case FxType::DELAY: {
            _delay_l.run(in_l.data(), out_l.data());
            _delay_r.run(in_r.data(), out_r.data());
            break;
        }

        case FxType::REVERB: {
            //_limiter_2.run(fx_in_l, fx_in_r, fx_in_l, fx_in_r);
            float *in[2] = {in_l.data(), in_r.data()};
            float *out[2] = {out_l.data(), out_r.data()};
            _reverb.run(in, out, BUFFER_SIZE);
            break;
        }

        default:
            out_l.fill(0);
            out_r.fill(0);
            break;
        }
    }

    FxType fx_type = FxType::NONE;
    FxSlot fx_slot = FxSlot::ONE;

  private:
    ChorusEngine &_chorus;
    BBDDelay &_delay_l;
    BBDDelay &_delay_r;
    MoniqueReverb &_reverb;
};

class MoniqueFx {
  public:
    int counter = 0;

    inline void run(audioBuffer &in_l, audioBuffer &in_r, audioBuffer &out_l, audioBuffer &out_r) {

        // Apply the send level
        float fx_send = std::clamp(_fx_send + _fx_send_level_mod, -1.f, 1.f);
        fx_send = fx_send * fx_send;
        float inc = (fx_send - _fx_send_smooth) / BUFFER_SIZE;
        for (int i = 0; i < BUFFER_SIZE; ++i) {
            _input_l[i] = in_l[i];
            _input_r[i] = in_r[i];
        }

        _expander.process(_input_l, _input_r);
        // to avoid clicks and pops, we duck the output when configuring fx
        if (_output_duck) {
            _fx_send_smooth = 0.f;
            inc = 0;
            _output_duck_smooth = 1.f;
            _duck_counter++;

            if (_duck_counter > (.1f * BUFFER_RATE)) {
                _duck_counter = 0;
                _output_duck = false;
                _delay_l.time_skip();
                _delay_r.time_skip();
            }
            if (_duck_counter < 5) {

                _delay_l.reset();
                _delay_r.reset();
            }
        }
        for (int i = 0; i < BUFFER_SIZE; ++i) {
            _input_l[i] = _input_l[i] * _fx_send_smooth;
            _input_r[i] = _input_r[i] * _fx_send_smooth;
            _fx_send_smooth += inc;
            out_l[i] = _input_l[i];
            out_r[i] = _input_r[i];
        }
        counter++;

        // preset the fx mix parameters which can be overwridden in the macro processing
        _fx_1_mix = _fx_1_mix_param;
        _fx_2_mix = _fx_2_mix_param;

        // Apply the FX Macro Level mod
        applyFxMacroLevelMod();

        // Parallel mode?
        if (_fx_1.fx_slot == _fx_2.fx_slot) {
            if ((counter % 100) == 0) {
                //    printf("\n %f %f", g1, g2);
            }

            // Parallel - both use dry signal
            _fx_1.run(_input_l, _input_r, _left_buffer, _right_buffer);
            _fx_2.run(_input_l, _input_r, out_l, out_r);
            float mix_1_inc = (_fx_1_mix - _fx_1_mix_smooth) / BUFFER_SIZE;
            float mix_2_inc = (_fx_2_mix - _fx_2_mix_smooth) / BUFFER_SIZE;
            // mix output buffers
            for (int i = 0; i < BUFFER_SIZE; ++i) {
                float g1 = _fx_1_mix_smooth + i * mix_1_inc;
                float g2 = _fx_2_mix_smooth + i * mix_2_inc;
                out_l.at(i) = out_l.at(i) * g2 + _left_buffer.at(i) * g1;
                out_r.at(i) = out_r.at(i) * g2 + _right_buffer.at(i) * g1;
            }
            _fx_1_mix_smooth = _fx_1_mix;
            _fx_2_mix_smooth = _fx_2_mix;

        } else {
            // Series - run them in sequence
            if (_fx_1.fx_slot < _fx_2.fx_slot) {

                float mix_1_inc = (_fx_1_mix - _fx_1_mix_smooth) / BUFFER_SIZE;
                float mix_2_inc = (_fx_2_mix - _fx_2_mix_smooth) / BUFFER_SIZE;

                // Run FX 2 on FX 1 output
                _fx_1.run(_input_l, _input_r, out_l, out_r);

                for (int i = 0; i < BUFFER_SIZE; ++i) {

                    float g1 = _fx_1_mix_smooth + i * mix_1_inc;
                    float dry = 1.0 - g1;
                    out_l.at(i) = _input_l.at(i) * dry + out_l.at(i) * g1;
                    out_r.at(i) = _input_r.at(i) * dry + out_r.at(i) * g1;
                }

                _fx_2.run(out_l, out_r, _slot_l, _slot_r);

                for (int i = 0; i < BUFFER_SIZE; ++i) {

                    float g2 = _fx_2_mix_smooth + i * mix_2_inc;
                    float dry = 1.0 - g2;
                    out_l.at(i) = out_l.at(i) * dry + _slot_l.at(i) * g2;
                    out_r.at(i) = out_r.at(i) * dry + _slot_r.at(i) * g2;
                }

                _fx_1_mix_smooth = _fx_1_mix;
                _fx_2_mix_smooth = _fx_2_mix;
            } else {
                float mix_1_inc = (_fx_1_mix - _fx_1_mix_smooth) / BUFFER_SIZE;
                float mix_2_inc = (_fx_2_mix - _fx_2_mix_smooth) / BUFFER_SIZE;

                // Run FX 1 on FX 2 output
                _fx_2.run(_input_l, _input_r, out_l, out_r);

                for (int i = 0; i < BUFFER_SIZE; ++i) {

                    float g1 = _fx_2_mix_smooth + i * mix_2_inc;
                    float dry = 1.0 - g1;
                    out_l.at(i) = _input_l.at(i) * dry + out_l.at(i) * g1;
                    out_r.at(i) = _input_r.at(i) * dry + out_r.at(i) * g1;
                }

                _fx_1.run(out_l, out_r, _slot_l, _slot_r);

                for (int i = 0; i < BUFFER_SIZE; ++i) {

                    float g2 = _fx_1_mix_smooth + i * mix_1_inc;

                    float dry = 1.0 - g2;
                    out_l.at(i) = out_l.at(i) * dry + _slot_l.at(i) * g2;
                    out_r.at(i) = out_r.at(i) * dry + _slot_r.at(i) * g2;
                }
                _fx_1_mix_smooth = _fx_1_mix;
                _fx_2_mix_smooth = _fx_2_mix;
            }
        }
        _output_duck_smooth -= 0.005 * _output_duck_smooth;
        float output_level = 1.f - _output_duck_smooth;
        for (int i = 0; i < BUFFER_SIZE; i++) {
            out_l.at(i) *= output_level;
            out_r.at(i) *= output_level;
        }
        if (_fx_1.fx_type == _fx_2.fx_type) {
            out_l.fill(0.f);
            out_r.fill(0.f);
        }

        // if a metrenome pulse has been triggered, then output a sine pulse at a defined freq for a defined period (matching the reaper metrenome)
        if (_trig_met) {
            for (int i = 0; i < BUFFER_SIZE; ++i) {
                _met_counter++;
                if (_met_counter > met_samples) {
                    _met_counter = 0;
                    _trig_met = false;
                }
                float met_freq_sel;
                if (_met_freq == 1) {
                    met_freq_sel = metrenome_freq_2;
                } else {
                    met_freq_sel = metrenome_freq;
                }
                float phase = ((float)_met_counter * met_freq_sel) / ((float)SAMPLE_RATE);
                float sample = met_volume * std::sin(phase * 2.f * M_PI);
                out_l.at(i) += sample;
                out_r.at(i) += sample;
            }
        }
    }

    void setModDests(float fx_send_level_mod, float fx_macro_level_mod) {
        _fx_send_level_mod = fx_send_level_mod;
        _fx_macro_level_mod = fx_macro_level_mod;
    }

    inline void applyFxMacroLevelMod() {
        switch (_macro_select_param) {

        case PresetCommonParameters::DELAY_FB:
            _delay_l.adjFB(_fx_macro_level_mod);
            _delay_r.adjFB(_fx_macro_level_mod);
            break;

        case PresetCommonParameters::DELAY_TIME: {
            _delay_l.adjTime(_fx_macro_level_mod);
            _delay_r.adjTime(_fx_macro_level_mod);
            break;
        }

        case PresetCommonParameters::DELAY_TONE:
            _delay_l.adjFilter(_fx_macro_level_mod);
            _delay_r.adjFilter(_fx_macro_level_mod);
            break;

        case PresetCommonParameters::REVERB_TIME:
            _reverb.adjDecay(_fx_macro_level_mod);
            break;

        case PresetCommonParameters::REVERB_PREDELAY:
            // no action, macro removed
            break;

        case PresetCommonParameters::REVERB_TONE:
            _reverb.adjTone(_fx_macro_level_mod);
            break;

        case PresetCommonParameters::REVERB_SHIMMER:
            _reverb.adjShimmer(_fx_macro_level_mod);
            break;

        case PresetCommonParameters::FX_1_MIX:
            _fx_1_mix = unipolarClip(_fx_1_mix_param + _fx_macro_level_mod);
            break;
        case PresetCommonParameters::FX_2_MIX:
            _fx_2_mix = unipolarClip(_fx_2_mix_param + _fx_macro_level_mod);
            break;
        default:
            break;
        }
    }

    int c2 = 0;

    inline void setTempo(float tempo) {
        _delay_l.setTempo(tempo);
        _delay_r.setTempo(tempo);
    }

    inline void processParamChange(const ParamDecoded &param, const uint param_id, const float value) {
        switch (param.p_type) {
        case ParamType::GLOBAL:
            switch (param.global_param) {
            case (int)GlobalParams::METRONOME_TRIGGER:
                if (value > 0.25) {
                    _met_freq = 1;
                } else {
                    _met_freq = 0;
                }
                _trig_met = true;
                _met_counter = 0;
                break;

            default:
                break;
            }
            break;

        case ParamType::PRESET_COMMON:
            switch (param.preset_common_param) {
            case PresetCommonParameters::FX_1_TYPE: {

                auto type = (FxType)(int)preset_common_from_normalised_float(PresetCommonParameters::FX_1_TYPE, value);
                _fx_1.fx_type = type;
                _output_duck = true;
                break;
            }

            case PresetCommonParameters::FX_2_TYPE: {
                c2++;

                auto type = (FxType)(int)preset_common_from_normalised_float(PresetCommonParameters::FX_2_TYPE, value);
                _fx_2.fx_type = type;
                _output_duck = true;
                break;
            }

            case PresetCommonParameters::FX_1_SLOT: {
                FxSlot slot = (FxSlot)(int)preset_common_from_normalised_float(PresetCommonParameters::FX_1_SLOT, value);
                _fx_1.fx_slot = slot;
                _left_buffer.fill(0.f);
                _right_buffer.fill(0.f);
                _output_duck = true;
                break;
            }

            case PresetCommonParameters::FX_2_SLOT: {
                FxSlot slot = (FxSlot)(int)preset_common_from_normalised_float(PresetCommonParameters::FX_2_SLOT, value);
                _fx_2.fx_slot = slot;
                _left_buffer.fill(0.f);
                _right_buffer.fill(0.f);
                _output_duck = true;
                break;
            }
            case PresetCommonParameters::FX_1_MIX: {
                _fx_1_mix_param = value;
                break;
            }
            case PresetCommonParameters::FX_2_MIX: {
                _fx_2_mix_param = value;
                break;
            }

            case PresetCommonParameters::FX_MACRO_SELECT: {
                const int macro_idx = (int)preset_common_from_normalised_float(PresetCommonParameters::FX_MACRO_SELECT, value);
                const auto fx_param_decoded = ParamDecoded(macro_params.at(macro_idx));
                if (fx_param_decoded.p_type == ParamType::PRESET_COMMON) {
                    _macro_select_param = fx_param_decoded.preset_common_param;
                }
                break;
            }

            case get_preset_common_id_from_mod_matrix(ModMatrixSrc::CONSTANT, PresetCommonModMatrixDst::FX_SEND_LEVEL): {
                _fx_send = value;
                break;
            }

            case PresetCommonParameters::CHORUS_MODE: {
                const ChorusModes mode = (ChorusModes)(int)preset_common_from_normalised_float(PresetCommonParameters::CHORUS_MODE, value);
                switch (mode) {
                case ChorusModes::I:
                    _chorus.setEnablesChorus(true, false);
                    break;

                case ChorusModes::II:
                    _chorus.setEnablesChorus(false, true);
                    break;

                case ChorusModes::I_II:
                    _chorus.setEnablesChorus(true, true);
                    break;

                default:
                    break;
                }
                break;
            }

            case PresetCommonParameters::DELAY_TIME: {
                _delay_l.setTime(value);
                _delay_r.setTime(value);
                break;
            }

            case PresetCommonParameters::DELAY_FB:
                _delay_l.setFB(value);
                _delay_r.setFB(value);
                break;

            case PresetCommonParameters::DELAY_MODE:
                // TODO
                break;

            case PresetCommonParameters::DELAY_TIME_SYNC: {
                const int val = preset_common_from_normalised_float(PresetCommonParameters::DELAY_TIME_SYNC, value);
                _delay_l.setTimeSync(val);
                _delay_r.setTimeSync(val);
            } break;

            case PresetCommonParameters::DELAY_SYNC:

                _delay_l.setTempoSync(value > PARAM_BOOL_COMP);
                _delay_r.setTempoSync(value > PARAM_BOOL_COMP);
                break;

            case PresetCommonParameters::DELAY_TONE:
                _delay_l.setFilter(value);
                _delay_r.setFilter(value);
                break;

            case PresetCommonParameters::REVERB_PRESET:
                _reverb.setPreset(preset_common_from_normalised_float(PresetCommonParameters::REVERB_PRESET, value));
                _output_duck = true;
                break;

            case PresetCommonParameters::REVERB_TIME:
                _reverb.setDecay(value);
                break;

            case PresetCommonParameters::REVERB_PREDELAY:
                _reverb.setPreDelay(value);
                _output_duck = true;
                break;

            case PresetCommonParameters::REVERB_TONE:
                _reverb.setTone(value);
                break;

            case PresetCommonParameters::REVERB_SHIMMER:
                _reverb.setShimmer(value);
                break;

            default:
                break;
            }
            break;

        default:
            break;
        }
    }

  private:
    bool _output_duck = false;
    int _duck_counter = 0;
    int _met_freq = 0;
    float _fx_1_mix = 1;
    float _fx_1_mix_param = 1;
    float _fx_2_mix_param = 1;
    float _output_duck_smooth = 1;
    float _fx_2_mix = 1;
    float _fx_1_mix_smooth = 1;
    float _fx_2_mix_smooth = 1;
    float _fx_send = 1.f;
    float _fx_send_smooth = 1.f;
    int _met_counter = 0;
    bool _trig_met = false;
    PresetCommonParameters _macro_select_param = PresetCommonParameters::FX_1_MIX;
    float _fx_send_level_mod = 0.f;
    float _fx_macro_level_mod = 0.f;
    ChorusEngine _chorus = ChorusEngine(SAMPLE_RATE);
    audioBuffer _left_buffer, _right_buffer;
    audioBuffer _slot_l, _slot_r;
    audioBuffer _input_l, _input_r;
    BBDDelay _delay_l;
    BBDDelay _delay_r;
    MoniqueReverb _reverb;
    FxType _fx_type_1, _fx_type_2;
    FxProcessor _fx_1 = FxProcessor(_chorus, _delay_l, _delay_r, _reverb);
    FxProcessor _fx_2 = FxProcessor(_chorus, _delay_l, _delay_r, _reverb);
    MoniqueExpander _expander;
};

} // namespace Monique
