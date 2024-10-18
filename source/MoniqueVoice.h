
#pragma once
#include "MoniqueAnalogModel.h"
#include "MoniqueCommon.h"
#include "MoniqueDriveCompensator.h"
#include "MoniqueEnvelope.h"
#include "MoniqueLfo.h"
#include "MoniqueMatrix.h"
#include "MoniqueOscMix.h"
#include "MoniqueOutputPanner.h"
#include "MoniqueSlaveEnvelope.h"
#include "MoniqueWavetable.h"
#include "common.h"

constexpr float VOICE_LOW_OUTPUT_THRESH = 0.1;

namespace Monique {

class MoniqueVoice {
  public:
    MoniqueVoice(const uint voice_num, VoiceData &voice_data, MoniqueSynthParameters &parameters, WavetableLoader &loader_a, WavetableLoader &loader_b, GlobalSynthParameters &gp, DriftSource &ds) :
        _voice_num(voice_num),
        _voice_data(voice_data),
        _parameters(parameters),

        _amp_vel_sense(_morphed_parameters.at((int)StateParameters::AMP_VEL)),
        _filt_vel_sense(_morphed_parameters.at((int)StateParameters::FILT_VEL)),
        _aux_vel_sense(_morphed_parameters.at((int)StateParameters::AUX_VEL)),
        _unison_detune_amt(_morphed_parameters.at((int)StateParameters::UNISON_DETUNE)),
        _glide_rate(_morphed_parameters.at((int)StateParameters::GLIDE_AMT)),
        _wavetable(WavetableOsc(loader_a, loader_b, _voice_morph, parameters.wt_interpolate, parameters.wt_slow_mode, _voice_data)),
        _midi_cc_1_src(_parameters.common_params.at((int)CommonParameters::MIDI_CC_1_MOD_SOURCE)),
        _midi_cc_2_src(_parameters.common_params.at((int)CommonParameters::MIDI_CC_2_MOD_SOURCE)),
        _gp(gp),
        _drift(ds),
        _osc_1_drift(ds.getOscSource()),
        _osc_2_drift(ds.getOscSource()),
        _osc_3_drift(ds.getOscSource()),
        _osc_unison_drift(ds.getOscSource())

    {

        // printf("\ninit voice %d", _voice_num);
        _matrixSetup();
        if (voice_num == 0) {
            _lfo_1.dump = true;
            _amp_env.print = true;
            _panner.print = true;
            _drive_compensator.print();
            _matrix.dump = true;
            _wavetable.dump = true;
        }
        _wavetable.setOutput(&_voice_data.wt_buffer);
        _voice_data.fm_1_1 = &fm_1_1;
        _voice_data.fm_1_2 = &fm_1_2;
        _voice_data.fm_2_1 = &fm_2_1;
        _voice_data.fm_2_2 = &fm_2_2;
        _drift_data.fill(0.f);
        if (voice_num == 0) {
            _lfo_1.dump = true;
        }
    }

    bool blacklisted() { return false; }

    void reloadCal() {
        _voice_model.reloadCal();
    }

    int setNoteOn(MidiNote note, float pan, float detune) {
        _pan_src = pan;
        setNoteOn(note);
        _unison_detune = detune;

        return _voice_num;
    }

    void paraModeSetup() {

        if (_parameters.para_mode == ParaphonicModes::DISABLE) {

            // if in disable mode, clear all the unison data
            _voice_data.unison_1_tmp = 1.f;
            _voice_data.unison_2_tmp = 0.f;
            _voice_data.unison_osc_enable = false;
            _unison_gate = false;
        } else {

            _voice_data.unison_osc_enable = true;
        }
    }

    bool realVoiceLastAlloc() {
        return _real_voice_is_last_alloc;
    }

    void setUnisonNote(MidiNote note) {
        _real_voice_is_last_alloc = false;
        float pitch = (midiNoteToCv(note.pitch + _parameters.octave_offset) / NOTE_GAIN);
        _unison_offset = pitch - _keyboard_pitch;
        _key_velocity_2 = note.velocity;
        _unison_gate = true;
        _unison_trigger = true;
        _current_unison = note;

        if (_parameters.para_mode == ParaphonicModes::DIST_RETRIG || _parameters.para_mode == ParaphonicModes::TIME_RETRIG) {
            if (_filt_env._env_state == AdsrState::REL || _filt_env._env_state == AdsrState::DEC) {
                _voice_trigger = true;
            }
        }
        // printf("unison note on %d \n", note.pitch);
    }

    void setUnisonNote(MidiNote note, const float detune) {
        _unison_voice_detune = detune;
        setUnisonNote(note);
    }

    int setNoteOn(MidiNote note) {
        if (!(_parameters._voice_mode == VoiceMode::LEGATO && _note_gate)) {
            _voice_trigger = true;
            _main_note_trigger = true;
            _buffer_counter = 0;
        }
        _note_gate = true;
        _real_voice_is_last_alloc = true;
        _pending_release = false;
        _key_velocity = note.velocity;
        _current_note = note;
        _keyboard_pitch = (midiNoteToCv(note.pitch + _parameters.octave_offset) / NOTE_GAIN);

        float pitch = (midiNoteToCv(_current_unison.pitch + _parameters.octave_offset) / NOTE_GAIN);
        _unison_offset = pitch - _keyboard_pitch;

        // if portamento glide is enabled, and we are currently in the release state, then we reset the glide filter to stop the voice from gliding
        if (_amp_env.inReleaseState() && _parameters.glide_mode == GlideModes::PORTAMENTO_LINEAR) {
            _keyboard_pitch_glide = (midiNoteToCv(note.pitch + _parameters.octave_offset) / NOTE_GAIN);
        }
        // printf("real note on %d \n", note.pitch);
        return _voice_num;
    }

    void getFxMod(float &send_mod, float &macro_mod) {
        send_mod = _fx_send_mod;
        macro_mod = _fx_macro_mod;
    }

    MidiNote &currentNote() {

        return _current_note;
    }

    MidiNote &currentUnison() {

        return _current_unison;
    }

    bool inReleased() {
        return _amp_env.inReleaseState();
    }

    bool pendingRelease() {
        return _pending_release;
    }

    bool pendingReleaseUni() {
        return _pending_release_uni;
    }

    void setUnalloc(bool set) {
    }

    void setPolyAT(const PolyPressure &pp) {
        // printf("\n %d %f", _voice_num, pp.pressure);
        _poly_at_level = pp.pressure;
    }

    bool probablyOff() {
        return _probably_off;
    }

    int setNoteOff(MidiNote note) {
        _release_velocity = note.note_off_vel;
        _pending_release = true;
        // printf(" free %d %d \n ", _voice_num, note.pitch);
        return _voice_num;
    }

    int setUnisonOff(MidiNote note) {
        _release_velocity_uni = note.note_off_vel;
        _pending_release_uni = true;
        //_release_velocity = note.note_off_vel;
        // printf(" free unison%d %d \n ", _voice_num, note.pitch);
        return _voice_num;
    }

    void setUnisonOff() {
        _release_velocity_uni = 0;
        _pending_release_uni = true;
        //_release_velocity = note.note_off_vel;
        // printf(" free unison%d %d \n ", _voice_num, note.pitch);
    }

    void setMpeData(const MpeChannelData *mpe_data) {
        _mpe_data = mpe_data;
    }

    bool isGateOn() { return _voice_gate_on; };

    bool mainGateOn() { return _note_gate; };

    bool unisonGateOn() { return _unison_gate; };

    bool isAtt() {
        return _amp_env._env_state == AdsrState::ATT;
    }

    int buffersSinceTrigger() {
        return _buffer_counter;
    }

    void forceNoteOff() {
        _pending_release = true;
        _pending_release_uni = true;
        _parameters.sustain_ped = false;
        _voice_gate_on = false;
        _note_gate = false;
        _unison_gate = false;
        // printf("\n force note off\n");
    }

    int counter = 0;

    void _parameter_morphing() {
        float morph = (_parameters.morph_value_smooth + _morph_dest);
        morph = unipolarClip(morph);
        _voice_morph = morph;
        const float inv_morph = 1 - morph;
        const float modwheel = _parameters.common_params.at(CommonParameters::MIDI_MODWHEEL);
        const auto &state = _parameters.layer_state;
        _modwheel = mpeParamSmooth(modwheel, _modwheel);
        for (int i = 0; i < TOTAL_STATE_PARAMS; ++i) {
            float a;
            float b;
            if (_parameters.state_morphed_blocked.at(i)) {
                if (state == State::A) {
                    a = _parameters.state_a_smooth.at(i);
                    b = a;
                } else {
                    b = _parameters.state_b_smooth.at(i);
                    a = b;
                }
            } else {
                a = _parameters.state_a_smooth.at(i);
                b = _parameters.state_b_smooth.at(i);
            }

            // perform weighted sum morph
            const float tmp = a * inv_morph + b * morph;

            if (i == get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::FX_SEND)) {
                // printf("\n %f", tmp);
            }
            switch (i) {
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_1, ModMatrixDst::OSC_2_PITCH):
                fm_1_2 = from_normalised_float(get_state_id_from_mod_matrix(ModMatrixSrc::OSC_1, ModMatrixDst::OSC_2_PITCH), tmp);
                _morphed_parameters.at(i) = 0;
                break;
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_2, ModMatrixDst::OSC_1_PITCH):
                fm_2_1 = from_normalised_float(get_state_id_from_mod_matrix(ModMatrixSrc::OSC_2, ModMatrixDst::OSC_1_PITCH), tmp);
                _morphed_parameters.at(i) = 0;
                break;
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_1, ModMatrixDst::OSC_1_PITCH):
                fm_1_1 = from_normalised_float(get_state_id_from_mod_matrix(ModMatrixSrc::OSC_1, ModMatrixDst::OSC_1_PITCH), tmp);
                _morphed_parameters.at(i) = 0;
                break;
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_2, ModMatrixDst::OSC_2_PITCH):
                fm_2_2 = from_normalised_float(get_state_id_from_mod_matrix(ModMatrixSrc::OSC_2, ModMatrixDst::OSC_2_PITCH), tmp);
                _morphed_parameters.at(i) = 0;
                break;
                // scaling for envelops  & LFO's must be performed at the actual envelop gen
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::DRIVE):
                _morphed_parameters.at(i) = driveAdjust(tmp);
                break;
            case (int)StateParameters::LFO_1_SYNC_RATE:
            case (int)StateParameters::LFO_2_SYNC_RATE:
            case (int)StateParameters::LFO_3_SYNC_RATE:
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_1_RATE):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_2_RATE):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_3_RATE):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::FILTER_EG_ATT)... get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::FILTER_EG_REL):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::AMP_EG_ATT)... get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::AMP_EG_REL):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::AUX_EG_ATT)... get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::AUX_EG_REL):
                _morphed_parameters.at(i) = (tmp);
                /* code */
                break;
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_1, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_2, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_3, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AMP_EG, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_3, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::FILTER_EG, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AFTERTOUCH, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::KEY_VELOCITY, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::PANPOSITION, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AUX_EG, ModMatrixDst::OSC_1_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::OFFSET, ModMatrixDst::OSC_1_PITCH):

            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_1, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_2, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_3, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AMP_EG, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_3, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::FILTER_EG, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AFTERTOUCH, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::KEY_VELOCITY, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::PANPOSITION, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AUX_EG, ModMatrixDst::OSC_2_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::OFFSET, ModMatrixDst::OSC_2_PITCH):

            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_1, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_2, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::LFO_3, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AMP_EG, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::OSC_3, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::FILTER_EG, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AFTERTOUCH, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::KEY_VELOCITY, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::PANPOSITION, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::AUX_EG, ModMatrixDst::OSC_3_PITCH):
            case get_state_id_from_mod_matrix(ModMatrixSrc::OFFSET, ModMatrixDst::OSC_3_PITCH):
                _morphed_parameters.at(i) = xAbsParameterCurve(tmp);
                break;

            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_1_LEVEL):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_2_LEVEL):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_3_LEVEL):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_4_LEVEL):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_1_GAIN):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_2_GAIN):
            case get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_3_GAIN):
                _morphed_parameters.at(i) = xSquaredParamCurve(tmp);
                break;

            default:
                _morphed_parameters.at(i) = state_param_from_normalised_float(i, tmp);
                break;
            }
        }
        // update osc tune vars
        constexpr float log2_440_over_32 = std::log2f(440.f / 32.f);
        constexpr float FUDGE = 0.00318;
        constexpr float interval2 = (log2_440_over_32 + ((float)14 - 9.f) / 12.f);
        constexpr float interval1 = (log2_440_over_32 + ((float)15 - 9.f) / 12.f);
        constexpr float interval = (6 + 1) * (interval1 - interval2);
        constexpr float gain = (1 * interval) / (NOTE_GAIN * FINE_CENTS_RANGE);
        constexpr float offset = -gain / 2 - 2.f / NOTE_GAIN + FUDGE;
        constexpr float lin = 1;
        constexpr float scale = 1. / (1 + lin);
        float fine;
        float coarse;
        fine = _morphed_parameters.at((int)StateParameters::OSC_1_FINE);
        coarse = _morphed_parameters.at((int)StateParameters::OSC_1_COARSE);
        fine = (fine * gain + offset);
        coarse = coarse / (NOTE_GAIN);
        if (_voice_num == 0) {
            // printf("\n %f %f ", _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::AMP_EG_REL)), 0);
        }
        _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_1_PITCH)) = fine + coarse;
        fine = _morphed_parameters.at((int)StateParameters::OSC_2_FINE);
        coarse = _morphed_parameters.at((int)StateParameters::OSC_2_COARSE);
        fine = (fine * gain + offset);
        coarse = coarse / NOTE_GAIN;
        _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_2_PITCH)) = fine + coarse;
        fine = _morphed_parameters.at((int)StateParameters::OSC_3_FINE);
        coarse = _morphed_parameters.at((int)StateParameters::OSC_3_COARSE);
        fine = (fine * gain + offset);
        coarse = coarse / NOTE_GAIN;
        _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::OSC_3_PITCH)) = fine + coarse;
        _voice_data.sub_osc = (int)state_param_from_normalised_float((int)StateParameters::SUB_OSC, _morphed_parameters.at((int)StateParameters::SUB_OSC));
        _noise_mode = (NoiseModes)(int)_morphed_parameters.at((int)StateParameters::OSC_4_MODE);
        _voice_data.aux_level = _parameters.ext_in_gain;
        const bool lp_2_pole = _morphed_parameters.at((int)StateParameters::FILTER_LP_MODE) < PARAM_BOOL_COMP;
        _voice_data.filter_2_slope = (FilterSlope)(!lp_2_pole);

        // update hard sync
        const bool hard_sync = _morphed_parameters.at((int)StateParameters::HARD_SYNC) > PARAM_BOOL_COMP;
        _voice_data.hard_sync = hard_sync;

        // update env slow mode
        _aux_slow_mode = _morphed_parameters.at((int)StateParameters::AUX_ENV_SLOW) > PARAM_BOOL_COMP;
        _aux_env.setSlowMode(_aux_slow_mode);

        // update envelope parameters
        _filt_env_reset = _morphed_parameters.at((int)StateParameters::FILT_ENV_RESET) > PARAM_BOOL_COMP;
        _amp_env_reset = _morphed_parameters.at((int)StateParameters::AMP_ENV_RESET) > PARAM_BOOL_COMP;
        // TODO: _vca_drone = _morphed_parameters.at((int)StateParameters::AMP_ENV_RESET) > PARAM_BOOL_COMP;
        _aux_env_reset = _morphed_parameters.at((int)StateParameters::AUX_ENV_RESET) > PARAM_BOOL_COMP;

        // update lfo unipolar state
        _lfo_1_unipolar = _morphed_parameters.at((int)StateParameters::LFO_1_UNIPOLAR) > PARAM_BOOL_COMP;
        _lfo_2_unipolar = _morphed_parameters.at((int)StateParameters::LFO_2_UNIPOLAR) > PARAM_BOOL_COMP;
        _lfo_2_unipolar = _morphed_parameters.at((int)StateParameters::LFO_3_UNIPOLAR) > PARAM_BOOL_COMP;
        _lfo_1_trigger = _morphed_parameters.at((int)StateParameters::LFO_1_RETRIGGER) > PARAM_BOOL_COMP;
        _lfo_2_trigger = _morphed_parameters.at((int)StateParameters::LFO_2_RETRIGGER) > PARAM_BOOL_COMP;
        _lfo_3_trigger = _morphed_parameters.at((int)StateParameters::LFO_3_RETRIGGER) > PARAM_BOOL_COMP;

        if (_parameters.lfo_1_sync) {
            _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_1_RATE)) = _morphed_parameters.at((int)StateParameters::LFO_1_SYNC_RATE);
        }
        if (_parameters.lfo_2_sync) {
            _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_2_RATE)) = _morphed_parameters.at((int)StateParameters::LFO_2_SYNC_RATE);
        }
        if (_parameters.lfo_3_sync) {
            _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::LFO_3_RATE)) = _morphed_parameters.at((int)StateParameters::LFO_3_SYNC_RATE);
        }
    }

    int c_v = 0;

    void runDriftVintage() {
        const float &dg = _parameters.layer_drift;
        const float &vg = _parameters.layer_vintage;
        _drift_data.fill(0.f);

        _drift_data.at((int)ModMatrixDst::OSC_1_PITCH) = _osc_1_drift * dg;
        _drift_data.at((int)ModMatrixDst::OSC_2_PITCH) = _osc_2_drift * dg;
        _drift_data.at((int)ModMatrixDst::OSC_3_PITCH) = _osc_3_drift * dg;
        _unison_drift_vintage = _osc_unison_drift * dg;
        //    if (_voice_num == 0 && (c_v++ % 1000 == 0))
        // printf("\n %f %f %f ", _osc_1_drift, dg, _drift_data.at(0));
    }

    void setOutMode() {
        switch (_parameters.output_mode) {
        case OutputRouting::ONE:
            _voice_data.mute_1 = false;
            _voice_data.mute_2 = true;
            break;
        case OutputRouting::TWO:
            _voice_data.mute_1 = true;
            _voice_data.mute_2 = false;
            break;
        case OutputRouting::ONE_TWO:
            _voice_data.mute_1 = false;
            _voice_data.mute_2 = false;
            break;
        default:
            break;
        }
    }

    void run() {
        _buffer_counter++;
        _voice_gate_on = _note_gate || _unison_gate;
        _parameter_morphing();
        runDriftVintage();
        auto &_glide_mode = _parameters.glide_mode;
        float g_tmp = _glide_rate * _glide_rate * _glide_rate;
        g_tmp = (1 - fastpow2(-(1.4 / ((float)BUFFER_RATE * 5. * g_tmp))));
        if (_glide_mode == GlideModes::LOG) {
            _keyboard_pitch_glide += (_keyboard_pitch - _keyboard_pitch_glide) * g_tmp;
        } else {
            const float glide_rate = g_tmp / 10;
            // linear glide mode
            float glide_tmp = _keyboard_pitch - _keyboard_pitch_glide;
            glide_tmp = glide_tmp > glide_rate ? glide_rate : glide_tmp;
            glide_tmp = glide_tmp < -glide_rate ? -glide_rate : glide_tmp;
            _keyboard_pitch_glide = glide_tmp + _keyboard_pitch_glide;
        }

        float at_sum = unipolarClip(_poly_at_level + _parameters.aftertouch);
        _aftertouch_src = mpeParamSmooth(at_sum, _aftertouch_src);
        float val = _morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc::CONSTANT, ModMatrixDst::FILTER_CUTOFF_LP));
        // printf("\n cut %f ", val);
        const float detune_scale = _unison_detune_amt * _unison_detune_amt;
        _keyboard_src = _keyboard_pitch_glide + _parameters.pitch_bend + _unison_detune * detune_scale / NOTE_GAIN + _parameters.master_detune;
        _voice_data.noise_mode = _noise_mode;

        // unison voice offset
        _voice_data.unison_offset = _unison_offset + _unison_voice_detune * detune_scale / NOTE_GAIN + _unison_drift_vintage - (_unison_detune * detune_scale / NOTE_GAIN);
        // TODO: MPE
        // TODO: Rel Vel
        // TODO: cv input

        counter++;
        if (counter > 300) {
            counter = 0;
            _osc_mixer_1.print();
            // printf("\n %f %f %f", _parameters.state_a.at((uint)StateParameters::Lfo1Slew),_parameters.state_b.at((uint)StateParameters::Lfo1Slew),_parameters.state_params.at((uint)StateParameters::Lfo1Slew));
            if (_voice_num == 0) {
                float val = _morphed_parameters.at(get_state_id_from_mod_matrix(Monique::ModMatrixSrc::CONSTANT, Monique::ModMatrixDst::OSC_1_PITCH));
                // printf("\n %f %f %f", _key_velocity, _release_velocity, 0);
            }
        }
        setOutMode();
        _matrix.run_slow_slots();
        _lfo_1.reCalculate();
        _lfo_2.reCalculate();
        _lfo_3.reCalculate();
        _amp_env.reCalculate();
        _filt_env.reCalculate();
        _aux_env.reCalculate();
        _unison_env.reCalculate();
        _drive_compensator.reCalculate();
        _panner.reCalculate();
        _osc_mixer_1.reCalculate();
        _osc_mixer_2.reCalculate();
        _voice_data.post_hp_gain = _drive_comp_level;
        for (int i = 0; i < CV_BUFFER_SIZE; ++i) {
            _cv_a = _gp._cv_1.at(i);
            _cv_b = _gp._cv_2.at(i);
            _matrix.run_fast_slots();

            _lfo_1.run();
            _lfo_2.run();
            _lfo_3.run();
            _amp_env.run();
            _aux_env.run();
            _filt_env.run();
            _unison_env.run();
            _drive_compensator.run();
            *_wavetable.getWtPitch() = _voice_data.osc_3_pitch_tmp;
            if (_voice_num == 0) {
            }
            *_wavetable.getWtVol() = _voice_data.osc_3_lev_tmp;
            *_wavetable.getWtPosition() = _voice_data.osc_3_shape_tmp;
            _wavetable.run();

            *_pan_vca_in = *_amp_env_output;
            _panner.run();
            _osc_mixer_1.run();
            _osc_mixer_2.run();
            _voice_data.copy_to_buffer(i);
            _osc_1_val = _osc_1_cv.at(i);
            _osc_2_val = _osc_2_cv.at(i);
            _voice_trigger = false;
            _unison_trigger = false;
            _main_note_trigger = false;
            _probably_off = (_panner._max_level < VOICE_LOW_OUTPUT_THRESH);
            if (!_parameters.sustain_ped) {
                if (_pending_release) {
                    _pending_release = false;
                    _note_gate = false;
                }
                if (_pending_release_uni) {
                    _pending_release_uni = false;
                    _unison_gate = false;
                }
            }
        }
        _wavetable.reCalculate();
    }

    inline auto &getVoiceModel() {
        return _voice_model;
    }

    void setLastReleased(bool is_last) {
        _voice_data.last_allocated = is_last;
    }

  private:
    // class data
    DriftSource &_drift;
    GlobalSynthParameters &_gp;
    int _buffer_counter = 0;
    float &_glide_rate;
    float _unison_voice_detune = 0;
    float _unison_offset = 0;
    float _keyboard_pitch_glide = 0;
    float &_unison_detune_amt;
    bool _aux_slow_mode = false;
    bool _probably_off = false;
    const uint _voice_num;
    VoiceData &_voice_data;
    MoniqueSynthParameters &_parameters;
    std::array<float, TOTAL_STATE_PARAMS> _morphed_parameters = state_param_init();
    WavetableOsc _wavetable;
    NoiseModes _noise_mode = NoiseModes::XOR;
    float _voice_morph = 0;
    MidiNote _current_note;
    MidiNote _current_unison;
    float _keyboard_pitch = 0;
    float _modwheel = 0;
    const MpeChannelData *_mpe_data;
    const float _key_offset = 51.f / 127.f;
    float _max_mix_out_level = 0;
    float fm_1_1 = 0;
    float fm_1_2 = 0;
    float fm_2_2 = 0;
    float fm_2_1 = 0;
    bool _voice_gate_on = false;
    bool _pending_release = false;
    bool _pending_release_uni = false;
    bool _note_gate = false;
    bool _voice_trigger = false;
    bool _unison_gate = false;
    bool _unison_trigger = false;
    float _key_velocity = 0;
    float _key_velocity_2 = 0;
    float _release_velocity = 0;
    float _release_velocity_uni = 0;
    float &_amp_vel_sense;
    float &_filt_vel_sense;
    float &_aux_vel_sense;
    bool _filt_env_reset = false;
    bool _amp_env_reset = false;
    bool _aux_env_reset = false;
    bool _vca_drone = false;
    bool _vcf_drone = false;
    bool _aux_drone = false;
    bool _real_voice_is_last_alloc = false;
    bool _lfo_1_trigger = false;
    bool _lfo_2_trigger = false;
    bool _lfo_3_trigger = false;
    float _compression_signal = 0;
    float _constant = 1.0f;
    float _zero = 0.f;
    bool _lfo_1_unipolar = false;
    bool _lfo_2_unipolar = false;
    bool _lfo_3_unipolar = false;
    bool _main_note_trigger = false;
    MoniqueVoiceModel _voice_model = MoniqueVoiceModel(_voice_num, _voice_data, _wavetable);
    float &_drive_comp_level = _drive_compensator.getMixerVcaCompensation();
    DriveCompensator _drive_compensator = DriveCompensator(_parameters.overdrive, _compression_signal, _max_mix_out_level, _parameters.layer_volume, _voice_data.overdrive);

    OscMixer _osc_mixer_2 = OscMixer(_voice_data.tri_2_lev_tmp, _voice_data.sqr_2_lev_tmp, _voice_data.osc_2_shape_tmp, _voice_data.osc_2_blend_tri, _voice_data.osc_2_blend_sqr);
    OscMixer _osc_mixer_1 = OscMixer(_voice_data.tri_1_lev_tmp, _voice_data.sqr_1_lev_tmp, _voice_data.osc_1_shape_tmp, _voice_data.osc_1_blend_tri, _voice_data.osc_1_blend_sqr);
    GenAdsrEnvelope _amp_env = GenAdsrEnvelope(_key_velocity, _amp_vel_sense, _amp_env_reset, _vca_drone, _note_gate, _main_note_trigger);
    GenAdsrEnvelope _filt_env = GenAdsrEnvelope(_key_velocity, _filt_vel_sense, _filt_env_reset, _vcf_drone, _voice_gate_on, _voice_trigger);
    GenAdsrEnvelope _aux_env = GenAdsrEnvelope(_key_velocity, _aux_vel_sense, _aux_env_reset, _aux_drone, _voice_gate_on, _voice_trigger);
    SlaveEnvelope _unison_env = SlaveEnvelope(_amp_env, _key_velocity_2, _unison_gate, _unison_trigger);
    NinaOutPanner _panner = NinaOutPanner(&_voice_data.vca_l_tmp, &_voice_data.vca_r_tmp, _drive_compensator.getMixbusCompensation(), _max_mix_out_level, _voice_data, _key_velocity, _key_velocity_2, _amp_env, _unison_env, _parameters);
    float test1, test2;
    MoniqueLfo _lfo_1 = MoniqueLfo(_parameters.lfo_1_a_shape, _parameters.lfo_1_b_shape, _parameters.morph_value_smooth, _morphed_parameters.at((uint)StateParameters::LFO_1_SLEW), _voice_trigger, _parameters.tempo, _lfo_1_unipolar, _parameters.lfo_1_sync, _parameters.lfo_1_global, _parameters.lfo_1_global_phase, _voice_num == 0, _lfo_1_trigger);
    MoniqueLfo _lfo_2 = MoniqueLfo(_parameters.lfo_2_a_shape, _parameters.lfo_2_b_shape, _parameters.morph_value_smooth, _morphed_parameters.at((uint)StateParameters::LFO_2_SLEW), _voice_trigger, _parameters.tempo, _lfo_2_unipolar, _parameters.lfo_2_sync, _parameters.lfo_2_global, _parameters.lfo_2_global_phase, _voice_num == 0, _lfo_2_trigger);
    MoniqueLfo _lfo_3 = MoniqueLfo(_parameters.lfo_3_a_shape, _parameters.lfo_3_b_shape, _parameters.morph_value_smooth, _morphed_parameters.at((uint)StateParameters::LFO_3_SLEW), _voice_trigger, _parameters.tempo, _lfo_3_unipolar, _parameters.lfo_3_sync, _parameters.lfo_3_global, _parameters.lfo_3_global_phase, _voice_num == 0, _lfo_3_trigger);
    float _keyboard_src = 0;
    float _aftertouch_src = 0;
    float _poly_at_level = 0;
    float _channel_at_level = 0;
    float _pan_src = 0;
    float _tmp_xor_level = 0;
    float *_pan_vca_in = _panner.getVcaIn();
    float *_amp_env_output = _amp_env.getOutput();
    float _morph_dest = 0;
    float _osc_1_val, _osc_2_val;
    float _unison_detune = 0.f;
    float &_midi_cc_1_src;
    float &_midi_cc_2_src;
    float _fx_send_mod = 0;
    float _fx_macro_mod = 0;
    float _cv_a = 0;
    float _cv_b = 0;

    // osc drift sources
    float &_osc_1_drift;
    float &_osc_2_drift;
    float &_osc_3_drift;
    float &_osc_unison_drift;
    float _unison_drift_vintage = 0.f;

    std::array<float, CV_BUFFER_SIZE> &_osc_1_cv = _voice_model.getCvOsc1();
    std::array<float, CV_BUFFER_SIZE> &_osc_2_cv = _voice_model.getCvOsc2();

    std::array<float *, TOTAL_MOD_DSTS> _getDsts() { return _matrix_data._matrix_dests; };

    std::array<float *, TOTAL_MOD_SRC> _getSrcs() {
        return _matrix_data._matrix_srcs;
    };

    std::array<float *, FAST_SRCS> _getFastSrcs() {
        return _matrix_data._fast_src;
    }

    std::array<float *, FAST_SRCS * FAST_DESTS> _getFastGains() {
        return _matrix_data._fast_gains;
    }

    void _matrixSetup() {

        std::array<bool, TOTAL_MOD_SRC> src_is_fast = {false};
        std::array<bool, TOTAL_MOD_DSTS> dst_is_fast = {false};
        for (auto &item : src_is_fast) {
            item = false;
        }
        for (auto &item : dst_is_fast) {
            item = false;
        }
        auto fast_src_it = src_is_fast.begin();
        auto fast_dst_it = dst_is_fast.begin();

        // find the address of all the fast sources in order and add them to the fast gains array
        // add fast sources to matrix and to fast matrix src array
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::LFO_1) = _lfo_1.getOutput();
        src_is_fast.at((uint)ModMatrixSrc::LFO_1) = true;
        _lfo_1.getOutput();
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::LFO_2) = _lfo_2.getOutput();
        src_is_fast.at((uint)ModMatrixSrc::LFO_2) = true;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::LFO_3) = _lfo_3.getOutput();
        src_is_fast.at((uint)ModMatrixSrc::LFO_3) = true;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::AUX_EG) = _aux_env.getOutput();
        src_is_fast.at((uint)ModMatrixSrc::AUX_EG) = true;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::AMP_EG) = _amp_env.getOutput();
        src_is_fast.at((uint)ModMatrixSrc::AMP_EG) = true;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::FILTER_EG) = _filt_env.getOutput();
        src_is_fast.at((uint)ModMatrixSrc::FILTER_EG) = true;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::OSC_3) = _wavetable.getWtCvOut();
        src_is_fast.at((uint)ModMatrixSrc::OSC_3) = true;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::OSC_1) = &_osc_1_val;
        src_is_fast.at((uint)ModMatrixSrc::OSC_1) = true;

        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::OSC_2) = &_osc_2_val;
        src_is_fast.at((uint)ModMatrixSrc::OSC_2) = true;
        // add slow sources
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::KEY_PITCH) = &_keyboard_src;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::KEY_VELOCITY) = &_key_velocity;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::AFTERTOUCH) = &_aftertouch_src;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::PANPOSITION) = &_pan_src;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::CONSTANT) = &_constant;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::OFFSET) = &_constant;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::MODWHEEL) = &_modwheel;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::EXPRESSION) = &_parameters.expression;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::MIDI_CC_1) = &_midi_cc_1_src;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::MIDI_CC_2) = &_midi_cc_2_src;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::CV_1) = &_cv_a;
        _matrix_data._matrix_srcs.at((uint)ModMatrixSrc::CV_2) = &_cv_b;

        // add dsts
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::LFO_1_GAIN) = _lfo_1.getGainIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::LFO_1_RATE) = _lfo_1.getPitchIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::LFO_2_GAIN) = _lfo_2.getGainIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::LFO_2_RATE) = _lfo_2.getPitchIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::LFO_3_GAIN) = _lfo_3.getGainIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::LFO_3_RATE) = _lfo_3.getPitchIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FX_MACRO) = &_fx_macro_mod;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FX_SEND) = &_fx_send_mod;

        dst_is_fast.at((uint)ModMatrixDst::OSC_1_PITCH) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_1_PITCH) = &_voice_data.osc_1_pitch_tmp;
        dst_is_fast.at((uint)ModMatrixDst::OSC_1_SHAPE) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_1_SHAPE) = _osc_mixer_1.getShapeInput();
        dst_is_fast.at((uint)ModMatrixDst::OSC_2_PITCH) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_2_PITCH) = &_voice_data.osc_2_pitch_tmp;
        dst_is_fast.at((uint)ModMatrixDst::OSC_2_SHAPE) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_2_SHAPE) = _osc_mixer_2.getShapeInput();
        dst_is_fast.at((uint)ModMatrixDst::OSC_1_LEVEL) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_1_LEVEL) = _osc_mixer_1.getGainInput();
        dst_is_fast.at((uint)ModMatrixDst::OSC_2_LEVEL) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_2_LEVEL) = _osc_mixer_2.getGainInput();

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_CUTOFF_HP) = &_voice_data.filt_1_cut_tmp;
        dst_is_fast.at((uint)ModMatrixDst::FILTER_CUTOFF_HP) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_RESONANCE_HP) = &_voice_data.filt_1_res_tmp;
        dst_is_fast.at((uint)ModMatrixDst::FILTER_RESONANCE_HP) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_CUTOFF_LP) = &_voice_data.filt_2_cut_tmp;
        dst_is_fast.at((uint)ModMatrixDst::FILTER_CUTOFF_LP) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_RESONANCE_LP) = &_voice_data.filt_2_res_tmp;
        dst_is_fast.at((uint)ModMatrixDst::FILTER_RESONANCE_LP) = true;

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_3_PITCH) = &_voice_data.osc_3_pitch_tmp;
        dst_is_fast.at((uint)ModMatrixDst::OSC_3_PITCH) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_3_LEVEL) = &_voice_data.osc_3_lev_tmp;
        dst_is_fast.at((uint)ModMatrixDst::OSC_3_LEVEL) = true;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_3_SHAPE) = &_voice_data.osc_3_shape_tmp;
        dst_is_fast.at((uint)ModMatrixDst::OSC_3_SHAPE) = true;

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::DRIVE) = _drive_compensator.getdriveInput();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_4_TONE) = &_voice_data.osc_4_tone_tmp;
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::OSC_4_LEVEL) = &_voice_data.xor_lev_tmp;
        dst_is_fast.at((uint)ModMatrixDst::OSC_4_LEVEL) = true;

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AMP_EG_ATT) = _amp_env.getAttIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AMP_EG_DEC) = _amp_env.getDecIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AMP_EG_SUS) = _amp_env.getSusIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AMP_EG_REL) = _amp_env.getRelIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AMP_EG_LEVEL) = _amp_env.getGainIn();

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_EG_ATT) = _filt_env.getAttIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_EG_DEC) = _filt_env.getDecIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_EG_SUS) = _filt_env.getSusIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_EG_REL) = _filt_env.getRelIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::FILTER_EG_LEVEL) = _filt_env.getGainIn();

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AUX_EG_ATT) = _aux_env.getAttIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AUX_EG_DEC) = _aux_env.getDecIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AUX_EG_SUS) = _aux_env.getSusIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AUX_EG_REL) = _aux_env.getRelIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::AUX_EG_LEVEL) = _aux_env.getGainIn();

        _matrix_data._matrix_dests.at((uint)ModMatrixDst::PAN) = _panner.getPanIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::SPIN) = _panner.getSpinIn();
        _matrix_data._matrix_dests.at((uint)ModMatrixDst::MORPH) = &_morph_dest;

        // add the fast dsts to the fast dst list. by doign it this way, they are in the same order that they are in the dst/src array
        int f_c = 0;
        for (int i = 0; i < TOTAL_MOD_DSTS; i++) {
            if (dst_is_fast.at(i)) {
                _matrix_data._fast_dst.at(f_c) = _matrix_data._matrix_dests.at(i);
                f_c++;
            }
        }
        f_c = 0;
        for (int i = 0; i < TOTAL_MOD_SRC; i++) {
            if (src_is_fast.at(i)) {
                _matrix_data._fast_src.at(f_c) = _matrix_data._matrix_srcs.at(i);
                f_c++;
            }
        }
        // Setup the matrix gain pointers

        for (uint s = 0; s < (uint)ModMatrixSrc::NUM_SRCS; s++) {
            for (uint d = 0; d < (uint)ModMatrixDst::NUM_DSTS; d++) {
                _matrix_data._gains.at(s * (uint)ModMatrixDst::NUM_DSTS + d) = &_morphed_parameters.at(get_state_id_from_mod_matrix(ModMatrixSrc(s), ModMatrixDst(d)));
            }
        }
        auto f_gain_it = _matrix_data._fast_gains.begin();
        for (auto &f_source : _matrix_data._fast_src) {
            uint src_i = std::distance(_matrix_data._matrix_srcs.begin(), std::find(_matrix_data._matrix_srcs.begin(), _matrix_data._matrix_srcs.end(), f_source));
            for (auto &f_dst : _matrix_data._fast_dst) {
                uint dst_i = std::distance(_matrix_data._matrix_dests.begin(), std::find(_matrix_data._matrix_dests.begin(), _matrix_data._matrix_dests.end(), f_dst));
                *f_gain_it = _matrix_data._gains.at((src_i)*TOTAL_MOD_DSTS + (dst_i));
                _matrix_data._gains.at((src_i)*TOTAL_MOD_DSTS + (dst_i)) = &_zero;
                f_gain_it++;
            }
        }
    }

    std::array<float, TOTAL_MOD_DSTS> _drift_data;
    MoniqueMatrix _matrix = MoniqueMatrix(_drift_data);
    ModMatrixData &_matrix_data = _matrix.getData();
};

} // namespace Monique
