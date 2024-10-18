/**
 * @file MoniqueLayer.h
 * @brief Monique Layer class definitions.
 *
 * @copyright Copyright (c) 2022-2024 Melbourne Instruments, Australia
 */
#pragma once
#include "MoniqueCommon.h"
#include "MoniqueDriftSource.h"
#include "MoniqueModel.h"
#include "MoniqueVoice.h"
#include "common.h"

namespace Monique {
const std::map<int, std::vector<float>> unison_pan_positions = {
    { 1,                                                                                                        {0.f}},
    { 2,                                                                                                      {-1, 1}},
    { 3,                                                                                                 {-1, 0.f, 1}},
    { 4,                                                                                       {-1, -0.33f, 0.33f, 1}},
    { 5,                                                                                      {-1, -0.5f, 0, 0.5f, 1}},
    { 6,                                                                        {-1.f, -0.6f, -0.2f, 0.2f, 0.6f, 1.f}},
    { 7,                                                                         {-1, -0.66f, -.33, 0, 0.33f, .66, 1}},
    { 8,                                                                {-1, -0.75f, -.5, -0.25, 0.25f, 0.5, 0.75, 1}},
    { 9,                                                                {-1, -0.75f, -0.5, -0.25, 0, .25, .5, .75, 1}},
    {10,                                                                 {-1, -0.8, -.6, -.4, -.2, .2, .4, .6, .8, 1}},
    {11,                                                             {-1, -0.8, -0.6, -.4, -.2, 0, .2, .4, .6, .8, 1}},
    {12, {-1.f, -0.8181f, -0.6363f, -0.45454f, -0.2787f, -0.0909f, 0.0909f, 0.2787f, 0.45454f, 0.6363f, 0.8181f, 1.f}}
};

class MoniqueLayer {

  public:
    MoniqueLayer(int layer_num, std::vector<ParamChange> &param_changes, std::array<int, NUM_VOICES> &num_voices_array, GlobalSynthParameters &gp) :
        _layer_number(layer_num),
        _param_changes(param_changes),
        _num_voices_array(num_voices_array),
        _gp(gp) {
#ifdef EN_DEBUG
        // printf("\n init layer %d ", layer_num);
        evaluateLayerVoices();
#endif
    }

    bool midiNoteFilter(const int channel, const int note_n) {
        if (note_n >= _layer_low_split && note_n <= _layer_high_split) {
            if (!_omni_mode && (channel == _layer_channel_filter)) {
                return true;
            } else if (_omni_mode) {
                return true;
            } else {
                return false;
            }
        }
        return false;
    }

    void evaluateLayerVoices() {
        int new_first_voice;
        int new_last_voice;
        if (_layer_number == 0) {
            new_first_voice = 0;
            new_last_voice = (_num_voices_array[(0)]);
            new_last_voice = std::min((int)NUM_VOICES, new_last_voice);
            _num_voices = new_last_voice;
        }
        if (_layer_number == 1) {
            new_first_voice = _num_voices_array[(0)];
            new_last_voice = std::max(0, _num_voices_array[(1)] + new_first_voice);
            new_last_voice = std::min((int)NUM_VOICES, new_last_voice);
            _num_voices = new_last_voice - new_first_voice;
        }
        _first_voice = new_first_voice;
        _last_voice = new_last_voice;

        // mute voices which are not allocated
        for (int i = 0; i < NUM_VOICES; i++) {
            if (i < new_first_voice || i > new_last_voice) {
                _voice_array[(i)].setUnalloc(true);
                _voice_array[(i)].getVoiceModel().setMute(true);
            } else {

                _voice_array[(i)].setUnalloc(false);
                _voice_array[(i)].getVoiceModel().setMute(false);
            }
        }
    };

    void allNotesOff() {

        // send a note off to every voice
        MidiNote note;
        note.pitch = 0;
        for (auto &voice : _voice_array) {
            last_released_voice = voice.setNoteOff(note);
            voice.forceNoteOff();
        }
        // clear mono note stacks
        _held_notes.clear();
    }

    bool findUnisonVoiceFree(int &pitch, int &unison_voice_found) {
        int distance = 127;
        bool found = false;
        switch (_parameters.para_mode) {
        case ParaphonicModes::TIME_NO_RETRIG:
        case ParaphonicModes::TIME_RETRIG: {
            for (int i = _first_voice; i < _last_voice; i++) {
                auto &voice = _voice_array.at(i);

                if (voice.realVoiceLastAlloc() && (!voice.unisonGateOn()) && (voice.isAtt() || (voice.buffersSinceTrigger() < SHORT_TRIG_TIME_THRESH) && voice.mainGateOn())) {

                    // find the closest note which is also within the time window
                    if (distance > (abs(pitch - (int)voice.currentNote().pitch))) {
                        unison_voice_found = i;
                        distance = abs(pitch - (int)voice.currentNote().pitch);
                    }
                    found = true;
                }
            }
        } break;
        case ParaphonicModes::DIST_NO_RETRIG:
        case ParaphonicModes::DIST_RETRIG: {
            for (int i = _first_voice; i < _last_voice; i++) {
                auto &voice = _voice_array.at(i);
                if (voice.realVoiceLastAlloc() && (!voice.unisonGateOn()) && !voice.probablyOff()) {

                    // find the closest note which is also within the time window
                    const int this_dist = abs(pitch - (int)voice.currentNote().pitch);
                    if (distance > (this_dist) && this_dist < PARAPHONIC_MAX_DIST) {
                        unison_voice_found = i;
                        distance = this_dist;
                        found = true;
                    }
                }
            }
        } break;

        default:
            break;
        }
        return found;
    }

    bool findUnisonVoice2nd(int &pitch, int &unison_voice_found) {
        int distance = 127;
        int voice_num = 0;
        bool found = false;
        switch (_parameters.para_mode) {
        case ParaphonicModes::TIME_NO_RETRIG:
        case ParaphonicModes::TIME_RETRIG: {
            int min_time = SHORT_TRIG_TIME_THRESH * 2;
            for (int i = _first_voice; i < _last_voice; i++) {
                auto &voice = _voice_array.at(i);
                const int d = abs(pitch - (int)voice.currentNote().pitch);
                const bool unison_is_free = !voice.unisonGateOn();
                int time = voice.buffersSinceTrigger();
                if ((time < 2 * SHORT_TRIG_TIME_THRESH) && voice.realVoiceLastAlloc()) {
                    if (time < min_time) {
                        min_time = time;
                    }
                    found = true;
                    unison_voice_found = i;
                }
            }
        } break;
        case ParaphonicModes::DIST_NO_RETRIG:
        case ParaphonicModes::DIST_RETRIG: {

            for (int i = _first_voice; i < _last_voice; i++) {
                auto &voice = _voice_array.at(i);
                const int d = abs(pitch - (int)voice.currentNote().pitch);
                const bool unison_is_free = !voice.unisonGateOn();
                if (d < distance && unison_is_free && voice.realVoiceLastAlloc()) {
                    distance = d;
                    voice_num = i;
                }
                if (distance < PARAPHONIC_MAX_DIST) {
                    found = true;
                    unison_voice_found = voice_num;
                }
            }
        } break;
        default:
            break;
        }
        return found;
    }

    void paraModeSetup() {

        for (auto &voice : _voice_array) {
            voice.paraModeSetup();
        }
    };

    volatile int fake_voice = 0;

    inline void handleNote(MidiNote note) {
        // check if the layer has voices
        if (_num_voices == 0) {
            return;
        }
        // check if the midi note is for this layer
        // printf(" l %d  note: %d %d \n", _layer_number, note.pitch, 0);
        if (!midiNoteFilter(note.channel, note.pitch)) {
            // printf("\n %d reject", _layer_number);
            return;
        }
        // printf(" %d ", _layer_number);

        if (_voice_mode == VoiceMode::POLY) {
            // printf("allocate poly ");
            if (note.note_type == MidiNote::type::ON) {

                int num_allocated = 0;
                if (_round_robin_voice >= _num_voices) {
                    _round_robin_voice = 0;
                }

                for (int i = _first_voice; i <= _last_voice; i++) {
                    if (num_allocated >= _num_unison) {
                        break;
                    }
                    // printf("\n check voice %d ", i);
                    //  Get the next voice to allocate
                    //  Note - wrap around if needed
                    int v = i + _round_robin_voice;
                    if (v >= _last_voice) {
                        v -= _num_voices;
                    }
                    int unison_voice;
                    if (findUnisonVoiceFree(note.pitch, unison_voice)) {
                        _voice_array.at(unison_voice).setUnisonNote(note, generate_detune(num_allocated));
                        ++num_allocated;
                        if (num_allocated >= _num_unison) {
                            break;
                        }
                    }
                    // Allocate this voice if it is  released
                    if (!_voice_array[v].blacklisted() && _voice_array[v].inReleased()) {
                        // printf("released ");

                        // allocate the voice
                        _last_allocated = _voice_array.at(v).setNoteOn(note, generate_pan(num_allocated), generate_detune(num_allocated));

                        // If we have allocated enough voices?
                        if (++num_allocated >= _num_unison) {
                            // Yes, set the round robin voice to be used in the next
                            // allocation
                            // Note - wrap around to the first voice if needed
                            _round_robin_voice = v + 1 - _first_voice;
                            if (_round_robin_voice >= _num_voices) {
                                _round_robin_voice = 0;
                            }
                            break;
                        }
                    }
                }
                if (_round_robin_voice >= _num_voices) {
                    _round_robin_voice = 0;
                }
                // here we try to allocate unison voices which meet some requirements according to the current paraphonic mode
                while ((num_allocated < _num_unison) && (num_allocated < _num_voices)) {
                    int unison_voice;

                    bool found = findUnisonVoice2nd(note.pitch, unison_voice);
                    if (!found) {
                        break;
                    } else {
                        _voice_array.at(unison_voice).setUnisonNote(note, generate_detune(num_allocated));
                        ++num_allocated;
                    }
                }
                if (_round_robin_voice >= _num_voices) {
                    _round_robin_voice = 0;
                }
                // printf("\n rr failed, %d ", num_allocated);

                //  Have we allocated enough voices AND there are other non-free voices available, here we steal voices with pending releases
                if ((num_allocated < _num_unison) && (num_allocated < _num_voices)) {
                    // Steal the oldest allocated voices until we have enough
                    for (int i = _first_voice; i <= _last_voice; i++) {
                        // printf("\n try steal %d ", i);
                        //  Get the next voice to allocate
                        //  Note - wrap around if needed
                        int v = i + _round_robin_voice;
                        if (v >= _last_voice) {
                            v -= _num_voices;
                        }
                        // Steal this voice if it is currently being used and there is a pending release
                        if (!_voice_array.at(v).blacklisted() && !_voice_array.at(v).inReleased() && _voice_array.at(v).pendingRelease()) {
                            // printf(" steal this ");
                            //  Allocate this voice
                            _last_allocated = _voice_array.at(v).setNoteOn(note, generate_pan(num_allocated), generate_detune(num_allocated));
                            // If we have allocated enough voices?
                            if (++num_allocated >= _num_unison) {
                                // Yes, set the round robin voice to be used in the next
                                // allocation
                                // Note - wrap around to the first voice if needed
                                _round_robin_voice = v + 1 - _first_voice;
                                if (_round_robin_voice >= _num_voices) {
                                    _round_robin_voice = 0;
                                }
                                break;
                            }
                        }
                    }
                }
                if (_round_robin_voice >= _num_voices) {
                    _round_robin_voice = 0;
                }

                // we still dont have enough voices, so now steal voices without pending releases
                if ((num_allocated < _num_unison) && (num_allocated < _num_voices)) {
                    // Steal the oldest allocated voices until we have enough
                    for (int i = _first_voice; i <= _last_voice; i++) {

                        // printf("\n try last %d ", i);
                        //  Get the next voice to allocate
                        //  Note - wrap around if needed
                        int v = i + _round_robin_voice;
                        if (v >= _last_voice) {
                            v -= _num_voices;
                        }

                        // Steal this voice if it is currently being used
                        if (!_voice_array[v].blacklisted() && !_voice_array[v].inReleased()) {
                            // Allocate this voice
                            _last_allocated = _voice_array.at(v).setNoteOn(note, generate_pan(num_allocated), generate_detune(num_allocated));
                            _voice_array.at(v).setUnisonOff();
                            // printf(" steal this");
                            //  If we have allocated enough voices?
                            if (++num_allocated >= _num_unison) {
                                // Yes, set the round robin voice to be used in the next
                                // allocation
                                // Note - wrap around to the first voice if needed
                                _round_robin_voice = v + 1 - _first_voice;
                                if (_round_robin_voice >= _num_voices) {
                                    _round_robin_voice = 0;
                                }
                                break;
                            }
                        }
                    }
                }
            }

            // handle poly note offs
            else {
                // Free each active voice with this note value
                for (auto &voice : _voice_array) {
                    // Get the voice MIDI note (if any), and check if it is the same as the
                    // passed note (pitch)
                    if ((voice.currentNote().pitch == note.pitch)) {
                        // Free this voice
                        last_released_voice = voice.setNoteOff(note);
                    }
                    if (_parameters.para_mode != ParaphonicModes::DISABLE) {
                        if (voice.currentUnison().pitch == note.pitch) {
                            voice.setUnisonOff(note);
                        }
                    }
                }
            }
        }

        // mono voice allocation
        else {
            if (note.note_type == MidiNote::type::ON) {
                // printf(" mono alloc  ");
                //  Keep track of the currently held mono notes so we can swap to held notes on
                //  release of the current note
                _held_notes.push_back(note);
                if (_voice_mode == VoiceMode::AUTO_LOW) {
                    // printf(" auto  ");
                    int lowest_note = 127;
                    for (auto &i_note : _held_notes) {
                        lowest_note = (int)i_note.pitch < lowest_note ? i_note.pitch : lowest_note;
                    }
                    if (lowest_note < note.pitch) {
                        return;
                    }
                }

                // Allocate the first unison number of voices
                int allocated = 0;
                int checked = 0;
                for (int i = _first_voice; i <= _last_voice; i++) {
                    // Is this voice not blacklisted?

                    int unison_voice;
                    if (findUnisonVoiceFree(note.pitch, unison_voice)) {
                    }
                    if (!_voice_array[i].blacklisted()) {
                        // Add this voice (always steal voices in unison mode)
                        allocated++;
                        float detune = generate_detune(allocated);
                        float pan = generate_pan(allocated);
                        _last_allocated = _voice_array.at(i).setNoteOn(note, pan, detune);
                    }

                    // If we have checked all possible voices, or we have allocated
                    // enough voices, exit the allocation loop
                    // Note this may be less than unison if there are not enough voices to allocate
                    if (++checked >= _num_voices) {
                        break;
                    }
                    if (allocated >= _num_unison) {
                        break;
                    }

                    if (_parameters.para_mode != ParaphonicModes::DISABLE) {
                        allocated++;
                        float detune = generate_detune(allocated);
                        _voice_array.at(i).setUnisonNote(note, (detune));
                        if (allocated >= _num_unison) {
                            break;
                        }
                    }
                }
            }
            // mono note off case
            else {
                // Legato or Mono-retrigger mode
                // Find all held notes with the same pitch and remove them
                for (auto it = _held_notes.begin(); it != _held_notes.end();) {
                    if (note.pitch == it->pitch) {
                        it = _held_notes.erase(it);
                    } else {
                        ++it;
                    }
                }
                // If there are no held notes remaining then all voices should be set off
                if (_held_notes.size() == 0) {
                    for (auto &voice : _voice_array) {
                        last_released_voice = voice.setNoteOff(note);
                        voice.setUnisonOff(note);
                    }
                } else {
                    // default: we just get the current newest note
                    MidiNote note = _held_notes.back();
                    if (_voice_mode == VoiceMode::AUTO_LOW) {
                        // find the lowest held note if we are in auto low mode
                        for (auto &n : _held_notes) {
                            if (n.pitch < note.pitch) {
                                note = n;
                            }
                        }
                    }
                    // There are held notes remaining, process the newest
                    for (MoniqueVoice &voice : _voice_array) {
                        if (!voice.inReleased()) {
                            if (voice.currentNote().pitch != note.pitch) {
                                _last_allocated = voice.setNoteOn(note);
                            }
                        }
                        if (voice.unisonGateOn()) {
                            if (voice.currentUnison().pitch != note.pitch) {
                                voice.setUnisonNote(note);
                            }
                        }
                    }
                }
            }
        }

        if (_parameters.all_notes_off && (note.note_type == MidiNote::ON)) {
            _fx_mod_voice_sel = _last_allocated;
        }

        // printf("\n");
    }

    inline float
    generate_detune(int number_allocated) {

        const int lookup = number_allocated % _num_unison;
        const std::vector<float> &detune = unison_pan_positions.at(_num_unison);
        return detune.at(lookup);
    }

    float generate_pan(int number_allocated) {
        number_allocated = _pan_counter++;
        if (_pan_counter >= _pan_number) {
            _pan_counter = 0;
        }
        float voice_pan = 0;

        // if pan number is zero, we should just return no pan
        if (_pan_number == 0 || _pan_number == 1) {
            return voice_pan;
        }
        if (_pan_mode == PanModes::SPREAD) {
            float span_range = 2.0f / (float)((_pan_number - 1));
            voice_pan = ((float)number_allocated * span_range) - 1.f;
            return voice_pan;
        }
        if (_pan_mode == PanModes::PINGPONG) {
            // ping pong odd case
            if (_pan_number % 2 != 0) {
                const float span_range = 2.0f / (float)((_pan_number - 1));
                const float mult = _pan_phase ? 0 + (float)((number_allocated + 2) / 2) : 0 - (float)((number_allocated + 1) / 2);
                _pan_phase = !_pan_phase;
                // printf("\n mult %f", mult);

                voice_pan = (mult * span_range);
                return voice_pan;
            }
            // ping pong even case
            else {
                const float span_range = 1.0f / (float)((_pan_number - 1));

                const float mult = _pan_phase ? 0 - (float)(number_allocated) : 0 + (float)(number_allocated + 1);
                // printf("\n %f %d", mult, _pan_phase);
                _pan_phase = !_pan_phase;
                voice_pan = (mult * span_range);
                return voice_pan;
            }
        }
        return voice_pan;
    }

    inline void polyPressure(PolyPressure &pp) {
        const int &pitch = pp.note_pitch;
        for (int i = _first_voice; i < _last_voice; i++) {
            if (_voice_array.at(i).currentNote().pitch == pitch) {
                _voice_array.at(i).setPolyAT(pp);
            }
        }
    }

    inline void reloadCal() {
        for (auto &voice : _voice_array) {
            voice.reloadCal();
        }
    }

    inline void updateParameter(const ParamDecoded &param, const float value) {
        float denormalised;
        // ParamNameGenerator name;
        //  printf("\n %d %s %d  %d %f", param_num, name.get_parameter_name(get_id_from_layer_param(param.layer_param)).c_str(), param.layer_num, param_num, value);
        switch (param.p_type) {
        case ParamType::GLOBAL:
            // denormlise the global parameter, then handle the specific cases;
            switch (param.global_param) {
            case GlobalParams::MASTER_TUNE:
                denormalised = state_param_from_normalised_float(GlobalParams::MASTER_TUNE, value);
                _master_detune_semitone = (denormalised - 0.5) / (3 * NOTE_GAIN);
                break;
            case GlobalParams::SELECTED_LAYER:
                _layer_select = std::round(value);
                _current_layer = _layer_select == _layer_number;
                break;

            default:
#ifdef EN_DEBUG
                // printf("\n unhandled global")
#endif
                break;
            }

            break;

        case ParamType::PRESET_COMMON:
            switch (param.preset_common_param) {
            case PresetCommonParameters::LAYER_1_NUM_VOICES: {
                denormalised = preset_common_from_normalised_float(PresetCommonParameters::LAYER_1_NUM_VOICES, value);
                if (true) {
                    _num_voices_array.at(0) = (int)denormalised;
                    evaluateLayerVoices();
                    allNotesOff();
                }

            } break;

            case PresetCommonParameters::LAYER_2_NUM_VOICES: {
                denormalised = preset_common_from_normalised_float(PresetCommonParameters::LAYER_2_NUM_VOICES, value);
                if (true) {
                    _num_voices_array.at(1) = (int)denormalised;
                    evaluateLayerVoices();
                    allNotesOff();
                }

            } break;
            case get_preset_common_id_from_mod_matrix((ModMatrixSrc)0, PresetCommonModMatrixDst::FX_SEND_LEVEL)... get_preset_common_id_from_mod_matrix(ModMatrixSrc::CONSTANT, PresetCommonModMatrixDst::FX_MACRO_LEVEL): {
                // only run this processing for layer 1
                if (_layer_number == 0) {
                    if (param.preset_common_param == get_preset_common_id_from_mod_matrix(ModMatrixSrc::CONSTANT, PresetCommonModMatrixDst::FX_MACRO_LEVEL)) {

                        int dst = (int)get_preset_common_dst_from_state_id(param.preset_common_param) + (int)ModMatrixDst::FX_SEND;
                        ModMatrixSrc src = get_preset_common_src_from_state_id(param.preset_common_param);
                        int id = get_state_id_from_mod_matrix(src, (ModMatrixDst)dst);
                        _parameters.state_a.at(id) = 0;
                        _parameters.state_b.at(id) = 0;
                    } else {
                        int dst = (int)get_preset_common_dst_from_state_id(param.preset_common_param) + (int)ModMatrixDst::FX_SEND;
                        ModMatrixSrc src = get_preset_common_src_from_state_id(param.preset_common_param);
                        int id = get_state_id_from_mod_matrix(src, (ModMatrixDst)dst);
                        _parameters.state_a.at(id) = value;
                        _parameters.state_b.at(id) = value;
                    }
                }
            };
                break;
            }
            break;

        case ParamType::LAYER:
            // skip parameter if its not for this layer
            if (param.layer_num != _layer_number) {
                return;
            }

            // denormlise the global parameter, then handle the specific cases;

            switch (param.layer_param) {
            case LayerParameters::OCTAVE_OFFSET:
                denormalised = from_normalised_float(LayerParameters::OCTAVE_OFFSET, value);
                _parameters.octave_offset = denormalised * 12.f;
                break;
            case LayerParameters::OUTPUT_ROUTING: {
                denormalised = from_normalised_float(LayerParameters::OUTPUT_ROUTING, value);
                OutputRouting mode = (OutputRouting)(int)denormalised;
                _parameters.output_mode = mode;
            } break;
            case LayerParameters::CHANNEL_FILTER: {
                denormalised = from_normalised_float(LayerParameters::CHANNEL_FILTER, value);
                int16_t setting = std::round(denormalised);
                // 0 means omni
                if (setting == 0) {
                    _omni_mode = true;
                    _layer_channel_filter = 0;
                } else {
                    _omni_mode = false;
                    if (true) {
                        allNotesOff();
                        _layer_channel_filter = setting - 1;
                    }
                }
                allNotesOff();
                // printf("\n set channel filter %d %d", _layer_number, _layer_channel_filter);
            } break;
                break;

            case LayerParameters::SPLIT_LOW_NOTE: {
                denormalised = from_normalised_float(LayerParameters::SPLIT_LOW_NOTE, value);
                int16_t setting = (int)(denormalised);
                if (true) {
                    _layer_low_split = setting;
                    allNotesOff();
                }

                // printf("\nlow %d %d", _layer_number, _layer_low_split);
            } break;
            case LayerParameters::EXT_IN_GAIN:
                _parameters.ext_in_gain = EXT_IN_GAIN_MULT * value;
                break;
            case LayerParameters::SPLIT_HIGH_NOTE: {
                denormalised = from_normalised_float(LayerParameters::SPLIT_HIGH_NOTE, value);
                int16_t setting = (int)(denormalised);
                if (true) {
                    _layer_high_split = setting;
                    allNotesOff();
                }
                // printf("\nhigh %d %d", _layer_number, _layer_high_split);
            } break;

            case LayerParameters::LAYER_VOLUME: {
                denormalised = from_normalised_float(LayerParameters::LAYER_VOLUME, value);

                _parameters.layer_volume = denormalised * denormalised;
                // printf("\nhigh %d %f", _layer_number, _parameters.layer_volume);
            } break;

            default:
#ifdef EN_DEBUG
                // printf("\n unhandled layer");
#endif
                break;
            }
            break;
        case ParamType::COMMON:

#ifdef EN_DEBUG
            // printf("\n common  %d %d %d %f", param_num, param.common_param, param.layer_num, value);
            // printf(" %s ", name.get_parameter_name(get_id_from_common_param(param.common_param)).c_str());
#endif
            // skip parameter if its not for this layer
            if (param.layer_num != _layer_number) {
                return;
            }

            // denormlise the global parameter, then handle the specific cases;

            switch (param.common_param) {
            case CommonParameters::GLIDE_MODE: {
                denormalised = from_normalised_float(CommonParameters::GLIDE_MODE, value);
                GlideModes mode = (GlideModes)(denormalised);
#ifdef EN_DEBUG
                // assert(mode < GlideModes::NUM_GLIDE_MODES);
#endif
                _glide_mode = mode;
            } break;

            case CommonParameters::POLY_MODE: {
                denormalised = from_normalised_float(CommonParameters::POLY_MODE, value);
                _voice_mode = (VoiceMode)(int)denormalised;
                // printf("\n vmode %d", (int)_voice_mode);
                if (_voice_mode >= VoiceMode::NUM_VOICE_MODES) {
                    // printf("\n mode %d", (int)_voice_mode);
                    _voice_mode = VoiceMode::POLY;
                }

#ifdef EN_DEBUG
                // assert(_voice_mode < VoiceMode::NUM_VOICE_MODES);
#endif
            } break;
            case CommonParameters::MORPH: {
                _parameters.common_params.at(CommonParameters::MORPH) = value;
                _refresh_morph = true;
                if (_current_layer) {
                    if (((value < 1) && (_layer_state == State::B)) ||
                        ((value > 0) && (_layer_state == State::A))) {
                        _currently_morphing = true;
                        _param_changes.insert(_param_changes.begin(), ParamChange(gen_param_id(GlobalParams::MORPHING), 1, true));
                    } else {
                        _exit_morph = _currently_morphing;
                        _currently_morphing = false;
                        _param_changes.insert(_param_changes.begin(), ParamChange(gen_param_id(GlobalParams::MORPHING), 0, true));
                    }
                }
            } break;
            case CommonParameters::STATE: {
                const State state = (State)std::round(value);
                // printf("\n set state %d", (int)state);
                _parameters.state_morphed_blocked.fill(false);
                _layer_state = state;
                if (_morph_mode_dance) {
                    _refresh_morph = true;
                    _parameters.common_params.at(CommonParameters::MORPH) = (float)state;
                    // printf(" morph %f", _parameters.common_params.at(CommonParameters::MORPH));
                    _exit_morph = true;
                    _currently_morphing = false;
                    _param_changes.insert(_param_changes.begin(), ParamChange(gen_param_id(GlobalParams::MORPHING), 0, true));
                }
            } break;
            case CommonParameters::WT_INTERPOLATE_MODE:
                denormalised = from_normalised_float(CommonParameters::WT_INTERPOLATE_MODE, value);
                _parameters.wt_interpolate = (int)denormalised;
                // printf("\n %d %f", _parameters.wt_interpolate, denormalised);
                break;
            case CommonParameters::VINTAGE_AMT:
                denormalised = from_normalised_float(CommonParameters::VINTAGE_AMT, value);
                _parameters.layer_vintage = VINTAGE_SCALE * denormalised;
                _parameters.layer_drift = DRIFT_SCALE * (denormalised * denormalised + DRIFT_OFFSET);
                break;
            case CommonParameters::PARAPHONIC_MODE:
                denormalised = from_normalised_float(CommonParameters::PARAPHONIC_MODE, value);
                _parameters.para_mode = (ParaphonicModes)(int)denormalised;
                paraModeSetup();
                break;
            case CommonParameters::LFO_1_GLOBAL:

                denormalised = from_normalised_float(CommonParameters::LFO_1_GLOBAL, value);
                _parameters.lfo_1_global = (int)denormalised;
                break;
            case CommonParameters::LFO_2_GLOBAL:

                denormalised = from_normalised_float(CommonParameters::LFO_2_GLOBAL, value);
                _parameters.lfo_2_global = (int)denormalised;
                break;
            case CommonParameters::LFO_3_GLOBAL:

                denormalised = from_normalised_float(CommonParameters::LFO_3_GLOBAL, value);
                _parameters.lfo_3_global = (int)denormalised;
                break;
            case CommonParameters::UNISON_VOICES:
                denormalised = from_normalised_float(CommonParameters::UNISON_VOICES, value);
                _num_unison = (int)denormalised;
                break;
            case CommonParameters::PAN_MODE:
                denormalised = from_normalised_float(CommonParameters::PAN_MODE, value);
                _pan_mode = PanModes::PINGPONG;
                break;
            case CommonParameters::PAN_NUM:
                denormalised = from_normalised_float(CommonParameters::PAN_NUM, value);
                _pan_number = (int)denormalised;
                break;
            case CommonParameters::PITCH_BEND_RANGE:
                denormalised = from_normalised_float(CommonParameters::PITCH_BEND_RANGE, value);
                _pb_range = denormalised;
                break;
            case CommonParameters::MIDI_MODWHEEL:
                denormalised = from_normalised_float(CommonParameters::MIDI_MODWHEEL, value);
                _parameters.common_params.at((int)CommonParameters::MIDI_MODWHEEL) = denormalised;
                break;
            case CommonParameters::MIDI_SUSTAIN:
                denormalised = from_normalised_float(CommonParameters::MIDI_SUSTAIN, value);
                _parameters.sustain_ped = denormalised > PARAM_BOOL_COMP;
                break;
            case CommonParameters::MIDI_PITCH_BEND:
                denormalised = from_normalised_float(CommonParameters::MIDI_PITCH_BEND, value);
                _pb_value = denormalised;
                break;
            case CommonParameters::MIDI_AFTERTOUCH:
                denormalised = from_normalised_float(CommonParameters::MIDI_AFTERTOUCH, value);
                _at_value = denormalised;
                break;
            case CommonParameters::MIDI_EXPRESSION:
                denormalised = from_normalised_float(CommonParameters::MIDI_EXPRESSION, value);
                _expression_value = denormalised;
                break;
            case CommonParameters::MIDI_CC_1_MOD_SOURCE:
                denormalised = from_normalised_float(CommonParameters::MIDI_CC_1_MOD_SOURCE, value);
                _parameters.common_params.at(CommonParameters::MIDI_CC_1_MOD_SOURCE) = value;
                break;
            case CommonParameters::MIDI_CC_2_MOD_SOURCE:
                denormalised = from_normalised_float(CommonParameters::MIDI_CC_2_MOD_SOURCE, value);
                _parameters.common_params.at(CommonParameters::MIDI_CC_2_MOD_SOURCE) = value;
                break;
            case CommonParameters::LFO_1_TEMPO_SYNC:
                denormalised = from_normalised_float(CommonParameters::LFO_1_TEMPO_SYNC, value);
                _parameters.lfo_1_sync = denormalised;
                break;
            case CommonParameters::LFO_2_TEMPO_SYNC:
                denormalised = from_normalised_float(CommonParameters::LFO_2_TEMPO_SYNC, value);
                _parameters.lfo_2_sync = denormalised;
                break;
            case CommonParameters::LFO_3_TEMPO_SYNC:
                denormalised = from_normalised_float(CommonParameters::LFO_3_TEMPO_SYNC, value);
                _parameters.lfo_3_sync = denormalised;
                break;
            default:
#ifdef EN_DEBUG
                // printf("\n unhandled common");
#endif
                break;
            }

            break;

        case ParamType::STATE: {
            if (param.layer_num != _layer_number) {
                return;
            }

#ifdef EN_DEBUG
#endif
            // printf("\n state %d %s %d %d %f", (int)param.state_param, name.get_parameter_name(get_id_from_state_param(param.state_param)).c_str(), (int)param.layer_num, (int)param.state, value);
            if (fabs(value - 0.5) > 0.01) {
                //     printf(" Here ");
            }
            // apply change to corresponding state
            switch (param.state) {
            case State::A:
                switch (param.state_param) {
                case StateParameters::WT_SELECT:
                    _wave_loader_a.loadWavetable(value);
                    /* code */
                    break;

                default:
                    _parameters.state_a.at((int)param.state_param) = value;
                    break;
                }
                break;

            case State::B:
                switch (param.state_param) {
                case StateParameters::WT_SELECT:
                    _wave_loader_b.loadWavetable(value);
                    /* code */
                    break;
                default:
                    _parameters.state_b.at((int)param.state_param) = value;
                    break;
                }
                break;

            // should never hit this case
            default:
#ifdef EN_DEBUG
                // assert(false);
#endif
                break;
            }
            if (_currently_morphing) {
                if (!(_parameters.state_morphed_blocked.at((int)param.state_param))) {

                    _parameters.state_a_smooth.at((int)param.state_param) = _parameters.state_a.at((int)param.state_param);
                    _parameters.state_b_smooth.at((int)param.state_param) = _parameters.state_b.at((int)param.state_param);
                }
                _parameters.state_morphed_blocked.at((int)param.state_param) = true;
            }
        } break;

        default:
#ifdef EN_DEBUG
            // printf("\n unhandled param");
#endif
            break;
        }
    }

    // run after updating all parameters for this buffer
    inline void paramRefresh() {

        // if the morph value is updated then we recalculate morph value and reflect back values if changed
        _current_layer = _layer_select == _layer_number;
        const float morph = _parameters.common_params.at(CommonParameters::MORPH);
        _parameters.morph_value = morph;
        if (_refresh_morph && _current_layer && _morph_mode_dance) {
            const State state = _layer_state;
            _refresh_morph = false;

            // if the morph value is non-zero then we are morphing and should reflect params
            if (_currently_morphing || _exit_morph) {
                _exit_morph = false;

                const float inv_morph = 1 - morph;
                // printf("\n refresh morph %f", morph);
                //  std::array<float, (uint)StateParameters::NumStateParams> morphed_params_return;
                for (uint i = 0; i < TOTAL_STATE_PARAMS; ++i) {
                    const float &a = _parameters.state_a.at(i);
                    const float &b = _parameters.state_b.at(i);
                    if (a != b) {
                        if (!_parameters.state_morphed_blocked.at(i)) {
                            // perform weighted sum morph
                            float tmp = a * inv_morph + b * morph;
                            // push the change to the vector for sending to the UI
                            _param_changes.push_back(ParamChange(gen_param_id((StateParameters)i, state, _layer_number), tmp, true));
                        } else {
                        }
                    }
                }
                // refresh LFO parameters
            }
        }
        _parameters.lfo_1_a_shape = (LfoOscShape)(int)state_param_from_normalised_float((int)StateParameters::LFO_1_SHAPE, _parameters.state_a.at((int)StateParameters::LFO_1_SHAPE));
        _parameters.lfo_1_b_shape = (LfoOscShape)(int)state_param_from_normalised_float((int)StateParameters::LFO_1_SHAPE, _parameters.state_b.at((int)StateParameters::LFO_1_SHAPE));

        _parameters.lfo_2_a_shape = (LfoOscShape)(int)state_param_from_normalised_float((int)StateParameters::LFO_2_SHAPE, _parameters.state_a.at((int)StateParameters::LFO_2_SHAPE));
        _parameters.lfo_2_b_shape = (LfoOscShape)(int)state_param_from_normalised_float((int)StateParameters::LFO_2_SHAPE, _parameters.state_b.at((int)StateParameters::LFO_2_SHAPE));

        _parameters.lfo_3_a_shape = (LfoOscShape)(int)state_param_from_normalised_float((int)StateParameters::LFO_3_SHAPE, _parameters.state_a.at((int)StateParameters::LFO_3_SHAPE));
        _parameters.lfo_3_b_shape = (LfoOscShape)(int)state_param_from_normalised_float((int)StateParameters::LFO_3_SHAPE, _parameters.state_b.at((int)StateParameters::LFO_3_SHAPE));
    }

    void checkMutes() {
    }

    void setTempo(float tempo) {
        _parameters.tempo = tempo;
    }

    void getFxMod(float &snd_mod, float &macro_mod) {
        _voice_array.at(_fx_mod_voice_sel).getFxMod(snd_mod, macro_mod);
        // printf("\n %d %f %F", _last_allocated, snd_mod, macro_mod);
    }

    void runSynth(ProcessAudioData &audio_data) {
        // midi smoothing and calcs
        const float current_pitchbend = _pb_range * _pb_value / NOTE_GAIN / 12.f;
        param_smooth(current_pitchbend, _parameters.pitch_bend);
        param_smooth(_at_value, _parameters.aftertouch);
        _parameters.expression = mpeParamSmooth(_expression_value, _parameters.expression);
        _parameters.morph_value_smooth = param_smooth(_parameters.morph_value, _parameters.morph_value_smooth);

        for (int i = 0; i < TOTAL_STATE_PARAMS; i++) {
            int idx = -1;
            if (!(i == idx)) {
                param_smooth(_parameters.state_a[i], _parameters.state_a_smooth[i]);
                param_smooth(_parameters.state_b[i], _parameters.state_b_smooth[i]);
            } else {
                _parameters.state_a_smooth[i] = _parameters.state_a[i];
                _parameters.state_b_smooth[i] = _parameters.state_b[i];
            }
        }
        for (int i = 0; i < (int)ModMatrixSrc::NUM_SRCS; i++) {
            _parameters.state_b_smooth.at(get_state_id_from_mod_matrix((ModMatrixSrc)i, ModMatrixDst::MORPH)) = _parameters.state_a_smooth.at(get_state_id_from_mod_matrix((ModMatrixSrc)i, ModMatrixDst::MORPH));
        }
        _parameters.all_notes_off = true;
        _drift_gen.recalculate();
        for (int v = _first_voice; v < _last_voice; v++) {
            _parameters.all_notes_off = !_voice_array.at(v).isGateOn() && _parameters.all_notes_off;
            if (v == last_released_voice) {
                _voice_array.at(v).setLastReleased(true);
            } else {
                _voice_array.at(v).setLastReleased(false);
            }
            _voice_array.at(v).run();
            _voice_array.at(v).getVoiceModel().generateVoiceBuffers(audio_data);
        }
    }

  private:
    int _last_allocated = 0;
    int _fx_mod_voice_sel = 0;
    GlobalSynthParameters &_gp;
    const int _layer_number;
    int _layer_select = 0;
    DriftSource _drift_gen;

    int16_t _layer_channel_filter = 0;
    int16_t _layer_low_split = 0;
    int16_t _layer_high_split = 127;
    int _num_voices = 3;
    std::vector<MidiNote> _held_notes;
    int _first_voice = 0;
    int _last_voice = _num_voices;
    int _num_unison = 1;
    int _round_robin_voice = 0;
    bool _omni_mode = true;
    bool _refresh_morph = false;
    bool _morph_mode_dance = true;
    bool _current_layer = false;
    bool _currently_morphing = false;
    bool _exit_morph = false;
    std::vector<ParamChange> &_param_changes;

    // panner
    PanModes _pan_mode = PanModes::OFF;
    int _pan_number = 0;
    int _pan_counter = 0;
    bool _pan_phase = false;

    // midi stuff
    float _pb_range = 1;
    float _pb_value = 0;
    float _expression_value = 0;

    int last_released_voice = 0;

    // parameters
    MoniqueSynthParameters _parameters;
    State &_layer_state = _parameters.layer_state;
    float &_master_detune_semitone = _parameters.master_detune;
    GlideModes &_glide_mode = _parameters.glide_mode;

    VoiceMode &_voice_mode = _parameters._voice_mode;
    float _at_value = 0;
    std::array<VoiceData, NUM_VOICES> _voice_data;
    std::array<int, NUM_VOICES> &_num_voices_array;
    WavetableLoader _wave_loader_a = WavetableLoader();
    WavetableLoader _wave_loader_b = WavetableLoader();
    std::array<MoniqueVoice, NUM_VOICES> _voice_array = {
        MoniqueVoice(0, _voice_data.at(0), _parameters, _wave_loader_a, _wave_loader_b, _gp, _drift_gen),
        MoniqueVoice(1, _voice_data.at(1), _parameters, _wave_loader_a, _wave_loader_b, _gp, _drift_gen),
        MoniqueVoice(2, _voice_data.at(2), _parameters, _wave_loader_a, _wave_loader_b, _gp, _drift_gen),
        MoniqueVoice(3, _voice_data.at(3), _parameters, _wave_loader_a, _wave_loader_b, _gp, _drift_gen),
        MoniqueVoice(4, _voice_data.at(4), _parameters, _wave_loader_a, _wave_loader_b, _gp, _drift_gen),
        MoniqueVoice(5, _voice_data.at(5), _parameters, _wave_loader_a, _wave_loader_b, _gp, _drift_gen)};
};

} // namespace Monique
