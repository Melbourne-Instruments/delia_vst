/**
 * @file MoniqueProcessor.cpp
 * @brief Monique Processor implementation.
 *
 * @copyright Copyright (c) 2022-2024 Melbourne Instruments, Australia
 */
#include "MoniqueProcessor.h"
#include "MoniqueController.h"
#include "base/source/fstreamer.h"
#include "pluginterfaces/base/ftypes.h"
#include "pluginterfaces/base/futils.h"
#include "pluginterfaces/base/ibstream.h"
#include "pluginterfaces/base/ustring.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"

#include <cassert>
#include <math.h>

//#define PRINTSIGNALS

// --- VST2 wrapper
namespace Steinberg {
namespace Vst {
namespace MoniqueSynth {

FUID Processor::uid(0x66F341CA, 0x882E464E, 0xB274407F, 0xBFC8C024);

Processor::Processor() {
    setControllerClass(Controller::uid);
    // start the GUI message thread & init scope vars/note storage
    _exit_gui_msg_thread = false;
    _gui_buffer_complete = false;
    _current_gui_samples.store(_gui_samples_1);
    _gui_msg_thread = new std::thread(&Processor::_processGuiMsg, this);

    // preallocate max notes * max midi channels
    _held_midi_notes.reserve(127 * 16);
#ifdef PRINTSIGNALS
    _exit_print_thread = false;
    _printThread = new std::thread(&Processor::debugPrinting, this, print_data1);
#endif
}

Processor::~Processor() {
    // Stop the Layer Manager synth processing - makes sure any threads
    // are tidied-up
    _layer_manager.stopSynth();

#ifdef PRINTSIGNALS
    // If the print thread is running
    if (_printThread) {
        // Stop the print thread
        _exit_print_thread = true;
        _printThread->join();
        delete _printThread;
    }
#endif
    // exit the scope thread
    if (_gui_msg_thread) {
        // Stop the GUI message thread
        _exit_gui_msg_thread = true;
        _gui_msg_thread->join();
    }
}

tresult PLUGIN_API Processor::initialize(FUnknown *context) {
    tresult result = AudioEffect::initialize(context);
    if (result == kResultTrue) {

        addAudioOutput(STR16("FX L"), SpeakerArr::kMono);
        addAudioOutput(STR16("FX R"), SpeakerArr::kMono);

        addAudioOutput(STR16("Voice 0 Audio"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 1 Audio"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 2 Audio"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 3 Audio"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 4 Audio"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 5 Audio"), SpeakerArr::kMono);

        addAudioOutput(STR16("Voice 0 CV"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 1 CV"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 2 CV"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 3 CV"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 4 CV"), SpeakerArr::kMono);
        addAudioOutput(STR16("Voice 5 CV"), SpeakerArr::kMono);

        // audio input for diagnostic audio and tuning measurements
        addAudioInput(STR16("Mix Left Input"), SpeakerArr::kMono);
        addAudioInput(STR16("Mix Right Input"), SpeakerArr::kMono);
        addAudioInput(STR16("External In 1"), SpeakerArr::kMono);
        addAudioInput(STR16("External In 2"), SpeakerArr::kMono);

        // MIDI event input bus, one channel
        addEventInput(STR16("Event Input"), 16);
#ifdef SIGDUMP
        buf.reserve(96000);
        for (int i = 0; i < 96000; i++) {
            buf.push_back(-2.0);
        }
#endif
    }

    return result;
}

tresult PLUGIN_API Processor::setState(IBStream *fileStream) {

    return kResultTrue;
}

tresult PLUGIN_API Processor::getState(IBStream *fileStream) {

    return kResultTrue;
}

tresult PLUGIN_API Processor::setBusArrangements(SpeakerArrangement *inputs,
    int32 numIns,
    SpeakerArrangement *outputs,
    int32 numOuts) {

    if (numIns == 4 && numOuts == 14 && outputs[0] == SpeakerArr::kMono) {
        return AudioEffect::setBusArrangements(inputs, numIns, outputs, numOuts);
    }
    return kResultFalse;
}

tresult PLUGIN_API Processor::canProcessSampleSize(int32 symbolicSampleSize) {

    // --- currently 32 bit only
    if (symbolicSampleSize == kSample32) {
        return kResultTrue;
    }
    return kResultFalse;
}

tresult PLUGIN_API Processor::setActive(TBool state) {
    if (state) {
    } else {
    }

    // call base class method
    return AudioEffect::setActive(state);
}

tresult PLUGIN_API Processor::process(ProcessData &data) {
    // Process all events
    processEvents(data.inputEvents);
    processParameterChanges(data.inputParameterChanges, data.outputParameterChanges);
    processAudio(data);

    return kResultTrue;
}

bool Processor::processParameterChanges(IParameterChanges *param_changes, IParameterChanges *out_param_changes) {
    // Is the param changes pointer specified?
    if (param_changes) {
        // Get the number of changes, and check if there are any to process
        int32 count = param_changes->getParameterCount();
        if (count > 0) {
            // printf("\n count %d \n", count);
            //  Process each param change
            for (int32 i = 0; i < count; i++) {
                // Get the queue of changes for this parameter, if no queue is
                // returned then skip this parameter change
                IParamValueQueue *queue = param_changes->getParameterData(i);
                if (!queue)
                    continue;

                // Get the param ID and if valid process
                ParamID paramId = queue->getParameterId();

                // Get the last (latest) point value in the queue.
                int32 sampleOffset;
                ParamValue value;
                queue->getPoint((queue->getPointCount() - 1), sampleOffset, value);
                _layer_manager.updateParameter(paramId, static_cast<float>(value));
                // if (paramId == 2182 || paramId == 134)
            }
            auto &changes = _layer_manager.paramRefresh();

            return true;
        }
    }
    auto &changes = _layer_manager.getChanges();
    int num_of_changes = (int)changes.size();

    // if there are alot of changes, we send them 256 at a time to avoid any buffer skips
    if (num_of_changes > Monique::MAX_PARAM_CHANGES) {
        num_of_changes = Monique::MAX_PARAM_CHANGES;
    }
    for (int i = 0; i < num_of_changes; i++) {
        const Monique::ParamChange &pc = changes.back();
        int32 index = 0;
        auto paramQueue = out_param_changes->addParameterData(pc.param_id, index);
        if (paramQueue != nullptr) {
            int32 index2 = 0;
            paramQueue->addPoint(0, pc.value, index2);
        }
        changes.pop_back();
    }
    return true;
}

static float phase = 0;

void Processor::processAudio(ProcessData &data) {
    float sum = 0;
    // assign data buffers
    _audio_data.output_l = (std::array<float, Monique::BUFFER_SIZE> *)data.outputs->channelBuffers32[0];
    _audio_data.output_r = (std::array<float, Monique::BUFFER_SIZE> *)data.outputs->channelBuffers32[1];

    // assign voice output buffers for each voice
    for (int i = 0; i < Monique::NUM_VOICES; i++) {
        _audio_data.voice_high_res.at(i) = (std::array<float, Monique::BUFFER_SIZE> *)data.outputs->channelBuffers32[2 + i];
        _audio_data.voice_cv.at(i) = (std::array<float, Monique::BUFFER_SIZE> *)data.outputs->channelBuffers32[2 + Monique::NUM_VOICES + i];
    }
    // assign input buffers
    _audio_data.input_l = (std::array<float, Monique::BUFFER_SIZE> *)data.inputs->channelBuffers32[0];
    _audio_data.input_r = (std::array<float, Monique::BUFFER_SIZE> *)data.inputs->channelBuffers32[1];
    _audio_data.input_1 = (std::array<float, Monique::BUFFER_SIZE> *)data.inputs->channelBuffers32[2];
    _audio_data.input_2 = (std::array<float, Monique::BUFFER_SIZE> *)data.inputs->channelBuffers32[3];

    // Process audio via the Layer Manager
    _layer_manager.setTempo((float)data.processContext->tempo);
    _layer_manager.runSynth(_audio_data);

#ifdef PRINTSIGNALS
    nina_buffer_counter++;
    // save the data for display
    if (nina_buffer_counter % 150 == 0) {
        for (int voice_num = 0; voice_num < Monique::NUM_VOICES; voice_num++) {
            for (int i = 0; i < Monique::BUFFER_SIZE; i++) {
                float sample = _audio_data.voice_cv.at(voice_num)->at(i);
                print_data1[voice_num][i] = sample;
            }
        }
        _print_data = true;
        nina_buffer_counter = 0;
    }
#endif

    // Get the GUI samples for the scope
    float *left = &(_audio_data.input_l->at(0));
    float *right = &(_audio_data.input_r->at(0));
    float *current_gui_samples = _current_gui_samples.load();

    // scope trigger filter setup
    static constexpr float x = 2 * M_PI * (SCOPE_HP_FREQ) / (float)Monique::SAMPLE_RATE;
    static constexpr float p_lp = (2 - std::cos(x)) - std::sqrt((2 - std::cos(x)) * (2 - std::cos(x)) - 1);

    if (!_gui_buffer_complete) {
        while (phase < (float)Monique::BUFFER_SIZE) {

            // calc time step for current window based on midi note length
            const float window = 1.0 / (common::GUI_NUM_SAMPLES * Monique::midiNoteToFreq(_current_low_midi_note.note));
            float phase_inc = (float)Monique::SAMPLE_RATE / (0.5 / window);

            // current sample to add to buffer
            const int samp = std::floor(phase);
            phase += phase_inc;

            // filter the LR samples
            float current_sample_l = 8 * (left[samp]);
            dc_filter_l = (1 - p_lp) * current_sample_l + p_lp * dc_filter_l;
            current_sample_l = current_sample_l - dc_filter_l;
            float current_sample_r = 8 * (right[samp]);
            dc_filter_r = (1 - p_lp) * current_sample_r + p_lp * dc_filter_r;
            current_sample_r = current_sample_r - dc_filter_r;
            float current_sample_sum = current_sample_l + current_sample_r;
            trigger_out = current_sample_sum;

            // if we have written at least 1/2 a buffer worth of data, then we can start looking for zero crossings
            if (!_zero_crossing_detect && (_current_samples_written > (uint)(common::GUI_NUM_SAMPLES / 2))) {

                const bool sign_1 = current_sample_sum < 0.f;
                const bool n_sign_2 = _previous_scope_trigger_sample > 0.f;
                const bool zero_x = sign_1 && n_sign_2;

                // if we find a zero crossing or the level is very low, then trigger a capture
                bool auto_trig = (_scope_auto_counter > common::GUI_NUM_SAMPLES / 2);
                if (zero_x || auto_trig) {
                    _zero_crossing_detect = true;
                    _scope_ring_end = _scope_ring_start + (common::GUI_NUM_SAMPLES / 2);
                    if (_scope_ring_end >= common::GUI_NUM_SAMPLES) {
                        _scope_ring_end -= common::GUI_NUM_SAMPLES;
                    }
                    // reset the auto counter
                    _scope_auto_counter = 0;
                }
                _scope_auto_counter++;
            }
            _previous_scope_trigger_sample = trigger_out;

            // write the samples into the ring buffer
            const int offset = _gui_buffer_counter++;
            _scope_ringbuffer_l[(_scope_ring_start)] = current_sample_l * _scope_dynamic_gain;
            _scope_ringbuffer_r[(_scope_ring_start)] = current_sample_r * _scope_dynamic_gain;
            _current_samples_written++;
            _scope_ring_start++;
            if (_scope_ring_start == common::GUI_NUM_SAMPLES) {
                _scope_ring_start = 0;
            }

            // if we have filled the ring buffer, then trigger the gui message transfer
            if (_zero_crossing_detect && (_scope_ring_start == _scope_ring_end)) {
                _zero_crossing_detect = false;
                _gui_buffer_complete = true;
                break;
            }
        }

        // reset the sample phase
        if (!_gui_buffer_complete) {
            phase -= 128.0;
        }
    }

    // All GUI buffers processed and ready to send samples?
    if (_gui_buffer_complete) {
        // printf("\nbuff send");
        //  If we haven't processed the previous GUI buffer yet, skip this one
        if (!_gui_samples_ready) {
            // Load the ring buffer samples into the GUI buffer in the ringbuffer order
            float *current_samples = _current_gui_samples;
            uint ring_counter = _scope_ring_end;
            for (uint i = 0; i < common::GUI_NUM_SAMPLES; i++) {
                current_samples[i * 2] = _scope_ringbuffer_l[(ring_counter)];
                current_samples[i * 2 + 1] = _scope_ringbuffer_r[(ring_counter++)];
                if (ring_counter == common::GUI_NUM_SAMPLES) {
                    ring_counter -= common::GUI_NUM_SAMPLES;
                }
            }

            // Swap the pointer to the other buffer to process the next samples
            if (_current_gui_samples == _gui_samples_1)
                _current_gui_samples = _gui_samples_2;
            else
                _current_gui_samples = _gui_samples_1;
            _gui_samples_ready = true;
            _zero_crossing_detect = false;
            _gui_buffer_complete = false;
            _current_samples_written = 0;
        } else {
        }
    }
}

void Processor::debugPrinting(float data1[Monique::NUM_VOICES][Monique::BUFFER_SIZE]) {
    printf("\nstart printing thread\n");
    while (!_exit_print_thread) {
        if (_print_data) {
            _print_data = false;
            static const int printvoices = Monique::NUM_VOICES;

            printf("\n\n\t\t");
            for (int i = 0; i < printvoices; i++) {
                printf("v:%d\t", i);
            }

            printf("\n Amp L\t\t");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data1[i][Monique::Cv1Order::Cv1AmpL];
                printf("%1.2f\t", tmp);
            }
            printf("\n Amp R\t\t");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data1[i][Monique::Cv1Order::Cv1AmpR];
                printf("%1.2f\t", tmp);
            }
            printf("\n F Cut\t\t");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data1[i][Monique::Cv1Order::Cv1FilterCut];
                printf("%1.2f\t", tmp);
            }
            printf("\n Res\t\t");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data1[i][Monique::Cv1Order::Cv1FilterQ];
                printf("%1.2f\t", tmp);
            }
            printf("\n bit array\t");
            for (int i = 0; i < printvoices; i++) {
                int tmp = data1[i][Monique::Cv1Order::Cv1BitArray + 10 * 8] * (1 / 1.19209303761637659268e-7f);
                printf("%x\t", tmp);
            }

            printf("\n");
        }
        struct timespec remaining, request = {1, 500 * 1000};
        nanosleep(&request, &remaining);
    }
    printf("\n exit printing thread\n");
}

void Processor::processEvents(IEventList *events) {
    // If events are specified
    if (events) {
        // Process the events
        auto event_count = events->getEventCount();
        for (int i = 0; i < event_count; i++) {
            // If the maximum number of notes have not been processed
            Event event;
            // Get the event
            auto res = events->getEvent(i, event);
            if (res == kResultOk) {
                // Parse the event type
                switch (event.type) {
                case Event::kNoteOnEvent: {
                    // printf("note on %d\n ", event.noteOn.pitch);
                    // Handle the MIDI note on event
                    const auto note = Monique::MidiNote(event.noteOn);
                    _layer_manager.handleNote(note);
                    // update scope note vars
                    _held_midi_notes.push_back(note.pitch);
                    if (note.pitch < _current_low_midi_note.note)
                        _current_low_midi_note.note = note.pitch;
                    _updateScopeGain();
                    if (!_current_low_midi_note.on) {
                        _current_low_midi_note.on = true;
                        _current_low_midi_note.note = note.pitch;
                        _updateScopeGain();
                    }
                    break;
                }

                case Event::kNoteOffEvent: {
                    // Handle the MIDI note off event
                    // printf("note off %d\n ", event.noteOn.pitch);
                    const auto note = Monique::MidiNote(event.noteOff);
                    uint i = 0;
                    while (i < _held_midi_notes.size()) {
                        if (_held_midi_notes.at(i) == note.pitch) {
                            _held_midi_notes.erase(_held_midi_notes.begin() + i);
                        } else {
                            i++;
                        }
                    }
                    _layer_manager.handleNote(Monique::MidiNote(event.noteOff));
                    if (note.pitch == _current_low_midi_note.note)
                        _current_low_midi_note.on = false;
                    break;
                }

                case Event::kPolyPressureEvent: {
                    _layer_manager.polyPressure(event.polyPressure);
                }

                default:
                    // Ignore all other events
                    break;
                }
            }
        }
    }
}

} // namespace MoniqueSynth
} // namespace Vst
} // namespace Steinberg
