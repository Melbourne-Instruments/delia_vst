#include "FactoryTestProcessor.h"
#include "FactoryTestController.h"
#include "base/source/fstreamer.h"
#include "pluginterfaces/base/ftypes.h"
#include "pluginterfaces/base/futils.h"
#include "pluginterfaces/base/ibstream.h"
#include "pluginterfaces/base/ustring.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include <chrono>
#include <math.h>

namespace Steinberg {
namespace Vst {
namespace MoniqueFactory {

// Constants
constexpr int MAX_MIDI_NOTES = 12;

#define PRINTSIGNALS

FUID FactoryTestProcessor::uid(0xC1A4E445, 0x720643C4, 0x9152F102, 0xFCE80559);

FactoryTestProcessor::FactoryTestProcessor() {
    setControllerClass(FactoryTestController::uid);
#ifdef PRINTSIGNALS
    _exit_print_thread = false;
    _printThread = new std::thread(&FactoryTestProcessor::debugPrinting, this, print_data1, print_data2, print_data3);
    for (int v = 0; v < NUM_VOICES; v++) {
        _params[voice_param_offset + num_voice_params * v + 1] = 0.5;
        _params[voice_param_offset + num_voice_params * v + 2] = 0.5;
    }
#endif
}

FactoryTestProcessor::~FactoryTestProcessor() {
#ifdef PRINTSIGNALS
    // If the print thread is running
    if (_printThread) {
        // Stop the print thread
        _exit_print_thread = true;
        _printThread->join();
        delete _printThread;
    }
#endif
}

tresult PLUGIN_API FactoryTestProcessor::initialize(FUnknown *context) {
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

tresult PLUGIN_API FactoryTestProcessor::setState(IBStream *fileStream) {

    return kResultTrue;
}

tresult PLUGIN_API FactoryTestProcessor::getState(IBStream *fileStream) {

    return kResultTrue;
}

tresult PLUGIN_API FactoryTestProcessor::setBusArrangements(SpeakerArrangement *inputs,
    int32 numIns,
    SpeakerArrangement *outputs,
    int32 numOuts) {

    if (numIns == 0 && numOuts == 14 && outputs[0] == SpeakerArr::kMono) {
        return AudioEffect::setBusArrangements(inputs, numIns, outputs, numOuts);
    }
    return kResultFalse;
}

tresult PLUGIN_API FactoryTestProcessor::canProcessSampleSize(int32 symbolicSampleSize) {

    // --- currently 32 bit only
    if (symbolicSampleSize == kSample32) {
        return kResultTrue;
    }
    return kResultFalse;
}

tresult PLUGIN_API FactoryTestProcessor::setActive(TBool state) {
    if (state) {
    } else {
    }

    // call base class method
    return AudioEffect::setActive(state);
}

tresult PLUGIN_API FactoryTestProcessor::process(ProcessData &data) {
    // Process all events
    processEvents(data.inputEvents);
    processParameterChanges(data.inputParameterChanges);
    processAudio(data);

    return kResultTrue;
}

int buf_counter = 0;

bool FactoryTestProcessor::processParameterChanges(IParameterChanges *param_changes) {
    // Is the param changes pointer specified?
    if (param_changes) {
        // Get the number of changes, and check if there are any to process
        int32 count = param_changes->getParameterCount();
        if (count > 0) {
            std::vector<ParamID> changed_param_ids(count);

            // Process each param change
            for (int32 i = 0; i < count; i++) {
                // Get the queue of changes for this parameter, if no queue is
                // returned then skip this parameter change
                IParamValueQueue *queue = param_changes->getParameterData(i);
                if (!queue)
                    continue;

                // Get the param ID and if valid process
                ParamID paramId = queue->getParameterId();

                // Get the last (latest) point value in the queue, any previous values
                // are not processed
                int32 sampleOffset;
                ParamValue value;
                queue->getPoint((queue->getPointCount() - 1), sampleOffset, value);

                // Set the param value and add to the vector of param IDs changed
                _params.at(paramId) = value;
                changed_param_ids.push_back(paramId);
            }

            // Update the parameters
            _updateParams(changed_param_ids);
            return true;
        }
    }
    return false;
};

std::array<float, NUM_VOICES> sin_phases = {0, 0, 0, 0, 0, 0};

void FactoryTestProcessor::processAudio(ProcessData &data) {
    float *high_res_audio_outputs[NUM_VOICES];
    float *cv_outputs[NUM_VOICES];
    float *fx_left = data.outputs[0].channelBuffers32[0];
    float *fx_right = data.outputs[0].channelBuffers32[1];
    for (int i = 0; i < NUM_VOICES; i++) {
        high_res_audio_outputs[i] = data.outputs[0].channelBuffers32[i + 2];
        cv_outputs[i] = data.outputs[0].channelBuffers32[i + 2 + NUM_VOICES];
    }
    int voice_num = 0;
    _noise_osc.setVolume(1);
    buf_counter++;
    if ((buf_counter % 1000) == 0) {
        _print_data = true;
    }
    auto input_itr = _voice_data_array.begin();
    int i = 0;
    int voicen = 0;
    for (int i = 0; i < BUFFER_SIZE; i++) {

        _fx_sine_phase += 500.f / 96000.f;
        _fx_sine_phase = _fx_sine_phase - (1.0f) * (_fx_sine_phase > 1.0f);
        fx_left[i] = std::sin(_fx_sine_phase * 2.0f * M_PI) * _fx_sine_level;
        fx_right[i] = std::sin(_fx_sine_phase * 2.0f * M_PI) * _fx_sine_level;
    }

    for (auto &voice : _voices) {

        // dodgy test of this method sending the actual pointers to the voice output. data is written directly into the output buffer like this
        voice.run(*(reinterpret_cast<std::array<float, BUFFER_SIZE> *>(cv_outputs[voice_num])));

        for (int i = 0; i < BUFFER_SIZE; i++) {
            const float _noise_osc_src = _noise_osc.getSample();
            sin_phases[voice_num] += 500.f / 96000.f;
            if (sin_phases[voice_num] > 1) {
                sin_phases[voice_num] -= 1;
            }

            high_res_audio_outputs[voice_num][i] = (sin_phases[voice_num] - 0.5) * _voice_settings[voice_num].noise * 2 + std::sin(2 * M_PI * sin_phases[voice_num]) * _voice_settings[voice_num].sine;

            print_data1[voice_num][i] = cv_outputs[voice_num][i];
            print_data2[voice_num][i] = cv_outputs[voice_num][i];
            if (voice_num == 0) {
                float tmp = std::sin(2 * M_PI * sin_phases[voice_num]);
                // printf("\n %f",tmp);
            }
        }
        voice_num++;
    }
};

void FactoryTestProcessor::debugPrinting(float data1[NUM_VOICES][BUFFER_SIZE], float data2[NUM_VOICES][BUFFER_SIZE], float data3[BUFFER_SIZE]) {
    printf("\nstart printing thread\n");
    while (!_exit_print_thread) {
        if (_print_data) {
            _print_data = false;
            static const int printvoices = NUM_VOICES;

            printf("\n\n");
            for (int i = 0; i < printvoices; i++) {
                printf("v:%d\t", i);
            }

            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data2[i][Monique::Cv1Order::Cv1AmpL];
                printf("%1.2f\t", tmp);
            }
            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data2[i][Monique::Cv1Order::Cv1AmpR];
                printf("%1.2f\t", tmp);
            }
            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data2[i][Monique::Cv1Order::Cv1FilterCut];
                printf("%1.2f\t", tmp);
            }
            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = data2[i][Monique::Cv1Order::Cv1FilterQ];
                printf("%1.2f\t", tmp);
            }
            printf("\n bit array\n");
            for (int i = 0; i < printvoices; i++) {
                int tmp = data2[i][Monique::Cv1Order::Cv1BitArray + 10 * 8] * (1 / 1.19209303761637659268e-7f);
                printf("%x\t", tmp);
            }
            const float count_scale = (2 << 22) / (73.75e6 / 2.0); //
            printf("\n freq\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = 1 / (data3[i * 4] * count_scale);
                printf("%1.0f\t", tmp);
            }
            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = 1 / (count_scale * data3[i * 4 + 1]);
                printf("%1.0f\t", tmp);
            }
            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = 1 / (count_scale * data3[i * 4 + 2]);
                printf("%1.0f\t", tmp);
            }
            printf("\n");
            for (int i = 0; i < printvoices; i++) {
                float tmp = 1 / (count_scale * data3[i * 4 + 3]);
                printf("%1.0f\t", tmp);
            }

            printf("\n");
        }
        struct timespec remaining, request = {1, 500 * 1000 * 1000};
        nanosleep(&request, &remaining);
    }
    printf("\n exit printing thread\n");
}

void FactoryTestProcessor::processEvents(IEventList *events) {
    int midi_note_count = 0;

    // If events are specified
    if (events) {
        // Process the events
        auto event_count = events->getEventCount();
        for (int i = 0; i < event_count; i++) {
            // If the maximum number of notes have not been processed
            if (midi_note_count < MAX_MIDI_NOTES) {
                Event event;

                // Get the event
                auto res = events->getEvent(i, event);
                if (res == kResultOk) {
                    // Parse the event type
                    switch (event.type) {
                    case Event::kNoteOnEvent: {
                        // Handle the MIDI note on event
                        handleMidiNoteOnEvent(MidiNote(event.noteOn));
                        midi_note_count++;
                        break;
                    }

                    case Event::kNoteOffEvent: {
                        // Handle the MIDI note off event
                        handleMidiNoteOffEvent(MidiNote(event.noteOff));
                        midi_note_count++;
                        break;
                    }

                    default: {
                        // Ignore all other events
                    }
                    }
                }
            }
        }
    }
}

void FactoryTestProcessor::handleMidiNoteOnEvent(const MidiNote &midi_note) {
    // Allocate voices for this note on
    if (midi_note.velocity > 0) {
    } else {
    }
}

void FactoryTestProcessor::handleMidiNoteOffEvent(const MidiNote &midi_note) {
    // Free the voices associated with this note off
}

auto constexpr noteGain = Monique::NOTE_GAIN;
auto constexpr CV_BUFFER_SIZE = Monique::CV_BUFFER_SIZE;

void FactoryTestProcessor::_updateParams(const std::vector<ParamID> &changed_param_ids) {

    for (int i = 0; i < NUM_VOICES; i++) {
        _voice_settings.at(i).noise = _params[voice_param_offset + num_voice_params * i + 0];

        std::fill_n(_voice_data_array.at(i).filt_2_res.begin(), CV_BUFFER_SIZE, _params[voice_param_offset + num_voice_params * i + 1]);
        std::fill_n(_voice_data_array.at(i).filt_2_cut.begin(), CV_BUFFER_SIZE, _params[voice_param_offset + num_voice_params * i + 2]);
        std::fill_n(_voice_data_array.at(i).vca_l.begin(), CV_BUFFER_SIZE, _params[voice_param_offset + num_voice_params * i + 4]);
        std::fill_n(_voice_data_array.at(i).vca_r.begin(), CV_BUFFER_SIZE, _params[voice_param_offset + num_voice_params * i + 5]);

        _voice_data_array.at(i).overdrive = 0.5 < _params[voice_param_offset + num_voice_params * i + 3];
        _voice_data_array.at(i).mute_1 = 0.5 < _params[voice_param_offset + num_voice_params * i + 6];
        _voice_data_array.at(i).mute_2 = 0.5 < _params[voice_param_offset + num_voice_params * i + 7];

        bool f_12db = 0.5 < _params[voice_param_offset + num_voice_params * i + 9];
        _voice_data_array.at(i).filter_2_slope = f_12db ? Monique::FilterSlope::MOOG_2_POLE : Monique::FilterSlope::MOOG_4_POLE;

        _voices.at(i).fxMuteLoopR(0.5 < _params[1]);
        _voices.at(i).fxMuteLoopL(0.5 < _params[0]);

        _voice_settings.at(i).sine = _params[voice_param_offset + num_voice_params * i + 8];
        _voices.at(i).setMute(_voice_data_array.at(i).mute_1);
        _fx_sine_level = _params[2];
    }
}

void FactoryTestProcessor::_processGuiMsg() {
}

} // namespace MoniqueFactory
} // namespace Vst
} // namespace Steinberg
