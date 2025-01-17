#pragma once

#include "MoniqueAnalogModel.h"
#include "MoniqueCommon.h"
#include "NoiseOscillator.h"
#include "common.h"
#include "pluginterfaces/vst/ivstevents.h"
#include "public.sdk/source/vst/vstaudioeffect.h"
#include "public.sdk/source/vst/vstparameters.h"
#include <algorithm>

#include <atomic>
#include <thread>

namespace Steinberg {
namespace Vst {
namespace MoniqueFactory {

// type imports
typedef Monique::MidiNote MidiNote;
auto constexpr SAMPLE_RATE = Monique::SAMPLE_RATE;
auto constexpr BUFFER_SIZE = Monique::BUFFER_SIZE;
auto constexpr NUM_VOICES = Monique::NUM_VOICES;

//#define DEBUGPRINTING
struct voice_settings {
    float noise = 0;
    float sine = 0;
};

constexpr int num_voice_params = 10;
constexpr int voice_param_offset = 3;

// Nina Processor class
class FactoryTestProcessor : public AudioEffect {
  public:
    FactoryTestProcessor();
    ~FactoryTestProcessor();

    /**
     * @brief
     *
     * @param context
     * @return tresult
     */
    tresult PLUGIN_API initialize(FUnknown *context);

    /**
     * @brief Set the Bus Arrangements object
     *
     * @param inputs speaker arrangement
     * @param numIns number of inputs
     * @param outputs speaker arrangement
     * @param numOuts
     * @return tresult
     */
    tresult PLUGIN_API setBusArrangements(SpeakerArrangement *inputs,
        int32 numIns,
        SpeakerArrangement *outputs,
        int32 numOuts);

    /**
     * @brief
     *
     * @param symbolicSampleSize
     * @return tresult
     */
    tresult PLUGIN_API canProcessSampleSize(int32 symbolicSampleSize);

    /**
     * @brief Set the Active object
     *
     * @param state
     * @return tresult
     */
    tresult PLUGIN_API setActive(TBool state);

    /**
     * @brief Set the State object
     *
     * @param fileStream
     * @return tresult
     */
    tresult PLUGIN_API setState(IBStream *fileStream);
    /**
     * @brief Get the State object
     *
     * @param fileStream
     * @return tresult
     */
    tresult PLUGIN_API getState(IBStream *fileStream);

    /**
     * @brief
     *
     * @param data
     * @return tresult
     */
    tresult PLUGIN_API process(ProcessData &data);

    /**
     * @brief Create a Instance object
     *
     * @return FUnknown*
     */
    static FUnknown *createInstance(void *) {
        return (IAudioProcessor *)new FactoryTestProcessor();
    }

    void debugPrinting(float data1[12][BUFFER_SIZE], float data2[12][BUFFER_SIZE], float data3[BUFFER_SIZE]);

    /**
     * @brief GUID for the processor
     *
     */
    static FUID uid;

  private:
    /**
     * @brief Voice Manager
     */

    /**
     * @brief Modulation Matrix Manager
     */

    std::array<Monique::VoiceData, NUM_VOICES> _voice_data_array;
    std::array<float, num_voice_params * NUM_VOICES + 4> _params;

    std::array<Monique::MoniqueAnalogModel, NUM_VOICES> _voices = {
        Monique::MoniqueAnalogModel(0, _voice_data_array.at(0)),
        Monique::MoniqueAnalogModel(1, _voice_data_array.at(1)),
        Monique::MoniqueAnalogModel(2, _voice_data_array.at(2)),
        Monique::MoniqueAnalogModel(3, _voice_data_array.at(3)),
        Monique::MoniqueAnalogModel(4, _voice_data_array.at(4)),
        Monique::MoniqueAnalogModel(5, _voice_data_array.at(5))};
    std::array<voice_settings, 12> _voice_settings;

    int _unison_voices = 1;
    float _unison_spread = 0;
    float _fx_sine_phase = 0;
    float _unison_pan = 0;
    bool _print_data = false;
    float _fx_sine_level = 0;
    float print_data1[12][BUFFER_SIZE], print_data2[12][BUFFER_SIZE], print_data3[BUFFER_SIZE];
    Monique::NoiseOscillator _noise_osc;

// use this to dump 1 second of the wavetable signal to a file, could be used to output any signal. very useful for realtime debugging
//#define SIGDUMP
#ifdef SIGDUMP
    int counter3 = 0;
    std::vector<float> buf;
#endif
    std::thread *_printThread;
    std::atomic<bool> _exit_print_thread;

    /**
     * @brief Process VST events (MIDI)
     */
    void processEvents(IEventList *events);

    /**
     * @brief
     *
     * @param param_changes
     * @return true
     * @return false
     */
    bool processParameterChanges(IParameterChanges *param_changes);
    /**
     * @brief Handle MIDI note ON event
     */
    void handleMidiNoteOnEvent(const MidiNote &midi_note);

    /**
     * @brief Handle MIDI note OFF event
     */
    void handleMidiNoteOffEvent(const MidiNote &midi_note);
    void processAudio(ProcessData &data);
    int counter_p = 0;
    void printDebugOutput();

    /**
     * @brief helper function to update all the parameter values if there has been
     * any changes
     */
    void _updateParams(const std::vector<ParamID> &changed_param_ids);

    /**
     * @brief Process GUI message thread function
     */
    void _processGuiMsg();
};

} // namespace MoniqueFactory
} // namespace Vst
} // namespace Steinberg