#pragma once
#include "GuiMsgQueue.h"
#include "MoniqueCommon.h"
#include "MoniqueLayerManager.h"
#include "common.h"
#include "pluginterfaces/vst/ivstevents.h"
#include "public.sdk/source/vst/vstaudioeffect.h"
#include "public.sdk/source/vst/vstparameters.h"
#include <atomic>
#include <thread>

#define stringPluginName "Monique"
#define stringOriginalFilename stringPluginName ".vst3"
#define stringFileDescription "Monique VST for the Melbourne Instruments Synth"
#define stringCompanyName "Melbourne Instruments"
#define stringLegalCopyright "© 2023 Melbourne Instruments"
#define stringLegalTrademarks ""

namespace Steinberg {
namespace Vst {
namespace MoniqueSynth {
struct lowMidiNote {
    bool on = false;
    uint note = 127;
};

class Processor : public AudioEffect {
  public:
    Processor();
    ~Processor();

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
        return (IAudioProcessor *)new Processor();
    }

    /**
     * @brief GUID for the processor
     *
     */
    static FUID uid;

    void _processGuiMsg() {
        Monique::GuiMsgQueue msg_queue;

        // Open the GUI message queue. we need to wait for the gui to finish creating it
        while (!msg_queue.open() && !_exit_gui_msg_thread) {

            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        if (msg_queue.open()) {
            // Run the thread in a loop processing messages until stopped
            while (!_exit_gui_msg_thread) {
                // Wait for a message to be received
                if (_gui_samples_ready) {
                    auto msg = Monique::GuiMsg();
                    float *current_gui_samples = _current_gui_samples.load();

                    // Get a pointer to the samples just processed
                    if (current_gui_samples == _gui_samples_1)
                        std::memcpy(msg.samples, _gui_samples_2, sizeof(msg.samples));
                    else
                        std::memcpy(msg.samples, _gui_samples_1, sizeof(msg.samples));

                    // Send the message
                    bool res = msg_queue.post_msg(msg);
                    if (res) {
                        _gui_samples_ready = false;
                    }
                }

                // Sleep 1ms and check again
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            // Close the GUI message queue
            msg_queue.close();
        }
    }

    void _updateScopeGain() {

        // dynamically adjust the scope gain based on the number of midi notes coming in
        int size = _held_midi_notes.size();
        if (size > 0) {
            if (size < 2)
                _scope_dynamic_gain = 1;
            else if (size < 4)
                _scope_dynamic_gain = .9;
            else if (size < 8)
                _scope_dynamic_gain = .8;
            else
                _scope_dynamic_gain = .6;
        }
    }

  private:
    Monique::ProcessAudioData _audio_data;
    Monique::LayerManager _layer_manager;

    void debugPrinting(float data1[12][Monique::BUFFER_SIZE]);
    bool _print_data = false;
    float print_data1[12][Monique::BUFFER_SIZE];
    std::atomic<bool> _exit_print_thread;
    std::thread *_printThread;
    int nina_buffer_counter = 0;

    // scope related variables
    float _gui_samples_1[common::GUI_SAMPLES_SIZE];
    float _gui_samples_2[common::GUI_SAMPLES_SIZE];
    std::array<float, common::GUI_NUM_SAMPLES> _scope_ringbuffer_l;
    std::array<float, common::GUI_NUM_SAMPLES> _scope_ringbuffer_r;
    uint _scope_ring_start = 0;
    uint _scope_ring_end = common::GUI_NUM_SAMPLES / 2;
    std::atomic<float *> _current_gui_samples;
    std::atomic<bool> _gui_samples_ready;

    std::thread *_gui_msg_thread;
    std::atomic<bool> _exit_gui_msg_thread;
    float dc_filter_l = 0;
    float dc_filter_r = 0;
    float phase_offset = 0;
    float phase = 0;
    uint _current_samples_written = 0;
    int _gui_buffer_counter = 0;
    bool _gui_buffer_complete;
    bool _zero_crossing_detect = false;
    lowMidiNote _current_low_midi_note;
    std::vector<uint> _held_midi_notes;
    float _scope_dynamic_gain = 1.0;
    float trigger_lp_var = 0;
    float trigger_hp_var = 0;
    float trigger_out = 0;
    float _previous_scope_trigger_sample = 0;
    int _scope_auto_counter = 0;
    static constexpr float SCOPE_HP_FREQ = 500.f;

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
    bool processParameterChanges(IParameterChanges *param_changes, IParameterChanges *out_param_changes);

    void processAudio(ProcessData &data);
};

} // namespace MoniqueSynth
} // namespace Vst
} // namespace Steinberg
