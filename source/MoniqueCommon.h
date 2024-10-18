#pragma once
#ifndef __x86_64__
#include "arm_neon.h"
#endif
#include "common.h"
#include "monique_synth_parameters.h"
#include "pluginterfaces/vst/ivstevents.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <sys/types.h>

#define EN_DEBUG

namespace Monique {

constexpr int SAMPLE_RATE = 96000;
constexpr int CV_SAMPLE_RATE = 12000;
constexpr int BUFFER_SIZE = 128;
constexpr float BUFFER_RATE = (float)SAMPLE_RATE / (float)BUFFER_SIZE;
constexpr int CV_BUFFER_SIZE = BUFFER_SIZE * CV_SAMPLE_RATE / SAMPLE_RATE;
constexpr int CV_TO_SR_MULT = SAMPLE_RATE / CV_SAMPLE_RATE;
constexpr float MIN_VEL = (1.0 / 127.0) / 2.0;
constexpr int OS_MULT_MODEL = 4;
static std::string cal_path = "/udata/delia/calibration/";
constexpr float EXT_IN_GAIN_MULT = 15.f;
constexpr float VINTAGE_SCALE = 1.f;
constexpr float DRIFT_SCALE = 1.f / 9000.f;
constexpr float DRIFT_OFFSET = 0.0f;
constexpr int PARAPHONIC_MAX_DIST = 10;
constexpr int SHORT_TRIG_TIME_THRESH = 0.1 * BUFFER_RATE;

typedef std::array<float, BUFFER_SIZE> audioBuffer;

struct VoiceOutput {
    std::array<float, BUFFER_SIZE> *output_1;
    std::array<float, BUFFER_SIZE> *output_2;

    VoiceOutput(std::array<float, BUFFER_SIZE> *a1, std::array<float, BUFFER_SIZE> *a2) {
        output_1 = a1;
        output_2 = a2;
    }
};

inline float _int16SampleToFloat(int16_t sample) {
    return static_cast<float>(sample) / static_cast<float>(std::numeric_limits<std::int16_t>::max());
};

inline const float _compareFloats(float f1, float f2) {
    return std::abs(f1 - f2) <= std::numeric_limits<float>::epsilon();
}

struct ProcessAudioData {
    // outputs
    std::array<float, BUFFER_SIZE> *output_l;
    std::array<float, BUFFER_SIZE> *output_r;

    // array of pointers to arrays for the voice data
    std::array<std::array<float, BUFFER_SIZE> *, NUM_VOICES> voice_high_res;
    std::array<std::array<float, BUFFER_SIZE> *, NUM_VOICES> voice_cv;

    // inputs
    std::array<float, BUFFER_SIZE> *input_l;
    std::array<float, BUFFER_SIZE> *input_r;
    std::array<float, BUFFER_SIZE> *input_1;
    std::array<float, BUFFER_SIZE> *input_2;
};

struct ParamChange {
    uint param_id;
    float value;

    ParamChange() :
        param_id(0), value(0.0f) {}

    ParamChange(uint p, float v, bool e = false) :
        param_id(p), value(v) {
        if (e) {
            // Set the MS nibble to indicate this should be exported to
            // any notification listeners
            param_id |= 0x80000000;
        }
    }
};

struct GlobalSynthParameters {
    CvInputMode cv_1_mode = CvInputMode::NEG_10_TO_10;
    CvInputMode cv_2_mode = CvInputMode::NEG_10_TO_10;
    float cv_1_gain_setting = 1.f;
    float cv_2_gain_setting = 1.f;
    float cv_1_offset_setting = 0.f;
    float cv_2_offset_setting = 0.f;
    std::array<float, CV_BUFFER_SIZE> _cv_1;
    std::array<float, CV_BUFFER_SIZE> _cv_2;
};

struct
    MoniqueSynthParameters {
    std::array<float, Monique::NUM_GLOBAL_PARAMS> global_params;
    std::array<float, Monique::NUM_PRESET_COMMON_PARAMS> preset_common_params;
    std::array<float, Monique::NUM_LAYER_PARAMS> layer_params;
    std::array<float, Monique::NUM_COMMON_PARAMS> common_params;
    std::array<float, TOTAL_STATE_PARAMS> state_a = normalised_state_param_init();
    std::array<float, TOTAL_STATE_PARAMS> state_b = normalised_state_param_init();
    std::array<float, TOTAL_STATE_PARAMS> state_a_smooth = normalised_state_param_init();
    std::array<float, TOTAL_STATE_PARAMS> state_b_smooth = normalised_state_param_init();
    std::array<bool, TOTAL_STATE_PARAMS> state_morphed_blocked = {false};
    std::array<float, TOTAL_STATE_PARAMS> state_params = normalised_state_param_init();
    State layer_state = State::A;
    float morph_value = 0;
    float morph_value_smooth = 0;
    float ext_in_gain = 1.f;
    float master_detune = 0;
    float octave_offset = 0;
    LfoOscShape lfo_1_a_shape = LfoOscShape::SINE;
    LfoOscShape lfo_1_b_shape = LfoOscShape::SINE;
    LfoOscShape lfo_2_a_shape = LfoOscShape::SINE;
    LfoOscShape lfo_2_b_shape = LfoOscShape::SINE;
    LfoOscShape lfo_3_a_shape = LfoOscShape::SINE;
    LfoOscShape lfo_3_b_shape = LfoOscShape::SINE;
    GlideModes glide_mode = GlideModes::LINEAR;
    ParaphonicModes para_mode = ParaphonicModes::DISABLE;
    OutputRouting output_mode = OutputRouting::ONE_TWO;
    float glide_rate = 1.0;
    float layer_volume = 1.0;
    float layer_drift = 0.f;
    float layer_vintage = 0.f;
    bool all_notes_off = true;
    bool sustain_ped = false;
    bool lfo_1_global = false;
    bool lfo_2_global = false;
    bool lfo_3_global = false;
    float lfo_1_global_phase = 0;
    float lfo_2_global_phase = 0;
    float lfo_3_global_phase = 0;
    bool lfo_1_sync = false;
    bool lfo_2_sync = false;
    bool lfo_3_sync = false;
    bool overdrive = false;
    bool wt_slow_mode = false;
    bool wt_interpolate = false;
    float pitch_bend = 0;
    float aftertouch = 0;
    float expression = 0;
    float tempo = 60;
    NoiseModes noise_mode = NoiseModes::WHITE;
    VoiceMode _voice_mode = VoiceMode::POLY;
};

struct PolyPressure {
    int16_t channel;
    int16_t note_pitch;
    float pressure;
    int note_id;

    PolyPressure(Steinberg::Vst::PolyPressureEvent e) {
        channel = e.channel;
        pressure = e.pressure;
        note_pitch = e.pitch;
        note_id = e.noteId;
    }
};

struct MidiNote {
    enum type {
        ON,
        OFF
    };

  public:
    int pitch;
    float velocity;
    float note_off_vel;
    int channel;
    type note_type;

    MidiNote() = default;
    ~MidiNote() = default;

    MidiNote(Steinberg::Vst::NoteOnEvent note) {
        // Set the note on attibutes
        pitch = note.pitch;
        velocity = note.velocity;
        channel = note.channel;
        note_off_vel = 0.0;
        note_type = ON;
    }

    MidiNote(Steinberg::Vst::NoteOffEvent note) {
        // Set the note off attibutes
        pitch = note.pitch;
        velocity = note.velocity;
        channel = note.channel;

        // Accoring to MIDI spec, if rel vel is not avaliable, it should be set to 0.5
        if (note.velocity > MIN_VEL) {
            note_off_vel = note.velocity;
        } else if (note.velocity > 1.0) {
            note_off_vel = 1.0;
        } else {
            note_off_vel = 0.5;
        }
        note_type = OFF;
    }
};

struct MpeChannelData {
    uint channel = 0;
    float mpe_x = 0.5;
    float mpe_y = .0;
    float mpe_z = .0;
};

enum class MoniqueFilterModes {
    STEREO_PARALLEL,
    MONOSERIES_LP_HP,
    PARALLEL_LP_HP,
    NUM_MODES
};

enum class FilterSlope {
    MOOG_2_POLE,
    MOOG_4_POLE,
    NUM_SLOPES
};
constexpr int CV_MUX_INC = BUFFER_SIZE / CV_BUFFER_SIZE;

enum VoiceBitMap {
    VOICE_MUTE_L,
    VOICE_MUTE_R,
    VOICE_MUTE_3,
    VOICE_MUTE_4,
    HARD_SYNC,
    DRIVE_EN_N,
    SUB_OSC_EN_N,
    MIX_MUTE_L,
    MIX_MUTE_R,
    FILTER_TYPE
};

struct MoniqueAnalogModelData {
    float res_high_clip = 0.f;
    float res_zero_offset = 1.f;
    float res_gain = -2;
    float res_low_clip = -1;

    float a = 2;
    float c = -1;

    double cut_a = -0.7409;
    double cut_b = 0.19969;
    double cut_c = .0008244;
    double cut_d = 11.25;
    double cut_e = -.4975;
    double cut_f = 3.40697330e-01;
    double cut_g = 8.79614120e-02;
    double cut_h = -1.61918083e-05;
    float cut_low_clip = -0.99;

    float vca_l_offset = 0;
    float vca_r_offset = 0;
    float voice_offset = 0;
    float voice_offset_od = 0;
};

inline double filter_model_function(const MoniqueAnalogModelData &m, double x, const bool two_pole) {
    constexpr float high_clip = .88;
    float two_pole_offset = -0.16;
    x = x + two_pole_offset * !two_pole;

    return std::clamp((float)(m.cut_a + fabs(m.cut_b) * x + m.cut_c * std::exp2(m.cut_d * x + m.cut_e)), m.cut_low_clip, high_clip);
}

enum Cv1Order {
    Cv1MixXor,
    Cv1Drive,
    Cv1FilterCut,
    Cv1FilterQ,
    Cv1AmpL,
    Cv1AmpR,
    Cv1Unused,
    Cv1BitArray
};

struct VoiceData {
    bool overdrive;
    bool hard_sync;
    bool sub_osc;
    bool mute_1;
    bool mute_2;
    bool last_allocated = false;
    bool unison_osc_enable = false;
    bool unison_select_alt_osc = false;
    MoniqueFilterModes filter_mode = MoniqueFilterModes::STEREO_PARALLEL;
    FilterSlope filter_1_slope = FilterSlope::MOOG_4_POLE;
    FilterSlope filter_2_slope = FilterSlope::MOOG_4_POLE;
    NoiseModes noise_mode = NoiseModes::WHITE;
    float aux_level = 1.f;
    float post_hp_gain = 1.f;
    float unison_offset = 0.f;
    float osc_1_blend_tri = 0;
    float osc_2_blend_tri = 0;
    float osc_1_blend_sqr = 0;
    float osc_2_blend_sqr = 0;
    float osc_1_pitch_tmp;
    float osc_1_shape_tmp;
    float osc_2_pitch_tmp;
    float osc_2_shape_tmp;
    float osc_3_pitch_tmp;
    float osc_3_shape_tmp;
    float tri_1_lev_tmp;
    float sqr_1_lev_tmp;
    float tri_2_lev_tmp;
    float sqr_2_lev_tmp;
    float osc_4_tone_tmp;
    float xor_lev_tmp;
    float noise_lev_tmp;
    float osc_3_lev_tmp;
    float filt_1_cut_tmp;
    float filt_1_res_tmp;
    float filt_2_cut_tmp;
    float filt_2_res_tmp;
    float unison_1_tmp;
    float unison_2_tmp;
    float vca_l_tmp;
    float vca_r_tmp;
    float *fm_1_1;
    float *fm_1_2;
    float *fm_2_1;
    float *fm_2_2;
    std::array<float, CV_BUFFER_SIZE> unison_1;
    std::array<float, CV_BUFFER_SIZE> unison_2;
    std::array<float, CV_BUFFER_SIZE> osc_1_pitch;
    std::array<float, CV_BUFFER_SIZE> osc_1_shape;
    std::array<float, CV_BUFFER_SIZE> osc_2_pitch;
    std::array<float, CV_BUFFER_SIZE> osc_2_shape;
    std::array<float, CV_BUFFER_SIZE> osc_3_pitch;
    std::array<float, CV_BUFFER_SIZE> osc_3_shape;
    std::array<float, CV_BUFFER_SIZE> tri_1_lev;
    std::array<float, CV_BUFFER_SIZE> sqr_1_lev;
    std::array<float, CV_BUFFER_SIZE> tri_2_lev;
    std::array<float, CV_BUFFER_SIZE> sqr_2_lev;
    std::array<float, CV_BUFFER_SIZE> xor_lev;
    std::array<float, CV_BUFFER_SIZE> noise_lev;
    std::array<float, CV_BUFFER_SIZE> osc_3_lev;
    std::array<float, CV_BUFFER_SIZE> filt_1_cut;
    std::array<float, CV_BUFFER_SIZE> filt_1_res;
    std::array<float, CV_BUFFER_SIZE> filt_2_cut;
    std::array<float, CV_BUFFER_SIZE> filt_2_res;
    std::array<float, CV_BUFFER_SIZE> vca_l;
    std::array<float, CV_BUFFER_SIZE> vca_r;
    std::array<float, BUFFER_SIZE> wt_buffer;

    void copy_to_buffer(const uint i) {
        unison_1[(i)] = unison_1_tmp;
        unison_2[(i)] = unison_2_tmp;
        osc_1_pitch[(i)] = osc_1_pitch_tmp;
        osc_1_shape[(i)] = osc_1_shape_tmp;
        osc_2_pitch[(i)] = osc_2_pitch_tmp;
        osc_2_shape[(i)] = osc_2_shape_tmp;
        osc_3_pitch[(i)] = osc_3_pitch_tmp;
        osc_3_shape[(i)] = osc_3_shape_tmp;
        tri_1_lev[(i)] = tri_1_lev_tmp;
        sqr_1_lev[(i)] = sqr_1_lev_tmp;
        tri_2_lev[(i)] = tri_2_lev_tmp;
        sqr_2_lev[(i)] = sqr_2_lev_tmp;
        xor_lev[(i)] = xor_lev_tmp;
        noise_lev[(i)] = noise_lev_tmp;
        osc_3_lev[(i)] = osc_3_lev_tmp;
        filt_1_cut[(i)] = filt_1_cut_tmp;
        filt_1_res[(i)] = filt_1_res_tmp;
        filt_2_cut[(i)] = filt_2_cut_tmp;
        filt_2_res[(i)] = filt_2_res_tmp;
        vca_l[(i)] = vca_l_tmp;
        vca_r[(i)] = vca_r_tmp;
    }
};

struct ModMatrixData {
    std::array<float *, FAST_DESTS> _fast_dst;
    std::array<float *, FAST_SRCS> _fast_src;
    std::array<float *, FAST_SRCS * FAST_DESTS> _fast_gains;
    std::array<float *, TOTAL_MOD_DSTS> _matrix_dests;
    std::array<float *, TOTAL_MOD_SRC> _matrix_srcs;
    std::array<float *, TOTAL_MOD_DSTS * TOTAL_MOD_SRC> _gains;
};
} // namespace Monique