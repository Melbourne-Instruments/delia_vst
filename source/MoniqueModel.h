#pragma once

#include "MoniqueAnalogModel.h"
#include "MoniqueCommon.h"
#include "MoniqueFilters.h"
#include "MoniqueWavetable.h"
#include "NoiseOscillator.h"
#include <array>
#include <iostream>
#include <sstream>

namespace Monique {

constexpr int OS_MULT = 2;
constexpr int OS_RATE = SAMPLE_RATE * OS_MULT;

struct ToneFilter {
    float _filter_1_state = 0;
    float _filter_2_state = 0;
    float _out_gain = 1;
    float lgain = 0;
    float hgain = 0;

    void update(float tone) {
        constexpr float amp_db = 6.f;

        constexpr float amp = amp_db / std::log(2.f);
        constexpr float gfactor = 5;
        const float gain = (tone - 1.f) * amp_db;
        float g1 = -gain;
        float g2 = gfactor * gain;

        lgain = fastexp(g1 / amp) - 1;
        hgain = fastexp(g2 / amp) - 1;
        constexpr float comp = 2.2;
        _out_gain = tone > 0.5 ? (1 + comp - (tone - 0.5) * comp * 2) : 1 + comp;
    }

    float run(const float signal) {
        float out;
        // stage 1 filter
        {
            constexpr float sr3 = 3 * OS_RATE;
            constexpr float f0 = 250;
            constexpr float omega = 2 * M_PI * f0;
            constexpr float n = 1 / (sr3 + omega);
            constexpr float a0 = 2 * omega * n;
            constexpr float b1 = (sr3 - omega) * n;
            float &lp_out = _filter_1_state;
            lp_out = a0 * signal + b1 * lp_out;
            out = signal + lgain * lp_out + hgain * (signal - lp_out);
        }

        // stage 2 filter
        {
            constexpr float sr3 = 3 * OS_RATE;
            constexpr float f0 = 5000;
            constexpr float omega = 2 * M_PI * f0;
            constexpr float n = 1 / (sr3 + omega);
            constexpr float a0 = 2 * omega * n;
            constexpr float b1 = (sr3 - omega) * n;
            float &lp_out = _filter_2_state;
            lp_out = a0 * out + b1 * lp_out;
            out = out + lgain * lp_out + hgain * (out - lp_out);
        }

        return fast_tanh(out / 4) * _out_gain * 4;
    }
};

struct filter6db {
    typedef double REAL;
    uint NBQ2 = 1;
    REAL biquada[1] = {-0.8360692637090134};
    REAL biquadb[1] = {1};
    REAL gain = 12.200274611405899;
    REAL xyv[6] = {0, 0, 0, 0, 0, 0};

    inline REAL applyfilter6db(REAL v) {
        int i, b, xp = 0, yp = 2, bqp = 0;
        REAL out = v / gain;
        for (i = 5; i > 0; i--) {
            xyv[i] = xyv[i - 1];
        }
        for (b = 0; b < NBQ2; b++) {
            int len = (b == NBQ2 - 1) ? 1 : 2;
            xyv[xp] = out;
            for (i = 0; i < len; i++) {
                out += xyv[xp + len - i] * biquadb[bqp + i] - xyv[yp + len - i] * biquada[bqp + i];
            }
            bqp += len;
            xyv[yp] = out;
            xp = yp;
            yp += len + 1;
        }
        return out;
    }
};

constexpr int UNISON = 2;
constexpr int UNISON_SHIFT = 1;
constexpr int NUM_OSC = 2;
constexpr int OSC_1 = 0;
constexpr int OSC_2 = UNISON - 1;

struct UnisonMoniqueOsc {
    static constexpr float reset = 1.0;
    bool sub = false;
    bool sync = false;
    std::array<float, NUM_OSC> gain_up;
    std::array<float, NUM_OSC> gain_down;
    std::array<float, NUM_OSC> shape;
    std::array<float, UNISON> xor_out;
    std::array<float, UNISON * NUM_OSC> osc_tri;
    std::array<float, UNISON * NUM_OSC> osc_sqr;
    std::array<float, UNISON * NUM_OSC> phase;
    std::array<float, UNISON * NUM_OSC> phase_inc;

    UnisonMoniqueOsc() {
        xor_out.fill(0.f);
        gain_up.fill(0.f);
        gain_down.fill(0.f);
        shape.fill(0.f);
        osc_tri.fill(0.f);
        osc_sqr.fill(0.f);
        phase.fill(0.f);
        phase_inc.fill(0.f);
    };

    inline void recalculate(const std::array<float, NUM_OSC * UNISON> freq_arr, const std::array<float, NUM_OSC> &shape_arr) {
        //  constants for calculations
        constexpr float n_gain = NOTE_GAIN;
        constexpr float fmax = 25000.f;
        static constexpr float max_osc_freq_l2 = std::log2(30000);
        float freq_store[NUM_OSC];
        for (int i = 0; i < NUM_OSC * UNISON; ++i) {
            float freq = (freq_arr[(i)] + 1) * n_gain;
            float freq_hz = fastpow2(freq);
            freq_hz = std::max(.020f, std::min(freq_hz, 9000.f));
            phase_inc[(i)] = (freq_hz) / (float)OS_RATE;
            if (i == OSC_1 * UNISON || i == OSC_2 * UNISON) {
                freq_store[i / UNISON] = freq_hz;
            }
        }
        for (int i = 0; i < NUM_OSC; ++i) {
            float sh_osc = shape_arr[(i)];
            const float sl = sh_osc;
            sh_osc = sl > 0 ? (1 - ((sl - 1) * (sl - 1))) : ((sl + 1) * (sl + 1) - 1);
            sh_osc = (sh_osc / 2) + 0.5;
            sh_osc = std::min(1.f - (freq_store[i] / fmax), std::max(sh_osc, freq_store[i] / fmax));
            gain_down[(i)] = 1.f / (1 - sh_osc);
            gain_up[(i)] = 1.f / sh_osc;
            shape[(i)] = sh_osc;
        }
    }

    inline void runOsc(const std::array<float, NUM_OSC * UNISON> &phase_mod) {

        std::array<bool, UNISON> sync_arr;
        // osc A code

        for (int i = 0; i < UNISON; ++i) {
            float phase_tmp = phase[i] + phase_inc[i];
            const bool wrap = phase_tmp > 1.0f;
            sync_arr[i] = wrap && sync;
            phase[i] = std::fmod(phase_tmp, 1.0f);
            phase_tmp = std::fmod(phase_tmp + phase_mod[i] + 20.f, 1.0f);
            phase_tmp > shape[OSC_1] ? (osc_tri[i] = (phase_tmp - ((shape[OSC_1] + 1) / 2)) * gain_down[(OSC_1)]) : (osc_tri[(i)] = ((shape[(OSC_1)] / 2) - phase_tmp) * gain_up[(OSC_1)]);

            if (!sub) {
                osc_sqr[(i)] = phase_tmp > shape[(OSC_1)] ? -1.0f : 1.f;
            } else {
                if (wrap) {
                    osc_sqr[(i)] = (osc_sqr[(i)] > 0 ? -1.0 : 1.0);
                }
            }
            xor_out[(i)] = osc_sqr[(i)] * osc_sqr[(i + UNISON)];
        }
        // osc B code
        for (int i = 2; i < UNISON * (OSC_2 + 1); ++i) {
            float phase_tmp = phase[(i)] + phase_inc[(i)];
            phase[(i)] = sync_arr[(i - UNISON)] ? 0 : std::fmod(phase_tmp, 1.0f);
            phase_tmp = std::fmod(phase_tmp + phase_mod[i] + 20.f, 1.0f);

            phase_tmp > shape[(OSC_2)] ? (osc_tri[(i)] = (phase_tmp - ((shape[(OSC_2)] + 1) / 2)) * gain_down[(OSC_2)]) : (osc_tri[(i)] = ((shape[(OSC_2)] / 2) - phase_tmp) * gain_up[(OSC_2)]);
            osc_sqr[(i)] = phase_tmp > shape[(1)] ? -1.0f : 1.f;
        }
    }
};

class MoniqueVoiceModel {

  public:
    // Constants
    MoniqueVoiceModel(uint voice_num, VoiceData &data, WavetableOsc &wavetable) :
        _voice_num(voice_num),
        _input(data),
        _wavetable(wavetable) {
        _analog_voice.fxMuteLoopL(false);
        _analog_voice.fxMuteLoopR(false);
        _noise_osc.setVolume(1.0);
    }

    ~MoniqueVoiceModel(){};

    void reloadCal() {
        _analog_voice.reloadCal();
    }

    bool blacklisted() { return false; };

    void setMuteDisable(bool mute_dis) {
    }

    void setBlacklist(bool bl) {
    }

    void reEvaluateMutes(){};

    void loadVoiceCalibration(){};

    void saveVoiceCalibration() {}

    inline void getNoiseSample(const NoiseModes &mode, const ProcessAudioData &audio_data, const float &ex_in_gain, const float &level, const int &i, const float &unison_1, const float &unison_2, float &noise_sample) {
        switch (mode) {
        case NoiseModes::WHITE:
            noise_sample = (_noise_osc.getSample() - 0.5) * level * unison_1;
            noise_sample += (_noise_osc.getSample() - 0.5) * level * unison_2;
            break;
        case NoiseModes::EXT_L:
            noise_sample = (*audio_data.input_1)[i] * level * ex_in_gain * (unison_1 + unison_2);
            break;
        case NoiseModes::EXT_R:
            noise_sample = (*audio_data.input_2)[i] * level * ex_in_gain * (unison_1 + unison_2);
            break;
        case NoiseModes::LOOP:
            noise_sample = 0;
        default:
            break;
        }
    };

    float phase = 0;
    int counter = 0;

    void generateVoiceBuffers(ProcessAudioData &audio_data) {
        auto &audio = *audio_data.voice_high_res.at(_voice_num);
        auto &cv = *audio_data.voice_cv.at(_voice_num);

        const auto filter_mode = _input.filter_mode;
        bool print = false;

        const bool od = _input.overdrive;
        bool sub = _input.sub_osc;
        _analog_osc.sync = _input.hard_sync;
        _analog_osc.sub = sub;

        const NoiseModes noise_mode = _input.noise_mode;
        const float &fm_1_1 = *_input.fm_1_1;
        const float &fm_1_2 = *_input.fm_1_2;
        const float &fm_2_1 = *_input.fm_2_1;
        const float &fm_2_2 = *_input.fm_2_2;

        const float dig_offset = _od_on ? _digital_voice_offset_od : _digital_voice_offset;

        float cv_osc_1, cv_osc_2;
        for (int cv_i = 0; cv_i < CV_BUFFER_SIZE; ++cv_i) {
            float osc_mod_osc_1_spl;
            float osc_mod_osc_2_spl;
            // unison
            float osc_1_unison = _input.unison_1.at(cv_i);
            float osc_2_unison = _input.unison_2.at(cv_i);

            // cv loop
            const float shape_osc_1 = _input.osc_1_shape[cv_i];
            const float shape_osc_2 = _input.osc_2_shape[cv_i];
            const float shape_osc_3 = _input.osc_3_shape[cv_i];

            // tone calc
            float tone = _input.osc_4_tone_tmp;
            _tone_filter.update(tone);
            const float pwm_osc_1 = cv_clip(-shape_osc_1);
            float vt = (shape_osc_1);
            const float pwm_osc_2 = cv_clip(-(shape_osc_2));
            const float unison_offset = _input.unison_offset;
            const float osc_1_freq = (_input.osc_1_pitch[cv_i]);
            const float osc_2_freq = _input.osc_2_pitch[cv_i];
            const float osc_3_freq = _input.osc_2_pitch[cv_i];

            // calcualate osc inc
            _analog_osc.recalculate({
                                        osc_1_freq,
                                        osc_1_freq + unison_offset,
                                        osc_2_freq,
                                        osc_2_freq + unison_offset,
                                    },
                {pwm_osc_1, pwm_osc_2});
            float level_osc_1_sqr, level_osc_2_sqr, level_osc_1_tri, level_osc_2_tri;
            level_osc_1_sqr = cv_clip(_input.sqr_1_lev[cv_i] * .61);
            level_osc_1_tri = cv_clip(_input.tri_1_lev[cv_i] * 1);
            level_osc_2_sqr = cv_clip(_input.sqr_2_lev[cv_i] * 0.61);
            level_osc_2_tri = cv_clip(_input.tri_2_lev[cv_i] * 1);
            const float level_left = cv_clip(_input.vca_l[cv_i]);
            const float level_right = cv_clip(_input.vca_r[cv_i]);
            const float level_xor = cv_clip(_input.xor_lev[cv_i]) * 0.302 * (float)(noise_mode == NoiseModes::XOR);
            const float noise_level = cv_clip(_input.xor_lev[cv_i]) * 0.61 * ((float)(noise_mode != NoiseModes::XOR));

            const float filter_cut_freq = fastpow2(5.1 * (2 * cv_clip(_input.filt_1_cut[cv_i])) + 4.5 + 0.554588852);
            const float filter_res = cv_clip(_input.filt_1_res[cv_i]);
            _highpass.Update(filter_cut_freq, filter_res);

            if (_voice_num == 0) {
                // printf("\n %f", level_osc_2_sqr);
            }

            bool print = false;
            float tmp = 0;
            std::array<float, CV_TO_SR_MULT * 4> osr_buffer;
            for (int i = 0; i < CV_TO_SR_MULT; ++i) {
                // 96k audio loop
                const int i_96 = cv_i * CV_TO_SR_MULT + i;
                float wt_sample = _input.wt_buffer[(i_96)];
                float noise_sample = 0;
                getNoiseSample(noise_mode, audio_data, _input.aux_level, noise_level, i_96, osc_1_unison, osc_2_unison, noise_sample);
                float loop_level = 0;
                if (noise_mode == NoiseModes::LOOP) {
                    loop_level = noise_level * 3;
                }
                noise_sample += loop_level * _prev_voice_output;
                const float wt_inc = (_prev_wt_sample - wt_sample) / OS_RATE;
                const float noise_inc = (_prev_noise_sample - noise_sample) / OS_RATE;
                float sample = 0;
                for (int j = 0; j < OS_MULT; ++j) {
                    // oversampled loop

                    phase_mod_0_a = (fm_2_1 * _analog_osc.osc_tri[2] * 10);
                    phase_mod_1_a = (fm_2_1 * _analog_osc.osc_tri[3] * 10);
                    phase_mod_0_b = (fm_1_2 * _analog_osc.osc_tri[0] * 10);
                    phase_mod_1_b = (fm_1_2 * _analog_osc.osc_tri[1] * 10);
                    _analog_osc.runOsc({phase_mod_0_a, phase_mod_1_a, phase_mod_0_b, phase_mod_1_b});
                    float xor_out = (_analog_osc.xor_out[0] * osc_1_unison + _analog_osc.xor_out[1] * osc_2_unison);
                    float mixer_output = xor_out * level_xor + _prev_noise_sample + j * noise_sample;
                    mixer_output = _tone_filter.run(mixer_output);
                    const float t_1_sum = (_analog_osc.osc_tri[0] * osc_1_unison + _analog_osc.osc_tri[1] * osc_2_unison) / 2;
                    const float s_1_sum = (_analog_osc.osc_sqr[0] * osc_1_unison + _analog_osc.osc_sqr[1] * osc_2_unison) / 2;
                    float sum_osc_1 = (t_1_sum)*level_osc_1_tri + (s_1_sum) * (level_osc_1_sqr);

                    const float t_2_sum = (_analog_osc.osc_tri[2] * osc_1_unison + _analog_osc.osc_tri[3] * osc_2_unison) / 2;
                    const float s_2_sum = (_analog_osc.osc_sqr[2] * osc_1_unison + _analog_osc.osc_sqr[3] * osc_2_unison) / 2;
                    float sum_osc_2 = (t_2_sum)*level_osc_2_tri + (s_2_sum) * (level_osc_2_sqr);

                    // update the mod variables
                    if (i == 0 && j == 0) {
                        cv_osc_1 = _analog_osc.osc_tri[0] * _input.osc_1_blend_tri + _analog_osc.osc_sqr[0] * _input.osc_1_blend_sqr;
                        cv_osc_2 = _analog_osc.osc_tri[2] * _input.osc_2_blend_tri + _analog_osc.osc_sqr[2] * _input.osc_2_blend_sqr;
                    }

                    mixer_output += sum_osc_1 + sum_osc_2;
                    mixer_output += _prev_wt_sample + j * wt_inc;
                    sample = aa_filter_1.applyfilter((mixer_output));
                    sample = _highpass.Process(sample);
                    sample = aa_filter_2.applyfilter(sample);
                }
                _prev_voice_output = sample * _input.post_hp_gain + dig_offset;
                (audio)[i_96] = _prev_voice_output;
                _prev_wt_sample = wt_sample;
                _prev_noise_sample = noise_sample;
            }
            _osc_1_cv.at(cv_i) = cv_osc_1;
            _osc_2_cv.at(cv_i) = cv_osc_2;
        }

        _analog_voice.run(cv);
    };

    void printStuff() {

        _print = true;
    }

    float getFilterCutOut() {
        return 0;
    }

    float getAveOscLevel() {
        return 0;
    }

    void setAllocatedToVoice(bool val) {
    }

    void set_filter_offset(float offset) {
        offset = 0.22;
        _filter_offset = offset * 0.2 - 0.1;
    }

    void setMute(bool mute) {
        _analog_voice.setMute(mute);
    }

    auto &getCvOsc1() {
        return _osc_1_cv;
    }

    auto &getCvOsc2() {
        return _osc_2_cv;
    }

  private:
    const uint _voice_num = 0;
    float _prev_wt_sample = 0;
    float _prev_noise_sample = 0;
    bool _sub_oscillator = false;
    bool _drive_boost = false;
    bool _voice_allocated_to_layer = false;
    bool _print = true;
    float _offset_osc_1 = 0;
    float _offset_osc_2 = 0;
    float _hp = 0;
    float _filter_offset = 0;
    float _prev_voice_output = 0;
    float phase_mod_0_a = 0;
    float phase_mod_1_a = 0;
    float phase_mod_0_b = 0;
    float phase_mod_1_b = 0;
    VoiceData &_input;
    NoiseOscillator _noise_osc;
    ToneFilter _tone_filter;
    PinkNoiseGen _pink_gen;
    WavetableOsc _wavetable;
    MoniqueAnalogModel _analog_voice = MoniqueAnalogModel(_voice_num, _input);
    float &_digital_voice_offset = _analog_voice.getOffset();
    float &_digital_voice_offset_od = _analog_voice.getOffsetOd();
    bool &_od_on = _analog_voice.getOdState();
    MoniqueFilter _filter;
    SimpleStateVariableFilter _highpass;
    UnisonMoniqueOsc _analog_osc;
    MoogAAFilter aa_filter_1, aa_filter_2, aa_filter_3;
    bool _series_connection = true;
    bool _stereo = false;
    float buf[128 * os_rate];
    std::array<float, CV_BUFFER_SIZE> _osc_1_cv, _osc_2_cv;

    // when voices are disabled the osc will run at appx ~500hz and with a squarish shape so they can still be tracked
};

} // namespace Monique
