#pragma once
#include "MoniqueCommon.h"
#include "SynthMath.h"

namespace Monique {

static constexpr int DELAY_SIZE = (int)(2.0 * SAMPLE_RATE);
static constexpr int TEMPO_SYNC_STEPS = TempoSyncMultipliers::NUM_SYNC_TEMPOS;
static constexpr float TIME_SLEW_CLIP = .35;
static constexpr float TAU = 1 / (2 * M_PI * 50.f);
static constexpr float DT = 1.f / (float)SAMPLE_RATE;
static constexpr float HP_A = 1.f - std::exp(-DT / TAU);

static int delay_counter = 0;

constexpr int NUM_DELAY_STAGES = 16;
constexpr float DELAY_FRAC_LEN_MULT = 1.f / (float)NUM_DELAY_STAGES;
constexpr int DELAY_LENGTH = DELAY_SIZE * DELAY_FRAC_LEN_MULT;

class FractionalDelay {

  public:
    FractionalDelay() {
        delay_line.reserve(DELAY_LENGTH);
        for (int i = 0; i < DELAY_LENGTH; i++) {
            delay_line.push_back(0.f);
        }
    }

    float runNoInterp(float sample_in) {
        delay_line.at(write_ptr) = sample_in;
        write_ptr++;
        write_ptr = write_ptr == DELAY_LENGTH ? 0 : write_ptr;
        const float sample_out = delay_line.at(read_ptr++);
        read_ptr = read_ptr == DELAY_LENGTH ? 0 : read_ptr;
        return sample_out;
    }

    float runInterpolate(float sample_in) {

        delay_line.at(write_ptr) = sample_in;
        float pos = (float)write_ptr - _length;
        pos = pos < 0.f ? pos += DELAY_LENGTH : pos;
        int pos_1 = (int)std::floor(pos);
        pos_1 = pos_1 >= DELAY_LENGTH ? pos_1 - DELAY_LENGTH : pos_1;
        const int pos_2 = pos_1 + 1 >= DELAY_LENGTH ? 0 : pos_1 + 1;
        float frac = pos - std::floor(pos);
        write_ptr++;
        write_ptr = write_ptr == DELAY_LENGTH ? 0 : write_ptr;
        const float sample_1 = delay_line.at((pos_1));
        const float sample_2 = delay_line.at(pos_2);
        return sample_1 + (sample_2 - sample_1) * frac;
    }

    void setLength(float length) {
        _length = length;
    }

    std::vector<float> delay_line;
    int read_ptr = 0;
    int write_ptr = 0;
    float _length = 0;
};

class BBDDelay {
  public:
    BBDDelay() {
        _delay_counter = delay_counter++;
    }

    ~BBDDelay() {}

    void mute() {
    }

    void setTime(float time, bool adj = false) {
        if (!adj) {
            _set_delay_time_setting = time;
        }
        _delay_time_setting = time;
    }

    void adjTime(float adj) {
        setTime(std::clamp(_set_delay_time_setting + adj, 0.f, 1.f), true);
    }

    void setTimeSync(float time) {
        _delay_time_setting_sync = time;
    }

    void setTempoSync(bool temposync) {
        _tempo_sync = temposync;
    }

    void setTempo(float tempo) {
        _tempo = tempo;
    }

    void reset() {
        for (auto &delay : _delay_buffers) {
            for (int i = 0; i < DELAY_LENGTH; i++) {
                delay.delay_line.at(i) = 0.f;
            }
        }
    }

    void time_skip() {
        _time_skip = true;
    }

    void setFB(float fb, bool adj = false) {
        if (!adj) {
            _set_feedback = fb;
        }
        _feedback = fb;
    }

    void adjFB(float adj) {
        _feedback = std::clamp(_set_feedback + adj, 0.f, 1.f);
    }

    void setFilterMix(float mix) {
        _filter_mix = mix;
    }

    void setFilter(float cf, bool adj = false) {
        if (!adj) {
            _set_cf = cf;
        }
        float cutoff = fastpow2(5 + cf * 9);
        float tau = 1.f / (2.f * M_PI * cutoff);
        constexpr float dT = 1.f / (float)SAMPLE_RATE;
        float a = 1 - expf32(-dT / tau);
        float b = 1.f - a;
        _filt_a0 = a;
        _filt_b1 = b;
    }

    void adjFilter(float adj) {
        setFilter(std::clamp(_set_cf + adj, 0.f, 1.f), true);
    }

    void setLFORate(float rate) {
        _lfo_rate = rate;
        _lfo_phase_inc = 628.31853f * (float)fastpow2(3.4 * ((float)(3.0f * _lfo_rate - 2.0f))) / (float)SAMPLE_RATE;
    }

    void setLFOGain(float gain) {
        _lfo_amount_set = gain;
    }

    void run(float *(input), float *output) {
        calc_time();
        _lfo_amount = 100 * _lfo_amount_set * (0.2 + 0.2 * _time_setting);
        _gain_smooth = _gain;

        // find the size of the linear step for the loop increment
        float wet = _wet_mix;
        _lfo_phase += _lfo_phase_inc;
        _lfo_phase = std::fmod(_lfo_phase, 2 * M_PI);

        for (int buf_i = 0; buf_i < BUFFER_SIZE; buf_i++) {

            // calc delay time
            const float lfo = _lfo_amount * (float)fastsin(_lfo_phase - M_PI);
            float next_read_pos = _time + lfo;
            // slew clip delay time
            float diff = next_read_pos - _read_position;
            if (!_time_skip) {
                diff = diff > TIME_SLEW_CLIP ? diff = TIME_SLEW_CLIP : diff;
                diff = diff < -TIME_SLEW_CLIP ? diff = -TIME_SLEW_CLIP : diff;
            }
            _time_skip = false;
            _read_position += diff;
            for (int d = 0; d < NUM_DELAY_STAGES; d++) {
                _delay_buffers.at(d).setLength(_read_position / (float)NUM_DELAY_STAGES);
            }

            float delay_in = input[buf_i] * _gain_smooth + _feedback * 1.05 * _prev_sample;
            _highpass_state += HP_A * (delay_in - _highpass_state);
            delay_in -= _highpass_state;
            const float filter_out = _filt_a0 * delay_in + _filt_b1 * _filt_Y1;
            _filt_Y1 = filter_out;
            delay_in = fasttanh(2 * filter_out) / 2;
            float delay_sample = delay_in;
            for (int d = 0; d < NUM_DELAY_STAGES; d++) {
                delay_sample = _delay_buffers.at(d).runInterpolate(delay_sample);
            }

            _prev_sample = delay_sample;

            // sum output into output buffer
            output[buf_i] = delay_sample * 2;
        }
    }

    void calc_time() {
        if (_tempo_sync) {
            float time_multiplier = _tempo_settings.getTempoSyncMultiplier(_delay_time_setting_sync);
            float time = 1.f / (_tempo * time_multiplier);
            float time_samples = time * SAMPLE_RATE;
            time_samples = time_samples > ((float)DELAY_SIZE * 0.8) ? ((float)DELAY_SIZE * 0.8f) : time_samples;
            time_samples = time_samples < 1 ? 1 : time_samples;
            _time = time_samples;
        } else {
            _time_setting = _delay_time_setting;
            float time = _time_setting;
            constexpr float time_offset = 0.18;
            constexpr float time_gain = 1;
            time = time_offset + (time_gain - time_offset) * time;
            _time = time * time * (float)DELAY_SIZE * 0.8;
        }
    }

  private:
    int _delay_counter;
    static constexpr float cutoff_frequency = 200.0;
    static constexpr float highpass_gain = cutoff_frequency / (float)(2.f * M_PI * SAMPLE_RATE);
    static constexpr float max_gain = 2.5;
    float _highpass_state = 0;
    float _filt_Y1 = 0;
    float _filt_a0 = 0;
    float _filt_b1 = 0;
    float _prev_sample = 0;
    float _read_position = 0;
    float _clipper_amp_dB = 8.6562;
    float _baseline_threshold_dB = -9.0;
    float _clipper_linear_coeff = 1.017;
    float _clipper_squared_coeff = -0.025;
    float _limit_dB = -6.9;
    float _filter_mix = 1;
    float _threshold_dB = _baseline_threshold_dB + _limit_dB;
    float _time_setting = 0;
    float _time = 0.8;
    float _feedback = .5;
    float _feedback_tone = 0.5;
    float _lfo_amount = 0.1;
    float _lfo_amount_set = 0;
    float _lfo_rate = 0.5;
    float _lfo_phase_inc = 0;
    float _lfo_phase = 0;
    float _loop_buffer = 0;
    float _env = 0;
    float _rel = 0;
    float _gain = 2.50;
    float _gain_smooth = 0;
    float _wet_mix = 1;
    int _buffer_ptr = 0;
    bool _tempo_sync = false;
    bool _time_skip = false;
    float _tempo = 1.f;
    float _delay_time_setting = 0;
    float _delay_time_setting_sync = 0;
    float _set_delay_time_setting = _delay_time_setting;
    float _set_feedback = _feedback;
    float _set_cf = 0;
    float _set_level = 0;
    TempoSyncMultipliers _tempo_settings;
    std::array<FractionalDelay, NUM_DELAY_STAGES> _delay_buffers;
};

} // namespace Monique