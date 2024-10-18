#include "FastApprox.h"
#include "MoniqueCommon.h"
#include "common.h"
#include <array>
#include <cmath>

#pragma once

namespace Monique {

inline float pitchShiftMultiplier(float dPitchShiftSemitones) {
    if (dPitchShiftSemitones == 0)
        return 1.0f;

    // 2^(N/12)
    //	return fastPow(2.0, dPitchShiftSemitones/12.0);
    return std::exp2f(dPitchShiftSemitones / 12.0f);
}

[[gnu::noinline]] static bool floatIsValid(float in) {
    constexpr float MAX_SIGNAL_FLOAT = 2.f;
    constexpr float MIN_SIGNAL_FLOAT = -2.f;

    bool res = in < MAX_SIGNAL_FLOAT && in > MIN_SIGNAL_FLOAT;
    if (!res) {
        // printf("\n %f \n", in);
    }
    return res;
}

inline float linearToLogPot(float val) {
    val = std::log2(val);
    return std::max(-200.f, val);
}

inline float logPotToLinear(float val) {
    val = std::exp2f(val);
    return val;
}

constexpr float NOTE_GAIN = 10.f;

inline constexpr float midiNoteToCv(int midi_n) {
    constexpr float cv_scale = NOTE_GAIN / 120.0;
    constexpr float cv_offset = -(60.f * cv_scale);
    return ((float)(midi_n)) * cv_scale + cv_offset;
}

inline float midiNoteToFreq(int midi_n) {
    return 440.f * fastpow2(((float)midi_n - 69.0f) / 12.f);
}

inline float midiSourceScaleBipolar(float value) {
    return value * 2.0f - 1.0f;
}

inline float midiSourceScaleUnipolar(float value) {
    return value * 2.0f;
}

// Special smoothing function for MPE data, it has a very low cutoff but will speed up with a bigger delta to avoid 'lag'
inline float mpeParamSmooth(float value, float &prev_value) {
    constexpr float filter_rate = 50.f / 750.f;
    constexpr float nonlinear_rate = 0.02;
    const float diff = value - prev_value;
    prev_value += ((diff) * (filter_rate + nonlinear_rate * std::fabs(diff)));
    return prev_value;
}

inline float mpeParamFall(float value, float &prev_value, float fall_time) {

    // if value is falling, then decay the value according to fall time
    if (((value) < (prev_value))) {
        prev_value += fall_time * (value - prev_value);
    } else {
        prev_value = value;
    }
    return value;
}

float inline xSquaredParamCurve(float value) {
    return value * value;
};

float inline xAbsParameterCurve(float value) {
    value = (value - .5) * 2;
    return value * std::abs(value);
};

inline float mpePbRangeConvert(uint semitones) {
    return 2.f * (((float)semitones) / 12.0) / NOTE_GAIN;
}

constexpr float MASTER_DETUNE_RANGE_SEMI_T = 2.f;

inline float calcMasterDetune(float param) {
    const float detune_semi_tones = (param - 0.5) * 2.0 * MASTER_DETUNE_RANGE_SEMI_T;
    return (detune_semi_tones / 12.0) / NOTE_GAIN;
}

inline float fastLog2(float x) {
    return fastlog2(x);
    // return std::log2(x);
    // return x;
}

inline float fastlog10(float num) { return std::log10(num); }

inline float volume_knob_cal(float input) {
    float output;
    output = (10.0f / 9.0f) * (exp10f32(input - 1.0f) - 1.0f / 10.0f);
    return output;
};

inline float fastCos(float input) { return std::cos(input); }

inline float fastSin(float input) { return std::sin(input); }

/**
 * @brief WORLDS WORST REALLY UNSAFE RINGBUFFER pls dont use. has an offset (for
 * filter with delay) and inits all values to 0
 *
 * @tparam N size of the ringbuffer
 * @tparam O starting offset point
 */

template <int N, int O>
class delayLineBuffer {
  public:
    delayLineBuffer() {
        for (float val : _array) {
            val = 0;
        }
    };

    ~delayLineBuffer() = default;
    ;

    void run(float number) {
        _array[_tail++] = number;
        _head++;
        _head = (_head > _max_size) ? _head : 0;
        _tail = (_tail >= _max_size) ? _tail : 0;
    }

    float read() {
        float val = _array[_head];

        return val;
    }

    int getSize() {
        int size;
        if (_head >= _tail) {
            size = _head - _tail;
        } else {
            size = _max_size + _head - _tail;
        }

        return size;
    }

  private:
    static constexpr const size_t _max_size = N;
    std::array<float, _max_size> _array;
    size_t _head = 0;
    size_t _tail = 0;
};

// 32 bit fast rng https://github.com/skeeto/hash-prospector
inline u_int32_t juicy_hash(u_int32_t x) {
    x ^= x >> 17;
    x *= 0xed5ad4bb;
    x ^= x >> 11;
    x *= 0xac4c1b51;
    x ^= x >> 15;
    x *= 0x31848bab;
    x ^= x >> 14;

    return x;
};

inline constexpr float driveAdjust(float drive) {
    constexpr float SWITCH_POINT = 0.622;
    constexpr float OUT_SWICH_POINT = 0.5;
    if (drive > SWITCH_POINT) {
        return (drive - SWITCH_POINT) * (OUT_SWICH_POINT / (1.f - SWITCH_POINT)) + OUT_SWICH_POINT;
    } else {
        return drive * (OUT_SWICH_POINT / SWITCH_POINT);
    }
};

static inline float param_smooth(const float input, float &filter_state) {
    float diff = input - filter_state;
    constexpr float DIFF_CLIP_L = 0.003;
    constexpr float DIFF_CLIP_S = 0.001;
    constexpr float DIFF_THRESH = 0.1;
    constexpr float PARAM_SMOOTH_COEFF = 15 / (float)BUFFER_RATE;
    float diff_c;
    if (std::abs(diff) > DIFF_THRESH) {

        diff_c = std::clamp(diff, -DIFF_CLIP_L, DIFF_CLIP_L);
    } else {

        diff_c = std::clamp(diff, -DIFF_CLIP_S, DIFF_CLIP_S);
    }
    filter_state += diff * PARAM_SMOOTH_COEFF + diff_c;
    return filter_state;
}

/**
 * @brief Clip the signal to the range 1.0->-1.0
 *
 * @param cv
 * @return float
 */

inline float cv_clip(float cv) {
    return std::max(std::min(cv, 1.0f), -1.0f);
};

inline float unipolarClip(float cv) {
    cv = (cv > 1.0f) ? 1.0f : cv;
    cv = (cv < -.0f) ? -.0f : cv;
    return cv;
};

class BiquadFilter {
  public:
    BiquadFilter(float a1, float a2, float b0, float b1, float b2);

    ~BiquadFilter() {}

    inline void reset() {
        // Reset the filter
        _m1 = 0;
        _m2 = 0;
        _dn = 1e-20f;
    }

    inline float process(float input) {
        float w = input - (_a1 * _m1) - (_a2 * _m2) + _dn;
        float out = (_b1 * _m1) + (_b2 * _m2) + (_b0 * w);
        _m2 = _m1;
        _m1 = w;
        _dn = -_dn;
        return out;
    }

  private:
    const float _a1, _a2, _b0, _b1, _b2;
    float _m1, _m2;
    float _dn;
};

inline int16_t _floatSampleToInt16(float sample) {
    sample = sample > 0.99 ? 0.99 : sample;
    sample = sample < -0.99 ? -0.99 : sample;
    int16_t val = std::round(sample * (std::numeric_limits<std::int16_t>::max()));
    return val;
};
#ifndef __x86_64__
inline void setFpStatusRegister(intptr_t fpsr) noexcept {

    asm volatile("msr fpcr, %0"
                 :
                 : "ri"(fpsr));
}

inline intptr_t getFpStatusRegister() noexcept {
    intptr_t fpsr = 0;

    asm volatile("mrs %0, fpcr"
                 : "=r"(fpsr));

    return fpsr;
}

inline void enableFlushToZeroMode(bool shouldEnable) noexcept {

    intptr_t mask = (1 << 24 /* FZ */);

    setFpStatusRegister((getFpStatusRegister() & (~mask)) | (shouldEnable ? mask : 0));
}
#else
inline void enableFlushToZeroMode(bool shouldEnable) {}
#endif

} // namespace Monique