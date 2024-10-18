/**
 * @file AnalogFiltGen.h
 * @brief Declaration of the analog filter calibrator class
 * @date 2022-07-06
 *
 * Copyright (c) 2023 Melbourne Instruments
 *
 */

#pragma once

#include "MoniqueCommon.h"
#include "SynthMath.h"
#include <array>
#include <cmath>
#include <stdint.h>
#include <vector>

namespace Monique {

#define MOOG_E 2.71828182845904523536028747135266250
#define MOOG_LOG2E 1.44269504088896340735992468100189214
#define MOOG_LOG10E 0.434294481903251827651128918916605082
#define MOOG_LN2 0.693147180559945309417232121458176568
#define MOOG_LN10 2.30258509299404568401799145468436421
#define MOOG_PI 3.14159265358979323846264338327950288
#define MOOG_PI_2 1.57079632679489661923132169163975144
#define MOOG_PI_4 0.785398163397448309615660845819875721
#define MOOG_1_PI 0.318309886183790671537767526745028724
#define MOOG_2_PI 0.636619772367581343075535053490057448
#define MOOG_2_SQRTPI 1.12837916709551257389615890312154517
#define MOOG_SQRT2 1.41421356237309504880168872420969808
#define MOOG_SQRT1_2 0.707106781186547524400844362104849039
#define MOOG_INV_PI_2 0.159154943091895

#define NO_COPY(C)         \
    C(const C &) = delete; \
    C &operator=(const C &) = delete
#define NO_MOVE(C)    \
    NO_COPY(C);       \
    C(C &&) = delete; \
    C &operator=(const C &&) = delete

#define SNAP_TO_ZERO(n)               \
    if (!(n < -1.0e-8 || n > 1.0e-8)) \
        n = 0;

// Linear interpolation, used to crossfade a gain table
inline float moog_lerp(float amount, float a, float b) {
    return (1.0f - amount) * a + amount * b;
}

inline float moog_min(float a, float b) {
    a = b - a;
    a += fabs(a);
    a *= 0.5f;
    a = b - a;
    return a;
}

// Clamp without branching
// If input - _limit < 0, then it really substracts, and the 0.5 to make it half the 2 inputs.
// If > 0 then they just cancel, and keeps input normal.
// The easiest way to understand it is check what happends on both cases.
inline float moog_saturate(float input) {
    float x1 = fabs(input + 0.95f);
    float x2 = fabs(input - 0.95f);
    return 0.5f * (x1 - x2);
}

// Imitate the (tanh) clipping function of a transistor pair.
// to 4th order, tanh is x - x*x*x/3; this cubic's
// plateaus are at +/- 1 so clip to 1 and evaluate the cubic.
// This is pretty coarse - for instance if you clip a sinusoid this way you
// can sometimes hear the discontinuity in 4th derivative at the clip point
inline float clip(float value, float saturation, float saturationinverse) {
    float v2 = (value * saturationinverse > 1 ? 1 : (value * saturationinverse < -1 ? -1 : value * saturationinverse));
    return (saturation * (v2 - (1. / 3.) * v2 * v2 * v2));
}

#define HZ_TO_RAD(f) (MOOG_PI_2 * f)
#define RAD_TO_HZ(omega) (MOOG_INV_PI_2 * omega)

#ifdef __GNUC__
#define ctz(N) __builtin_ctz(N)
#else
template <typename T>
inline int ctz(T x) {
    int p, b;
    for (p = 0, b = 1; !(b & x); b <<= 1, ++p)
        ;
    return p;
}
#endif

inline double fast_tanh(double x) {
    double x2 = x * x;
    return x * (27.0 + x2) / (27.0 + 9.0 * x2);
}

class LadderFilterBase {
  public:
    LadderFilterBase(float sampleRate) :
        sampleRate(sampleRate) {}

    virtual ~LadderFilterBase() {}

    virtual void Process(float *samples, uint32_t n) = 0;
    virtual void SetResonance(float r) = 0;
    virtual void SetCutoff(float c) = 0;

    float GetResonance() { return resonance; }

    float GetCutoff() { return cutoff; }

  protected:
    float cutoff;
    float resonance;
    float sampleRate;
};

/*
Copyright 2012 Stefano D'Angelo <zanga.mail@gmail.com>
Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.
THIS SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

/*
This model is based on a reference implementation of an algorithm developed by
Stefano D'Angelo and Vesa Valimaki, presented in a paper published at ICASSP in 2013.
This improved model is based on a circuit analysis and compared against a reference
Ngspice simulation. In the paper, it is noted that this particular model is
more accurate in preserving the self-oscillating nature of the real filter.
References: "An Improved Virtual Analog Model of the Moog Ladder Filter"
Original Implementation: D'Angelo, Valimaki
*/

// Thermal voltage (26 milliwats at room temperature)

inline double vox_fasttanh2(const double x) {
    const double ax = fabs(x);
    const double x2 = x * x;

    return (x * (2.45550750702956 + 2.45550750702956 * ax + (0.893229853513558 + 0.821226666969744 * ax) * x2) /
            (2.44506634652299 + (2.44506634652299 + x2) *
                                    fabs(x + 0.814642734961073 * x * ax)));
}

struct MoogAAFilter {
    typedef float REAL;
#define NBQ 2
    REAL biquada[4] = {0.5885934246248039, -1.1168900804200317, 0.19773812910557706, -0.8167265333195526};
    REAL biquadb[4] = {1.0000000000000002, 0.595240468137887, 1, 1.660089401681661};
    REAL gain = 52.85208803784947;
    float xyv[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    /**
     * @brief 4th order chyby filter at 0.022 normalised cutoff
     *
     * @param v
     * @return float
     */
  public:
    inline REAL applyfilter(REAL v) {
        int i, b, xp = 0, yp = 3, bqp = 0;
        REAL out = v / gain;
        for (i = 8; i > 0; --i) {
            xyv[i] = xyv[i - 1];
        }
        for (b = 0; b < NBQ; ++b) {
            int len = 2;
            xyv[xp] = out;
            for (i = 0; i < len; ++i) {
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

constexpr int os_rate = 4;

class ImprovedMoog {
  public:
    ImprovedMoog() {
        memset(V, 0, sizeof(V));
        memset(dV, 0, sizeof(dV));
        memset(tV, 0, sizeof(tV));

        drive = 1.0f;
        Update(1000, .1);
    }

    ~ImprovedMoog() {}

#define moog_tanh fast_tanh

    // Template param ADD will sum into the samples_out array rather than overwrite
    template <bool NONLINEAR>
    inline float Process(const float input) {
        constexpr float hp_gain = 0.997;
        float dV0, dV1, dV2, dV3;
        const bool hp = _hp;

        dV0 = -g * (std::tanh((input + offset + resonance * V[3]) / (2.0 * VT)) + tV[0]);
        V[0] += (dV0 + dV[0]) / (2.0 * sampleRate);
        dV[0] = dV0;
        if (NONLINEAR) {
            tV[0] = moog_tanh(V[0] / (2.0 * VT));
        } else {
            tV[0] = (V[0] / (2.0 * VT));
        }

        dV1 = g * (tV[0] - tV[1]);
        V[1] += (dV1 + dV[1]) / (2.0 * sampleRate);
        dV[1] = dV1;
        if (NONLINEAR) {
            tV[1] = moog_tanh(V[1] / (2.0 * VT));
        } else {
            tV[1] = (V[1] / (2.0 * VT));
        }

        dV2 = g * (tV[1] - tV[2]);
        V[2] += (dV2 + dV[2]) / (2.0 * sampleRate);
        dV[2] = dV2;
        if (NONLINEAR) {
            tV[2] = moog_tanh(V[2] / (2.0 * VT));
        } else {
            tV[2] = (V[2] / (2.0 * VT));
        }

        dV3 = g * (tV[2] - tV[3]);
        V[3] += (dV3 + dV[3]) / (2.0 * sampleRate);
        dV[3] = dV3;
        if (NONLINEAR) {
            tV[3] = moog_tanh(V[3] / (2.0 * VT));
        } else {
            tV[3] = (V[3] / (2.0 * VT));
        }

        // printf("\n off %f", offset);

        if (hp) {
            return input * hp_gain + V[3];
        } else {
            return V[3];
        }
    }

    void print() {
        printf("\n filter: %f %f", resonance, VT);
    }

    void filterOffset(float off) {
        offset = off;
    }

    void vt(float vt) {
        VT = vt;
    }

    void Update(float c, float r) {
        if (r > 1.0)
            r = 1.0;
        if (r < 0.001)
            r = 0.001;
        resonance = 6.844399855 * r;
        cutoff = c;
        if (cutoff > 22000) {
            cutoff = 22000;
        }
        if (cutoff < 20)
            x = (MOOG_PI * cutoff) / sampleRate;
        g = 4.0 * MOOG_PI * VT * cutoff * (1.0 - x) / (1.0 + x);
    }

    void hp(bool is_hp) {
        _hp = is_hp;
    }

    bool log = false;

  private:
    bool _hp;
    float V[4];
    float dV[4];
    float tV[4];
    float VT = .4;
    float cutoff;
    float resonance;
    float hp_compensation = 1.0;
    float hp_res_comp = 0.0;
    float _hp_gain = 0;
    static constexpr float sampleRate = SAMPLE_RATE * os_rate;
    float x;
    float g;
    float drive;
    float offset = 0;
};

struct SimpleStateVariableFilter {
    float _band = 0;
    float _high = 0;
    float _low = 0;
    static constexpr float sample_rate = SAMPLE_RATE * os_rate;
    float _cut_coeff_old = 1.f;
    float _cutoff_old = 20;
    float _cut_inc = 0.f;
    float _res_set = .5;

    inline void Update(float cutoff_f, float res) {
        const float cut_coeff_new = cutoff_f / sample_rate * 2.02 * M_PI;
        const float cut_coeff_old = _cutoff_old / sample_rate * 2.02 * M_PI;
        _cut_inc = (cut_coeff_old - cut_coeff_new) / (float)(OS_MULT_MODEL * BUFFER_SIZE);
        _cut_coeff_old = cut_coeff_old;
        _res_set = res;
        _cutoff_old = cutoff_f;
    }

    inline float Process(const float in) {
        _cut_coeff_old += _cut_inc;

        float res = 1.4142135 - 2.0f * _res_set;
        float band2 = fast_tanh(_band); // clip the fb so we can have self resonance
        _low = _low + _cut_coeff_old * band2;
        _high = in - _low - res * band2;
        _band = _cut_coeff_old * _high + band2;

        return fast_tanh(_high);
    }
};

struct SvfHighpass {

    float _k, _A, _Asqrt, _a1, _a2, _a3, _g, _m0, _m1, _m2, _v1, _v2, _v3, _ic1eq, _ic2eq;
    float _cutoff = 20;
    float _resonance = 0;

    inline void update() {

        float q = _resonance + 0.01;
        float g = tan((_cutoff / ((float)SAMPLE_RATE) * M_PI * os_rate));
        float k = 1 / (q * _A);
        _A = 1;
        float _ASqrt = sqrt(_A);
        _a1 = 1 / (1 + g * (g + k));
        _a2 = g * _a1;
        _a3 = g * _a2;
        _m0 = 1;
        _m1 = -k;
        _m2 = -1;
    }

    inline float processSample(float in) {
        _v3 = (in - _ic2eq);
        _v1 = 1 * _a1 * _ic1eq + 1 * (_a2 * _v3);
        _v2 = _ic2eq + _a2 * _ic1eq + _a3 * _v3;
        _ic1eq = (2 * _v1 - _ic1eq);
        _ic2eq = 2 * _v2 - (1 * _ic2eq);

        return _m0 * in + _m1 * _v1 + _m2 * _v2;
    }
};

class MoniqueFilter {
  public:
    void setup(float cut, float cut_2, float res, float res_2, bool hp) {
        _filter_1.hp(false);
        _filter_1.Update(cut, res);
        _filter_2._cutoff = cut_2;
        _filter_2._resonance = res_2;
    }

    void run(float *in, float *out_l, float *out_r) {
        const bool series_connection = _series_connection;
        const bool stereo = _stereo;
        if (series_connection) {

            for (int i = 0; i < BUFFER_SIZE; i++) {
                float sample = in[i];
                float out;
                for (int i_2 = 0; i_2 < os_rate; i_2++) {
                    out = 10 * aa_filter_1.applyfilter(sample);
                    const float lp_out = _filter_1.Process<false>(out);
                    out = aa_filter_3.applyfilter(lp_out);
                    out = aa_filter_2.applyfilter(out);
                    sample = 0.f;
                }
                out_l[i] = out;
            }
        } else {
            for (int i = 0; i < BUFFER_SIZE; i++) {
                float sample = in[i];
                float s_l, s_r;
                for (int i_2 = 0; i_2 < os_rate; i_2++) {
                    const float upsampled_input = 10 * aa_filter_1.applyfilter(sample);
                    const float lp_out_l = _filter_1.Process<false>(upsampled_input);
                    s_l = aa_filter_2.applyfilter(lp_out_l);
                    float lp_out_r;
                    s_r = aa_filter_2.applyfilter(lp_out_r);
                    sample = 0.f;
                }
                out_l[i] = s_l;
                out_r[i] = s_r;
            }
        }
    }

    bool log = false;

  private:
    ImprovedMoog _filter_1;
    SvfHighpass _filter_2;
    MoogAAFilter aa_filter_1, aa_filter_2, aa_filter_3;
    bool _series_connection = true;
    bool _stereo = false;
    float buf[128 * os_rate];
};

} // namespace Monique