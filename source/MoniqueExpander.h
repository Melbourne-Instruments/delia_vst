#pragma once
#include "MoniqueCommon.h"
#include "SynthMath.h"

namespace Monique {

constexpr float THRESHOLD = -6.5;
constexpr float LIN_THRESH = std::pow(10.f, THRESHOLD / 20.f);
constexpr float RATIO = 5.0;
constexpr float MAX_GR = -4.0f;
constexpr float RMS_REL = 0.1;
constexpr float RMS_SAMPLE = RMS_REL * (float)SAMPLE_RATE;
constexpr float RMS_TIME_CONST = std::pow(37.f / 100.f, 1.f / RMS_SAMPLE);
static_assert(RMS_TIME_CONST > 0);
constexpr float DC_FILTER = 0.00002;

class MoniqueExpander {
    float _rms_1 = 0;
    float _rms_2 = 0;
    int counter = 0;
    float _dc_1 = 0;
    float _dc_2 = 0;

  public:
    MoniqueExpander() {
    }

    void process(audioBuffer &input_1, audioBuffer &input_2) {
        float s1, s2;
        for (uint i = 0; i < BUFFER_SIZE; ++i) {
            float &sample_1 = input_1[i];
            float &sample_2 = input_2[i];
            sample_1 -= _dc_1;
            sample_2 -= _dc_2;
            _dc_1 += sample_1 * DC_FILTER;
            _dc_2 += sample_2 * DC_FILTER;
            // calulate the square sum of the 2 channels, then integrate with a simple LPF and square, this approximates the response of a RMS detector
            const float ave = sample_1 * sample_1 + sample_2 * sample_2;
            _rms_1 = ave + RMS_TIME_CONST * _rms_1;
            float rms = sqrtf32(_rms_1);

            // calculate a dB based gain reduction depending on the amount the RMS level is lower than the threshold
            float db_under_thresh = fastlog2(rms / LIN_THRESH);
            db_under_thresh = db_under_thresh > 0 ? 0.f : db_under_thresh;
            db_under_thresh = db_under_thresh < MAX_GR ? MAX_GR : db_under_thresh;
            const float gain_reduction = db_under_thresh * (RATIO - 1.f) / RATIO;
            float gain_linear = fastpow2(gain_reduction);

            // apply the gain reduction
            sample_1 *= gain_linear;
            sample_2 *= gain_linear;
            s1 = gain_linear;
            s2 = rms / LIN_THRESH;
        }
    }
};

} // namespace Monique