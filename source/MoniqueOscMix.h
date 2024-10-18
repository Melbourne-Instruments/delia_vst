

#include "SynthMath.h"
#include "common.h"

#pragma once

namespace Monique {

constexpr float saw_pos = 0.25;
constexpr float square_pos = 0.75;

class OscMixer {
  public:
    OscMixer(float &tri, float &sqr, float &pwm, float &blend_tri, float &blend_sqr) :
        _tri_lev(tri), _sqr_lev(sqr), _pwm(pwm), _blend_tri(blend_tri), _blend_sqr(blend_sqr){};
    ~OscMixer(){};

    static inline float shapeToPwm(const float shape) {
        static constexpr float slope = 3.f;
        static constexpr float offset_2 = 1.f;
        static constexpr float section_1_end = 1.f / 3.f;
        float out;
        if (shape < section_1_end) {
            out = slope * shape;
        } else {
            out = std::abs(-slope * (shape - section_1_end) + 1.f);
        }
        return out;
    }

    static inline void shapeToLevels(float shape, float &sqr, float &tri) {
        static constexpr float m = -1.f / (square_pos - saw_pos);
        static constexpr float c = -m * (square_pos);
        float tmp = m * shape + c;
        tmp = std::max(0.f, std::min(1.f, tmp));
        tri = tmp;
        sqr = 1.0f - tri;
    }

    float *getShapeInput() {
        return &_shape;
    }

    float *getGainInput() {
        return &_gain;
    }

    void reCalculate() {
    }

    void print() {
        printb = true;
    }

    void run() {

        const float gain = _gain;

        // get OSC pwm and sqr/tri mix from the shape control

        _pwm = shapeToPwm(_shape);
        shapeToLevels(_shape, _sqr_lev, _tri_lev);
        _blend_tri = _tri_lev;
        _blend_sqr = _sqr_lev;
        _sqr_lev *= gain;
        _tri_lev *= gain;
    }

  private:
    static constexpr float max_vol = 0.5f;
    float _shape = 0.0f;
    float _gain = 1.0f;
    float &_tri_lev;
    float &_sqr_lev;
    float &_blend_tri;
    float &_blend_sqr;
    float &_pwm;
    bool printb = false;
};

} // namespace Monique