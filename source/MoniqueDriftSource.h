
#pragma once

#include "MoniqueCommon.h"
#include "SynthMath.h"
#include <array>

constexpr float DRIFT_DC_FILTER = 4000;
constexpr int OSC_SOURCES = 24;

namespace Monique {
class DriftSource {

  public:
    DriftSource(){};

    ~DriftSource(){};

    inline void generateOscDrift(float &sample) {
        const float noise = int(juicy_hash(_rnj_state++) & 0x7fffffff) * (1.f / 0x7fffffff) - 0.5;
        sample += noise;
        sample -= sample / DRIFT_DC_FILTER;
    }

    void recalculate() {
        for (int i = 0; i < OSC_SOURCES; i++) {
            generateOscDrift(_osc_detune_sources.at(i));
        }
    };

    float &getOscSource() {
        _osc_detune_sources.at(len) = 0.f;
        return _osc_detune_sources.at(len++);
    }

    int len = 0;
    int _rnj_state = static_cast<int>(*("DELIA"));
    std::array<float, OSC_SOURCES> _osc_detune_sources;
};
}; // namespace Monique