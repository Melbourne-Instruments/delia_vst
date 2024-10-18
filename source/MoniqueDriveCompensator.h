

#include "SynthMath.h"
#include "common.h"

#pragma once

namespace Monique {

volatile static float fac = 0.3;

class DriveCompensator {
  public:
    DriveCompensator(bool &overdrive, float &compresson_signal, float &main_out_max_level, float &patch_vol, bool &od_en_out) :
        _overdrive(overdrive), _compression_signal(compresson_signal), _main_out_max_level(main_out_max_level), _patch_volume(patch_vol), _od_en_out(od_en_out) {
    }

    ~DriveCompensator(){};

    float &getMixbusCompensation() {
        return _main_vca_gain;
    }

    float &getMixerVcaCompensation() {
        return _mix_gain;
    }

    float *getdriveInput() {
        return &_drive_level;
    }

    void reCalculate() {

        constexpr float od_gain_db = 19.2f;
        constexpr float drive_comp_ratio = 2.5f;

        // range of the drive knob
        static constexpr float drive_range_db = 10.f;
        static constexpr float od_gain = std::pow(10.f, od_gain_db / 20.f);
        static constexpr float od_min_drive = std::pow(10.f, -od_gain_db / 20.f);
        static constexpr float drive_range = std::pow(10.f, drive_range_db / 20.f) / 2;
        static constexpr float od_comp_start = 0.52;
        constexpr float od_actual_min_gain = (od_min_drive) * (1 / (1 + od_min_drive));

        // experimentally found drive comp factor for main output vca
        static constexpr float compensation_factor = .3;
        static constexpr float min_drive_gain = std::pow(10.f, (-drive_range_db) / 20.f);
        static constexpr float max_vol = 0.6f;
        constexpr float factor = 2.5f;
        constexpr float actual_min_gain = (min_drive_gain) * (1 / (1 + min_drive_gain));

        float drive_level = unipolarClip(_drive_level);

        // calculate the level reduction for main out mix muting
        float thresh = 0.02;
        // thresh = fac;
        float mix_mult = 1.0;
        if (_main_out_max_level < thresh) {
            mix_mult = mix_mult * (1.f / thresh) * _main_out_max_level;
        }
        float comp_factor;

        // OD condition when drive is greater than .5, enable drive. perform compensation and adjust main mix output level
        if (drive_level > 0.5) {
            _od_en_out = true;
            float drive_level_od = drive_level * 2 - 1;

            comp_factor = od_comp_start - drive_level_od * .3;
            _pre_compression_vca_gain = comp_factor;
            _pre_compression_vca_gain *= _patch_volume;
            _mix_gain = (od_min_drive + (1 - od_min_drive) * drive_level_od);

        }

        // OD off condition
        else {
            drive_level = drive_level * 2;
            _od_en_out = false;
            comp_factor = 1 - drive_level / factor;
            _pre_compression_vca_gain = comp_factor;
            _pre_compression_vca_gain *= _patch_volume;
            float total_gain = min_drive_gain + (drive_range) * (drive_level);
            _mix_gain = (min_drive_gain + (1 - min_drive_gain) * drive_level);
        }
        // duck the mix VCA's to improve nulling/bleed on main outputs
        _mix_gain *= mix_mult;
        if (_print) {
            // printf("\ncomp: %f %f %f", mix_mult, _mix_gain, _main_out_max_level);
        }
    }

    void print() {
        _print = true;
    }

    void run() {
        _main_vca_gain = _pre_compression_vca_gain - _compression_signal;
    }

  private:
    // amount of gain applied by the overdrive button
    bool &_overdrive;
    float &_compression_signal;
    float &_main_out_max_level;
    float &_patch_volume;

    float _pre_compression_vca_gain = 1;

    float _main_vca_gain = 1;
    float _mix_gain = 1;
    float _drive_level = 0;
    bool _print = false;
    bool &_od_en_out;
};

} // namespace Monique