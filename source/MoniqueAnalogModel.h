
#pragma once

#include "MoniqueCommon.h"
#include "SynthMath.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

namespace Monique {

constexpr float MUTE_THRESHOLD = 0.001;

class MoniqueAnalogModel {

  public:
    bool loadCalibration() {
        std::fstream file_handler;
        std::stringstream fpath;
        fpath << cal_path << "voice_" << _voice_n << ".cal";
        file_handler.open(fpath.str(), std::ios::in);
        if (file_handler.is_open()) {
            std::string line;
            std::getline(file_handler, line);
            std::istringstream iss(line);
            float main_vca_cal_1 = 0;
            float main_vca_cal_2 = 0;
            float offset = 0;
            float offset_od = 0;
            if (!(iss >> main_vca_cal_1 >> main_vca_cal_2 >> offset >> offset_od)) {
                printf("\nfailed voice cal load %d ", _voice_n);
            }
            _model_data.vca_l_offset = main_vca_cal_1;
            _model_data.vca_r_offset = main_vca_cal_2;
            _model_data.voice_offset = offset;
            _model_data.voice_offset_od = offset_od;
            // printf("\n loaded cal for voice %d %f %f %f", _voice_n, main_vca_cal_1, main_vca_cal_2, offset);
            return true;
        } else {
            return false;
        }
        return false;
    }

    void reloadCal() {
        bool res = loadFilterCal();
        bool res_2 = loadCalibration();
        printf("\nload cal %d %d", res, res_2);
    }

    bool loadFilterCal() {
        std::fstream file_handler;
        std::stringstream fpath;
        fpath << cal_path << "voice_" << _voice_n << "_filter.model";
        file_handler.open(fpath.str(), std::ios::in);
        if (file_handler.is_open()) {
            std::string line;
            std::getline(file_handler, line);
            std::istringstream iss(line);
            float a, b, c, d, e, f, g, h, temp;
            if (!(iss >> a >> b >> c >> d >> e >> f >> g >> h >> temp)) {
                printf("\nfailed filter voice cal load %d ", _voice_n);
                return false;
            }
            _model_data.cut_a = a;
            _model_data.cut_b = b;
            _model_data.cut_c = c;
            _model_data.cut_d = d;
            _model_data.cut_e = e;
            _model_data.cut_f = f;
            _model_data.cut_g = g;
            _model_data.cut_h = h;
            _model_data.cut_low_clip = -.99;
            // printf("\n %f %f %f %f %f %f ", a, b, c, d, e, f, temp);
            printf("\n loaded filter cal for voice %d", _voice_n);
            return true;
        }
        printf("\nfailed voice filter cal load %d ", _voice_n);

        // reset filter model to init state
        MoniqueAnalogModelData init_model;
        _model_data.cut_a = init_model.cut_a;
        _model_data.cut_b = init_model.cut_b;
        _model_data.cut_c = init_model.cut_c;
        _model_data.cut_d = init_model.cut_d;
        _model_data.cut_e = init_model.cut_e;
        _model_data.cut_f = init_model.cut_f;
        _model_data.cut_g = init_model.cut_g;
        _model_data.cut_h = init_model.cut_h;
        _model_data.cut_low_clip = -.99;
        return false;
    }

    bool saveCalibration() {
        std::fstream file_handler;
        std::stringstream fpath;
        fpath << cal_path << "voice_" << _voice_n << ".cal";
        file_handler.open(fpath.str(), std::ios::out | std::ios::trunc);
        if (file_handler.is_open()) {
            file_handler << _model_data.vca_l_offset << " "
                         << _model_data.vca_r_offset << " "
                         << _model_data.voice_offset << " "
                         << _model_data.voice_offset_od << " ";
            file_handler.close();
            // printf("\n saved voice %d cal to file %f %f %f", _voice_n, _model_data.vca_l_offset, _model_data.vca_r_offset, _model_data.voice_offset);
            return true;
        }
        return false;
    }

    void setMainLVcaOffset(float offset) {
        _model_data.vca_l_offset = offset;
    }

    void setMainRVcaOffset(float offset) {
        _model_data.vca_r_offset = offset;
    }

    void setVoiceOffset(float offset) {

        _model_data.voice_offset = offset;
    }

    void setVoiceOffsetOd(float offset) {

        _model_data.voice_offset_od = offset;
    }

    bool &getOdState() {
        return _voice_data.overdrive;
    }

    float &getOffsetOd() {
        return _model_data.voice_offset_od;
    }

    float &getOffset() {
        return _model_data.voice_offset;
    }

    MoniqueAnalogModelData &getModelData() {
        return _model_data;
    }

    std::array<float, 2> getFilterOut() {
        return {_last_cut_value, _last_res_value};
    }

    MoniqueAnalogModel(int voice_n, VoiceData &voice_data) :
        _voice_n(voice_n),
        _voice_data(voice_data) {
        printf("\n init analog voice %d", _voice_n);
        loadCalibration();
        bool res = loadFilterCal();
        _voice_data.filt_1_cut.fill(1);
        _voice_data.filt_2_cut.fill(1);
    };

    void setMute(bool mute) {
        _force_mute = mute;
    }

    void fxMuteLoopL(bool mute) {
        _loop_mute_l = mute;
    }

    void fxMuteLoopR(bool mute) {
        _loop_mute_r = mute;
    }

    void run(std::array<float, BUFFER_SIZE> &cv_buffer) {
        float max_l = 0;
        float max_r = 0;
        bool en_mute = (_mutes_allowed && !(_voice_data.last_allocated)) && (_voice_mute_en);
        en_mute = !((_voice_data.last_allocated));
        for (int cv_i = 0; cv_i < CV_BUFFER_SIZE; cv_i++) {
            const float &vca_l = std::abs(_voice_data.vca_l[cv_i]);

            const float &vca_r = std::abs(_voice_data.vca_r[cv_i]);

            max_l = max_l < vca_l ? vca_l : max_l;
            max_r = max_r < vca_l ? vca_l : max_r;
        }
        if (max_l < MUTE_THRESHOLD) {
            //_voice_data.mute_1 = en_mute;
        } else {
        }
        if (max_r < MUTE_THRESHOLD) {
            //_voice_data.mute_2 = en_mute;
        } else {
        }
        _voice_data.mute_1 |= _force_mute;
        _voice_data.mute_2 |= _force_mute;
        // generate the voice bit array which is const for each buffer
        bool filter_12db = !(_voice_data.filter_2_slope == FilterSlope::MOOG_2_POLE);
        int bit_array = ((int)(_voice_data.mute_1) << VoiceBitMap::VOICE_MUTE_L) +
                        ((int)(_voice_data.mute_2) << VoiceBitMap::VOICE_MUTE_R) +
                        ((int)!_voice_data.overdrive << VoiceBitMap::DRIVE_EN_N) +
                        ((int)filter_12db << VoiceBitMap::FILTER_TYPE) +
                        ((int)_loop_mute_l << VoiceBitMap::MIX_MUTE_L) +
                        ((int)_loop_mute_r << VoiceBitMap::MIX_MUTE_R);
        const auto float_bit_array = (float)((((double)bit_array) + 0.1) / (double)((2 << 22) - 1));
        int sample = 0;
        const auto &cal = _model_data;
        // if (_voice_n == 0)
        //     printf(" \n ");
        for (int cv_i = 0; cv_i < CV_BUFFER_SIZE; cv_i++) {
            {
                float tmp_c;
                // Run filter calibration

                const float &item = _voice_data.filt_2_cut[cv_i];
                tmp_c = item < 1.f ? item : 1.f;
                tmp_c = tmp_c > -0 ? tmp_c : -0;

                tmp_c = filter_model_function(cal, tmp_c, filter_12db);
                cv_buffer[sample + Cv1Order::Cv1FilterCut] = tmp_c;
                _last_cut_value = tmp_c;
                if (_voice_n == 0) {
                    // printf("\n filter %f %f ", item, tmp_c);
                }
            }

            {
                const float &item = unipolarClip(_voice_data.filt_2_res[cv_i]);
                const float res_cal = 2 * (item - item * item / 2);
                float tmp_res = (res_cal)*cal.res_gain + cal.res_zero_offset;
                if (_voice_n == 0) {
                    //       printf("\n %f ", tmp_res);
                }
                tmp_res = tmp_res < cal.res_high_clip ? tmp_res : cal.res_high_clip;
                tmp_res = tmp_res > cal.res_low_clip ? tmp_res : cal.res_low_clip;
                cv_buffer[sample + Cv1Order::Cv1FilterQ] = tmp_res;
                _last_res_value = tmp_res;
            }

            // run L & R vca calibration
            const float &vca_l = _voice_data.vca_l[cv_i];
            cv_buffer[sample + Cv1Order::Cv1AmpL] = vca_l + cal.vca_l_offset;
            const float &vca_r = _voice_data.vca_r[cv_i];
            cv_buffer[sample + Cv1Order::Cv1AmpR] = vca_r + cal.vca_r_offset;
            max_l = max_l < vca_l ? vca_l : max_l;
            max_r = max_r < vca_l ? vca_l : max_r;
            // if (_voice_n == 0)
            //     printf(" %1.1f", _voice_data.vca_l[cv_i]);

            // set the bit array sample
            cv_buffer[sample + Cv1BitArray] = float_bit_array;

            sample += CV_MUX_INC;
        }
    }

    MoniqueAnalogModelData _model_data;
    const int _voice_n;
    bool _voice_mute_en = true;
    bool _force_mute = true;
    bool _loop_mute_l = false;
    bool _loop_mute_r = false;
    bool _mutes_allowed = true;
    float _last_cut_value, _last_res_value;
    VoiceData &_voice_data;
};
} // namespace Monique