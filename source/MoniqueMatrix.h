/**
 * @file NinaMatrix.h
 * @brief Declaration of Nina's matrix
 * @date 2022-07-18
 *
 * Copyright (c) 2023 Melbourne Instruments
 *
 */
#pragma once

#include "MoniqueCommon.h"
#include "common.h"

//#define MATRIX_DBG

namespace Monique {

class MoniqueMatrix {
  public:
    int voicen = 0;
    bool dump = false;
    MoniqueMatrix(std::array<float, TOTAL_MOD_DSTS> &drift_data) :
        _drift_data(drift_data){};

    ~MoniqueMatrix(){};

    void add_source(int source_num, const float *const source, bool is_fast = false);
    void add_dst(int destination_num, float *const destination, bool is_fast = false);

    void set_gain_addr(int source_num, int dest_num, float *gain);

    /**
     * @brief calculate the sum of mod_gain * slot_dst_pairs for each dst and apply. also add the cached slow sum
     *
     */
    void run_fast_slots() {
        const float buffer_counter_mult = (float)_buffer_counter;
        for (uint dst = 0; dst < FAST_DESTS; ++dst) {
            float sum = 0;
            for (uint src = 0; src < FAST_SRCS; ++src) {
                sum += *(_data._fast_src[src]) * (*_data._fast_gains[(src * FAST_DESTS + dst)]);
#ifdef MATRIX_DBG
                assert(floatIsValid(*(_data._fast_src[src])));
#endif
            }

            // linear interp on slow sources as a dezip method
            float slow_cache_smooth = _dst_cache[dst] + _dst_cache_delta[dst] * buffer_counter_mult;
            *(_data._fast_dst[dst]) = sum + slow_cache_smooth;
        }
        _buffer_counter++;
    }

    /**
     * @brief calculate the sum of all slow gain * slot_dst_pairs for each dst and cache the value for the fast slot function to use
     *
     */
    void run_slow_slots() {

        int fast_dst_counter = 0;

        // TODO: test caching the sources here
        for (int i = 0; i < TOTAL_MOD_DSTS; ++i) {
#ifdef MATRIX_DBG
            bool printing = (i == (uint)ModMatrixDst::FX_SEND);
#endif
            bool fast_dst = false;
            if (_data._matrix_dests[i] == _data._fast_dst[fast_dst_counter]) {
                fast_dst = true;
            }
            float sum = 0.0f;
            for (int j = 0; j < TOTAL_MOD_SRC; ++j) {
                sum += ((*(_data._matrix_srcs[(j)])) * (*_data._gains[(j * TOTAL_MOD_DSTS + i)]));
#ifdef MATRIX_DBG
                if (!floatIsValid((*(_data._matrix_srcs[(j)])))) {
                    printf("\n %d %d %f %f\n", i, j, ((*(_data._matrix_srcs[(j)]))), ((*(_data._gains[(j * TOTAL_MOD_DSTS + i)]))));
                    printf("\n %d %d %f %f\n", i, j, ((*(_data._matrix_srcs[(j)]))), ((*(_data._gains[(j * TOTAL_MOD_DSTS + i)]))));
                }
                assert(floatIsValid((*(_data._matrix_srcs[(j)]))));
                assert(floatIsValid(sum));

                if (i == (int)ModMatrixDst::FX_SEND && printing) {
                    if (fabs(*_data._gains.at(j * (uint)ModMatrixDst::NUM_DSTS + i)) > 0.001) {
                        printf("\n");
                        printf(src_names[(j)].c_str());
                        printf("\t\t\t");
                        printf(dst_names[i].c_str());
                        printf("\t\t%f * %f  %f ", *_data._matrix_srcs.at(j), *_data._gains.at(j * (uint)ModMatrixDst::NUM_DSTS + i), sum);
                    }
                    if (j == (int)ModMatrixSrc::LFO_1 && printing) {
                        printf("\n filt lfo gain %f", *_data._fast_gains.at((int)ModMatrixSrc::LFO_1 * (uint)FAST_DESTS + (int)ModMatrixDst::AMP_EG_SUS));
                    }
                }
#endif
            }

#ifdef MATRIX_DBG
            if (i == (uint)ModMatrixDst::FX_SEND) {
                printf("\n dst %f", sum);
            }
#endif
            if (dump && fabs((_drift_data.at(i)) > 0.0001)) {
                // printf("\n dst %d %f", i, _drift_data.at(i));
            }
            if (fast_dst) {
                const float current_val = _dst_cache[fast_dst_counter];
                _dst_cache[fast_dst_counter] = sum + _drift_data.at(i);
                _dst_cache_delta[fast_dst_counter] = (sum - current_val) / ((float)BUFFER_SIZE);
                if (fast_dst_counter < (FAST_DESTS - 1)) {
                    fast_dst_counter++;
                }
            } else {

                *(_data._matrix_dests[i]) = sum + _drift_data.at(i);
            }
        }
#ifdef MATRIX_DBG
        float v = *_data._matrix_dests.at((uint)ModMatrixDst::FX_SEND);
#endif
        print = false;

        // reset the slow buffer counter which is used for slow source dezip
        _buffer_counter = 0;
    }

    bool print = false;

    ModMatrixData &getData() {
        return _data;
    }

  private:
    float zero = 0.0f;
    int _buffer_counter = 0;
    ModMatrixData _data;
    std::array<float, TOTAL_MOD_DSTS> &_drift_data;
    std::array<float, FAST_DESTS> _dst_cache = {0.0};
    std::array<float, FAST_DESTS> _dst_cache_delta = {0.0};
};

} // namespace Monique