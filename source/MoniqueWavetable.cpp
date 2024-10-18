/**
 * @file WavetableOsc.cpp
 * @brief
 *
 * @copyright Copyright (c) 2022-2024 Melbourne Instruments, Australia
 */
#include "MoniqueWavetable.h"

// Local static functions
static void *_monitor_wt(void *data);

namespace Monique {

// Constants
constexpr uint MIN_FREQ_HZ = 10;
constexpr uint MAX_FREQ_HZ = 20000;
constexpr int WT_MONITOR_RT_THREAD_PRIORITY = 1;
constexpr int WT_MONITOR_RT_THREAD_CPU_AFINITY = 2;

/**
 *-----------------------------------------------------------------------------
 * WavetableOsc class
 *-----------------------------------------------------------------------------
 */

void WavetableOsc::reCalculate() {
    // If a wavetable has been loaded
    _current_wavetable_a = _wavetable_loader_a.getCurrentWavetable();
    _current_wavetable_b = _wavetable_loader_b.getCurrentWavetable();
    const auto &wavetable_a = _current_wavetable_a->_samples;
    const auto &wavetable_b = _current_wavetable_b->_samples;
    // printf("\n here %d %d", _current_wavetable_a, _current_wavetable_b);
    if (_current_wavetable_a && _current_wavetable_b) {
        const bool interpolate = _interpolate;
        _wt_loaded = true;
        _num_waves_a = _current_wavetable_a->getNumWaves();
        _num_waves_b = _current_wavetable_b->getNumWaves();
        _output_position = 0;

        // precalc gain a  & b for factors that are updated at buffer rate
        _gain_a = fastCos((M_PIf32 / 2) * (_morph)) * 1;
        _gain_b = fastSin((M_PIf32 / 2) * (_morph)) * 1;
        for (uint cv_i = 0; cv_i < CV_BUFFER_SIZE; ++cv_i) {
            const float pitch = _pitch_cv_buffer.at(cv_i) + 1;
            float position = _shape_cv_buffer.at(cv_i);
            float scaled_pitch = pitch * NOTE_GAIN + ((float)((int)_slow_mode) * SLOW_WAVE_SUB);
            const float unison_pitch = scaled_pitch + _unison_osc_pitch * NOTE_GAIN;

            _phase_advance = (std::exp2f(scaled_pitch)) / (float)SAMPLE_RATE;
            const float phase_adv_unison = (std::exp2f(unison_pitch)) / (float)SAMPLE_RATE;

            float unison_osc_1 = _unison_envelope_1.at(cv_i);
            float unison_osc_2 = _unison_envelope_2.at(cv_i);

            // scale and clip the position signal
            position = std::min(1.0f, std::max(0.f, position));

            // calculate final gain based on current wt volume
            const float gain_a = _gain_a;
            const float gain_b = _gain_b;
            const float wt_out_vol = _wt_vol;

            // filter the position signal for sound quality and then calculate the wave number based on the size of each wavetable
            _position_filter_state += WT_POS_SMOOTH_COEFF * (position - _position_filter_state);
            int mipmap_index;
            int mipmap_index_uni;
            mipmap_index = _mipmap_index;
            mipmap_index_uni = _mipmap_index_uni;
            if (interpolate) {
                mipmap_index = (int)std::floor(scaled_pitch - 4.4) - 2;

                mipmap_index = std::min(7, std::max(0, mipmap_index));

                int mipmap_index_uni = (int)std::floor(unison_pitch - 4.4) - 2;
                mipmap_index_uni = std::min(7, std::max(0, mipmap_index_uni));
            }
            uint mipmap_offset = getMipmapOffset(mipmap_index);
            // calculate mipmap vars
            uint mipmap_offset_a = mipmap_offset * _num_waves_a;
            uint mipmap_offset_b = mipmap_offset * _num_waves_b;
            uint phase_mult_mip_level = std::min(mipmap_index, 4);
            uint phase_mult = (WAVE_LENGTH >> phase_mult_mip_level);
            uint phase_mask = phase_mult - 1;
            _mipmap_index = mipmap_index;

            // unison variables
            uint mipmap_offset_uni = getMipmapOffset(mipmap_index_uni);
            uint mipmap_offset_a_uni = mipmap_offset_uni * _num_waves_a;
            uint mipmap_offset_b_uni = mipmap_offset_uni * _num_waves_b;
            uint phase_mult_mip_level_uni = std::min(mipmap_index_uni, 4);
            uint phase_mult_uni = (WAVE_LENGTH >> phase_mult_mip_level_uni);
            uint phase_mask_uni = phase_mult_uni - 1;

            if (interpolate) {

                // calculate floor and ceil of position, this is used to fetch the samples. we then perform linear interp on the samples.
                const uint pos_a = std::floor(_position_filter_state * ((float)_num_waves_a - 1));
                const uint pos_b = std::floor(_position_filter_state * ((float)_num_waves_b - 1));
                const uint pos_a_low = pos_a * phase_mult;
                const uint pos_b_low = pos_b * phase_mult;
                const uint pos_a_high = (pos_a + 1) * phase_mult;
                const uint pos_b_high = (pos_b + 1) * phase_mult;
                const float morph_frac_a = _position_filter_state * ((float)_num_waves_a - 1) - (float)pos_a;
                const float morph_frac_b = _position_filter_state * ((float)_num_waves_b - 1) - (float)pos_b;

                // unison vars
                const uint pos_a_low_uni = pos_a * phase_mult_uni;
                const uint pos_b_low_uni = pos_b * phase_mult_uni;
                const uint pos_a_high_uni = (pos_a + 1) * phase_mult_uni;
                const uint pos_b_high_uni = (pos_b + 1) * phase_mult_uni;

                for (int i = 0; i < WT_BUFFER_SIZE; ++i) {
                    // increment and wrap the current phase.
                    _current_phase += _phase_advance;
                    _current_phase = std::fmod(_current_phase, 1.f);

                    // unison
                    _current_phase_uni += phase_adv_unison;
                    _current_phase_uni = std::fmod(_current_phase_uni, 1.f);

                    // calculate the position of sample 1
                    float pos = _current_phase * (float)phase_mult;
                    const uint sample_pos_1 = ((uint)std::floor(pos)) & phase_mask;

                    // wrap sample 2 position
                    const uint sample_pos_2 = (sample_pos_1 + 1) & phase_mask;
                    float frac_part = pos - (float)sample_pos_1;

                    // unison
                    float pos_uni = _current_phase_uni * (float)phase_mult_uni;
                    const uint sample_pos_1_uni = ((uint)std::floor(pos_uni)) & phase_mask_uni;

                    // wrap sample 2 position
                    const uint sample_pos_2_uni = (sample_pos_1_uni + 1) & phase_mask_uni;
                    float frac_part_uni = pos_uni - (float)sample_pos_1_uni;

                    _wt_out_a = getInterpSample(morph_frac_a, wavetable_a, mipmap_offset_a, pos_a_low, pos_a_high, sample_pos_1, sample_pos_2, frac_part);
                    _wt_out_b = getInterpSample(morph_frac_b, wavetable_b, mipmap_offset_b, pos_b_low, pos_b_high, sample_pos_1, sample_pos_2, frac_part);

                    // unison
                    float wt_a_uni = getInterpSample(morph_frac_a, wavetable_a, mipmap_offset_a_uni, pos_a_low_uni, pos_a_high_uni, sample_pos_1_uni, sample_pos_2_uni, frac_part_uni);
                    float wt_b_uni = getInterpSample(morph_frac_b, wavetable_b, mipmap_offset_b_uni, pos_b_low_uni, pos_b_high_uni, sample_pos_1_uni, sample_pos_2_uni, frac_part_uni);

                    // mix the A and B wavetable outputs
                    float morphed_sum = gain_a * _wt_out_a + gain_b * _wt_out_b;
                    float morphed_sum_uni = (gain_a * wt_a_uni + gain_b * wt_b_uni) * 2;

                    _wt_cv_out = morphed_sum * 2.f;

                    _audio_buffer[cv_i * WT_BUFFER_SIZE + i] = wt_out_vol * (_wt_cv_out * unison_osc_1 + unison_osc_2 * morphed_sum_uni);
                }
            } else {
                uint wave_pos_offset_a = _wave_pos_a * phase_mult;
                uint wave_pos_offset_b = _wave_pos_b * phase_mult;

                // unison
                uint wave_pos_offset_a_uni = _wave_pos_a * phase_mult_uni;
                uint wave_pos_offset_b_uni = _wave_pos_b * phase_mult_uni;
                for (int i = 0; i < WT_BUFFER_SIZE; ++i) {

                    // increment the phase and wrap. if the value wraps, we also update the wavetable number, so wavetables only change when the phase wraps and there is a zero crossing.
                    _current_phase += _phase_advance;
                    if (_current_phase >= 1.f) {
                        _current_phase = std::fmod(_current_phase, 1.f);

                        mipmap_index = (int)std::floor(scaled_pitch - 4.4) - 2;

                        // update the mipmap vars on a zerocrossing
                        mipmap_index = std::min(7, std::max(0, mipmap_index));
                        mipmap_offset = getMipmapOffset(mipmap_index);
                        mipmap_offset_a = mipmap_offset * _num_waves_a;
                        mipmap_offset_b = mipmap_offset * _num_waves_b;
                        phase_mult_mip_level = std::min(mipmap_index, 4);
                        phase_mult = (WAVE_LENGTH >> phase_mult_mip_level);
                        phase_mask = phase_mult - 1;
                        _mipmap_index = mipmap_index;

                        // update the wave position if the phase has wrapped to avoid clicks
                        _wave_pos_a = ((uint)std::round(position * ((float)_num_waves_a - 1)));
                        _wave_pos_b = ((uint)std::round(position * ((float)_num_waves_b - 1)));
                        wave_pos_offset_a = _wave_pos_a * phase_mult;
                        wave_pos_offset_b = _wave_pos_b * phase_mult;
                    }

                    // unison

                    _current_phase_uni += phase_adv_unison;
                    if (_current_phase_uni >= 1.f) {
                        _current_phase_uni = std::fmod(_current_phase_uni, 1.f);

                        mipmap_index_uni = (int)std::floor(unison_pitch - 4.4) - 2;

                        // update the mipmap vars on a zerocrossing
                        mipmap_index_uni = std::min(7, std::max(0, mipmap_index_uni));
                        mipmap_offset_uni = getMipmapOffset(mipmap_index_uni);
                        mipmap_offset_a_uni = mipmap_offset_uni * _num_waves_a;
                        mipmap_offset_b_uni = mipmap_offset_uni * _num_waves_b;
                        phase_mult_mip_level_uni = std::min(mipmap_index_uni, 4);
                        phase_mult_uni = (WAVE_LENGTH >> phase_mult_mip_level_uni);
                        phase_mask_uni = phase_mult_uni - 1;

                        // update the wave position if the phase has wrapped to avoid clicks
                        wave_pos_offset_a_uni = _wave_pos_a * phase_mult_uni;
                        wave_pos_offset_b_uni = _wave_pos_b * phase_mult_uni;
                        _mipmap_index_uni = mipmap_index_uni;
                    }

                    // get the position of sample 1
                    float pos = _current_phase * (float)(phase_mult);
                    const uint sample_pos_1 = (uint)std::floor(pos) & phase_mask;

                    // wrap sample 2 position
                    const uint sample_pos_2 = (sample_pos_1 + 1) & (phase_mask);
                    float frac_part = pos - (float)sample_pos_1;

                    // Get the two samples for each wt
                    _wt_out_a = getSample(wavetable_a, frac_part, mipmap_offset_a, wave_pos_offset_a, sample_pos_1, sample_pos_2);
                    _wt_out_b = getSample(wavetable_b, frac_part, mipmap_offset_b, wave_pos_offset_b, sample_pos_1, sample_pos_2);
                    // printf("\n %f %f", _wt_out_b, gain_b);

                    // unison
                    //  get the position of sample 1
                    float pos_uni = _current_phase_uni * (float)(phase_mult_uni);
                    const uint sample_pos_1_uni = (uint)std::floor(pos_uni) & phase_mask_uni;

                    // wrap sample 2 position
                    const uint sample_pos_2_uni = (sample_pos_1_uni + 1) & (phase_mask_uni);
                    float frac_part_uni = pos_uni - (float)sample_pos_1_uni;

                    // Get the two samples for each wt
                    float a_uni = getSample(wavetable_a, frac_part_uni, mipmap_offset_a_uni, wave_pos_offset_a_uni, sample_pos_1_uni, sample_pos_2_uni);
                    if (dump) {
                        // printf("\n %d %d %d %f  %d %d %d %f", mipmap_offset_a, wave_pos_offset_a, sample_pos_1, _current_phase, mipmap_offset_a_uni, wave_pos_offset_a_uni, sample_pos_1_uni, _current_phase_uni);
                    }
                    float b_uni = getSample(wavetable_b, frac_part_uni, mipmap_offset_b_uni, wave_pos_offset_b_uni, sample_pos_1_uni, sample_pos_2_uni);

                    // mix the A and B wavetable outputs
                    float morphed_sum = gain_a * _wt_out_a + gain_b * _wt_out_b;
                    _wt_cv_out = morphed_sum * 2.f;

                    // unison output
                    float morphed_sum_uni = (gain_a * a_uni + gain_b * b_uni) * 2;

                    // write the output to the array
                    _audio_buffer[cv_i * WT_BUFFER_SIZE + i] = wt_out_vol * (_wt_cv_out * unison_osc_1 + unison_osc_2 * morphed_sum_uni);
                }
            }
        }
    }

    else {
        _wt_loaded = false;
        _audio_buffer.fill(0);
    }

    _buffer_position = 0;
}

float WavetableOsc::run() {

    const uint &bp = _buffer_position;
    _pitch_cv_buffer[bp] = _pitch;
    _shape_cv_buffer[bp] = _position;

    for (int i = 0; i < WT_BUFFER_SIZE; ++i) {
        uint bp_audio = WT_BUFFER_SIZE * bp + i;
        (*_output).at(bp_audio) = _audio_buffer[bp_audio];
    }

    _buffer_position++;
    return _audio_buffer[BUFFER_SIZE - 1];
}

void WavetableOsc::reset() {
    // Reset oscillator data
    _current_value = 0;
    _current_phase = 0;
    _phase_advance = 0;
    _mipmap_index = 0;
}

/**
 *-----------------------------------------------------------------------------
 * WavetableLoader class
 *-----------------------------------------------------------------------------
 */
static int wt_loader_ctr = 0;

WavetableLoader::WavetableLoader() {
    // Kick-off a worker thread to monitor wavetable changes
    _current_wavetable = 0;
    _exit_thread = false;
    _load_wavetable = false;
    _wt_select_num = -1.0;
    _wt_ctr = wt_loader_ctr++;

    // Create a RT thread (not Xenomai) for monitoring wavetable changes
    // Initialise the RT thread attributes
    pthread_attr_t task_attributes;
    pthread_attr_init(&task_attributes);
    pthread_attr_setdetachstate(&task_attributes, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setinheritsched(&task_attributes, PTHREAD_EXPLICIT_SCHED);
    pthread_attr_setschedpolicy(&task_attributes, SCHED_FIFO);

    // Set the RT thread affinity to CPU core 2
    cpu_set_t cpus;
    CPU_ZERO(&cpus);
    CPU_SET(WT_MONITOR_RT_THREAD_CPU_AFINITY, &cpus);
    pthread_attr_setaffinity_np(&task_attributes, sizeof(cpu_set_t), &cpus);

    // Create and spawn the RT thread
    struct sched_param rt_params = {.sched_priority = WT_MONITOR_RT_THREAD_PRIORITY};
    pthread_attr_setschedparam(&task_attributes, &rt_params);
    pthread_create(&_load_wavetable_thread, &task_attributes, _monitor_wt, this);
}

WavetableLoader::~WavetableLoader() {
    // Wait for any worker threads to finish
    _exit_thread = true;
    pthread_join(_load_wavetable_thread, nullptr);
}

const Wavetable *WavetableLoader::getCurrentWavetable() const {
    // Return the current wavetable
    return _current_wavetable.load();
}

void WavetableLoader::loadWavetable(float select) {
    if (select != _wt_select_num) {
        // Save the wavetable select number to load, and let the monitor
        // wavetable thread know a new wavetable needs loading
        _wt_select_num = select;
        _load_wavetable = true;
    }
}

void WavetableLoader::monitorWavetable() {
    // Do forever or until exited
    while (!_exit_thread) {
        // Is there a wavetable ready to load?
        if (_load_wavetable) {
            try {
#ifdef _NINA_UNIT_TESTS
                auto start = std::chrono::steady_clock::now();
#endif
                std::vector<std::string> filenames;
                AudioFile<float> file;
                struct dirent **dirent = nullptr;
                int num_files;

                // Scan the Sushi wavetable folder
                num_files = ::scandir(common::MONIQUE_WT_DIR, &dirent, 0, ::versionsort);
                if (num_files > 0) {
                    // Process each file in the folder
                    for (uint i = 0; i < num_files; i++) {
                        // If we've not found the max number of wavetables yet and this a normal file
                        if ((filenames.size() < common::MAX_NUM_WAVETABLE_FILES) && (dirent[i]->d_type == DT_REG)) {
                            // If it has a WAV file extension
                            auto name = std::string(dirent[i]->d_name);
                            if (name.substr((name.size() - (sizeof(".wav") - 1))) == ".wav") {
                                // Add the filename
                                filenames.push_back(dirent[i]->d_name);
                            }
                        }
                        ::free(dirent[i]);
                    }
                }
                if (dirent) {
                    ::free(dirent);
                }

                // Are there any wavetables to process?
                if (filenames.size() > 0) {
                    // Load the wavetable
                    auto current_wt_select_num = _wt_select_num;
                    auto filename = filenames.at((uint)std::round(((float)(filenames.size()) * _wt_select_num)));
                    if (!file.load(MONIQUE_WT_FILE(filename)))
                        throw std::invalid_argument("Wavetable does not exist");

                    // Check the number of samples is valid
                    if (file.getNumChannels() == 0 ||
                        (file.samples[0].size() % WAVE_LENGTH))
                        throw std::runtime_error(
                            "Wavetable number of channels/samples is invalid");

                    // Get the number of waves and check it is valid
                    auto num_waves = file.samples[0].size() / WAVE_LENGTH;
                    if (num_waves > MAX_NUM_WAVES)
                        throw std::runtime_error("Wavetable number of waves is invalid");

                    // Get a pointer to the free wavetable slot
                    Wavetable *wavetable = &_wavetable_slot1;
                    if (_current_wavetable == &_wavetable_slot1)
                        wavetable = &_wavetable_slot2;

                    // Reset the wavetable
                    wavetable->reset();

                    // Process the wavetable samples
                    const float *samples = file.samples[0].data();
                    wavetable->processWaves(num_waves, samples);

                    // Has the selected wavetable changed during the load?
                    // If so, loop again to load the new wavetable
                    // printf("\nloaded wt: %s\n", filename.c_str());
                    if (current_wt_select_num == _wt_select_num) {
                        // The wavetable has been loaded
                        _load_wavetable = false;
                        _current_wavetable.store(wavetable);
                    }
                    printf("\n %d loaded %s \n", _wt_ctr, filename.c_str());
                    fflush(stdout);
#ifdef _NINA_UNIT_TESTS
                    auto finish = std::chrono::steady_clock::now();
                    float tt = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
                    printf("\nLoad wavetable: time taken = %d ms", (int)tt);
                    fflush(stdout);
#endif
                }
            } catch (...) {
                // Catch all if any error happens during processing
                // Just ignore any exceptions for now, wavetable is
                // not loaded
                _load_wavetable = false;
            }
        }

        // Sleep for 50ms
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
}

/**
 *-----------------------------------------------------------------------------
 * Wavetable class
 *-----------------------------------------------------------------------------
 */

void Wavetable::reset() {
    // Reset the filters
    _num_waves = 0;
}

uint Wavetable::getNumWaves() const {
    // Return the number of waves
    return _num_waves;
}

void Wavetable::processWaves(uint num_waves, const float *samples) {
    // Process the wave samples
    std::array<WavetableLowpassFilter, NUM_LOWPASS_IN_FILTERS + 1> lowpass_in =
        {
            // first filter does nothing
            WavetableLowpassFilter(20001u),
            WavetableLowpassFilter(10000u),
            WavetableLowpassFilter(5000u),
            WavetableLowpassFilter(2500u),
            WavetableLowpassFilter(1250u),
            WavetableLowpassFilter(625u),
            WavetableLowpassFilter(313u),
            WavetableLowpassFilter(157u),
        };

    std::array<float, WAVE_LENGTH * 3> wave_process_buffer;
    wave_process_buffer.fill(0);
    int dec_sample = 0;
    _num_waves = num_waves;

    for (int mipmap = 0; mipmap < NUM_MIPMAPS; mipmap++) {
        const float *smp_ptr = samples;
        int mipmap_lim = mipmap;
        if (mipmap_lim > 4)
            mipmap_lim = 4;
        int phase_inc = 1 << (mipmap_lim);
        for (int wave_num = 0; wave_num < num_waves; wave_num++) {

            // Process and filter the samples for the current wave

            // run the filter in the forward direction
            for (int i = 0; i < WAVE_LENGTH * 3; i++) {
                float sample = smp_ptr[i & 0x7FF];

                wave_process_buffer.at(i) = lowpass_in.at(mipmap).process(sample);
            }

            // run the filter backwards to remove the phase offset and further filter
            for (int i = WAVE_LENGTH * 3 - 1; i >= 0; i--) {
                wave_process_buffer.at(i) = lowpass_in.at(mipmap).process(wave_process_buffer.at(i));
            }

            // decimate and copy samples to the sample buffer
            for (int sample_n = 0 + (phase_inc / 2); sample_n < WAVE_LENGTH; sample_n += phase_inc) {
                _samples.at(dec_sample++) = _floatSampleToInt16(wave_process_buffer.at(WAVE_LENGTH + sample_n));
            }
            // move the sample pointer to the next wave
            smp_ptr += WAVE_LENGTH;
        }
    }
}

/**
 *-----------------------------------------------------------------------------
 * WavetableWave class
 *-----------------------------------------------------------------------------
 */

inline float WavetableWave::getSample(uint sample_num, uint mipmap_index) const {
    // Get the sample from the specified wave
    return _int16SampleToFloat(_wavetable_wave[((sample_num << 3) + mipmap_index)]);
}

inline const int16_t *WavetableWave::getSampleAddr(uint sample_num, uint mipmap_index) const {
    return &(_wavetable_wave[0]);
}

/**
 *-----------------------------------------------------------------------------
 * WavetableLowpassFilter class
 *-----------------------------------------------------------------------------
 */

inline void WavetableLowpassFilter::reset() {
    // Reset the biquad filters
    _biquad_1.reset();
    _biquad_2.reset();
    _biquad_3.reset();
    _biquad_4.reset();
}

inline float WavetableLowpassFilter::process(float input) {
    // Process the biquad filters in sequence

    // handle bypass processing
    if (bypass) {
        return input;
    }

    // handle normal processing
    return _biquad_4.process(_biquad_3.process(_biquad_2.process(_biquad_1.process(input))));
}

/**
 *-----------------------------------------------------------------------------
 * BiquadFilter class
 *-----------------------------------------------------------------------------
 */

BiquadFilter::BiquadFilter(float a1, float a2, float b0, float b1, float b2) :
    _a1(a1), _a2(a2), _b0(b0), _b1(b1), _b2(b2) {
    // Reset the filer
    reset();
}
} // namespace Monique

//----------------------------------------------------------------------------
// _monitor_wt
//----------------------------------------------------------------------------
static void *_monitor_wt(void *data) {
    // printf("\nHere1\n");
    auto wt_loader = static_cast<Monique::WavetableLoader *>(data);
    wt_loader->monitorWavetable();

    // To suppress warnings
    return nullptr;
}
