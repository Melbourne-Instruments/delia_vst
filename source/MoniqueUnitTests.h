
#include "MoniqueCommon.h"
#include "MoniqueLayerManager.h"
#include "common.h"
#include <cstdio>
#include <iostream>
#include <string_view>
#include <vector>

void saveVectorToWave(std::vector<float> samples, std::string name) {
#ifdef MONIQUE_NATIVE
    std::string file_base_dir = "./";
#else
    std::string file_base_dir = "/udata/";
#endif

    AudioFile<float> audio_data;
    std::fstream file_handler;
    audio_data.shouldLogErrorsToConsole(true);
    audio_data.setNumChannels(1);
    audio_data.setNumSamplesPerChannel(samples.size());
    int size = samples.size();

    printf("\n wt size %d\n", size);
    for (int i = 0; i < size; i++) {
        audio_data.samples[0][i] = (float)samples.at(i) / 4.f;
    }
    audio_data.setSampleRate(96000);
    audio_data.setBitDepth(24);
    audio_data.save(file_base_dir + name + ".wav");
}

void saveVectorToBinary(std::vector<float> samples, std::string name) {

#ifdef MONIQUE_NATIVE
    std::string file_base_dir = "./";
#else
    std::string file_base_dir = "/udata/";
#endif
    std::fstream file;
    std::ios_base::iostate exceptionMask = file.exceptions() | std::ios::failbit;
    file.exceptions(exceptionMask);
    std::string path = file_base_dir + name + ".bin";
    std::cout << path << std::endl;
    printf("\n samples %d ", (int)samples.size());
    file.open(path.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    file.write(reinterpret_cast<char *>(samples.data()), (samples.size()) * sizeof(float));
    file.close();
}

void loadAndRunLayer() {
    using namespace Monique;
    std::vector<int> data;
    int num_1 = std::rand();
    int num_2 = std::rand();

    // comment
    printf("\nhello world\n");
    int size = 10;
    for (int i = 0; i < size; i++) {
        data.push_back(0);
    }
    int i = 0;
    while (true) {
        num_1 = num_1 ^ num_2;
        data[i++] = num_1;
        if (i > size) {
            i = 0;
            break;
        }
    }
    std::array<float, Monique::BUFFER_SIZE> test_buffer;
    ProcessAudioData _audio_data;
    _audio_data.input_1 = &test_buffer;
    _audio_data.input_2 = &test_buffer;
    _audio_data.input_l = &test_buffer;
    _audio_data.input_r = &test_buffer;
    _audio_data.output_l = &test_buffer;
    _audio_data.output_r = &test_buffer;
    _audio_data.voice_cv.at(0) = &test_buffer;
    _audio_data.voice_cv.at(1) = &test_buffer;
    _audio_data.voice_cv.at(2) = &test_buffer;
    _audio_data.voice_cv.at(3) = &test_buffer;
    _audio_data.voice_cv.at(4) = &test_buffer;
    _audio_data.voice_cv.at(5) = &test_buffer;

    _audio_data.voice_high_res.at(0) = &test_buffer;
    _audio_data.voice_high_res.at(1) = &test_buffer;
    _audio_data.voice_high_res.at(2) = &test_buffer;
    _audio_data.voice_high_res.at(3) = &test_buffer;
    _audio_data.voice_high_res.at(4) = &test_buffer;
    _audio_data.voice_high_res.at(5) = &test_buffer;

    LayerManager layer_manager;
    layer_manager.paramRefresh();
    layer_manager.runSynth(_audio_data);
}

void parameterSmoothing() {
    using namespace Monique;
    std::vector<float> in_data;
    in_data.push_back(0);
    in_data.push_back(0);
    float level = 0;
    for (int i = 0; i < 100; i++) {
        if (i % 10 == 0) {
            level += .1;
        }
        in_data.push_back(1.f);
    }
    float filter_state = 0;
    for (int i = 0; i < (int)in_data.size() - 1; i++) {
        param_smooth(in_data.at(i), filter_state);
        in_data.at(i) = filter_state;
    }
    saveVectorToBinary(in_data, "param_smooth_test");
}

void runAllUnitTests() {
    loadAndRunLayer();
    parameterSmoothing();
}