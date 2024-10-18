
#pragma once

#include "common.h"
#include <cstring>
#include <mqueue.h>

namespace Steinberg {
namespace Vst {
namespace MoniqueFactory {
constexpr int CHANNELS = 4;
constexpr int SAMPLER_MSG_SIZE = 1;

// GUI message
struct SamplerMsg {
    int samples;

    // Constructor
    SamplerMsg() {
        // Set the samples buffer to all zeros
        samples = 1;
    }
};

class SamplerMsgQueue {
  public:
    SamplerMsgQueue();
    ~SamplerMsgQueue();

    bool open();
    bool post_msg();
    void close();

  private:
    mqd_t _sampler_mq_desc;
};

} // namespace MoniqueFactory
} // namespace Vst
} // namespace Steinberg