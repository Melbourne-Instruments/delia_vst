/**
 * @file GuiMsgQueue.h
 * @brief
 * @version 0.1
 * @date 2021-10-27
 *
 * @copyright Copyright (c) 2023 Melbourne Instruments, Australia
 */

#pragma once

#include "MoniqueCommon.h"
#include "common.h"
#include <cstring>
#include <mqueue.h>

namespace Monique {
// Constants
constexpr char GUI_MSG_QUEUE_NAME[] = "/delia_samples_msg_queue";
constexpr uint GUI_MSG_QUEUE_SIZE = common::GUI_NUM_BUFFERS;

// GUI message
struct GuiMsg {
    float samples[common::GUI_SAMPLES_SIZE];

    // Constructor
    GuiMsg() {
        // Set the samples buffer to all zeros
        std::memset(samples, 0, sizeof(samples));
    }
};

class GuiMsgQueue {
  public:
    GuiMsgQueue() {
        // Initialise data
        _gui_mq_desc = (mqd_t)-1;
    }

    ~GuiMsgQueue() {
        close();
    }

    bool open() {
        mq_attr attr;

        // Open the GUI Message Queue
        std::memset(&attr, 0, sizeof(attr));
        attr.mq_maxmsg = GUI_MSG_QUEUE_SIZE;
        attr.mq_msgsize = sizeof(GuiMsg);
        _gui_mq_desc = ::mq_open(
            GUI_MSG_QUEUE_NAME, (O_WRONLY),
            (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH), &attr);
        return _gui_mq_desc != (mqd_t)-1;
    }

    bool post_msg(const GuiMsg &msg) {
        // If the GUI Message Queue is valid
        if (_gui_mq_desc != (mqd_t)-1) {
            // Write the message to the GUI Message Queue
            auto res = (::mq_send(_gui_mq_desc, (char *)&msg, sizeof(msg), 0));
            if (res == 0)
                return true;
        }
        return false;
    }

    void close() {
        // Is the GUI Message Queue open?
        if (_gui_mq_desc != (mqd_t)-1) {
            // Close the GUI FIFO
            ::mq_close(_gui_mq_desc);
            ::mq_unlink(GUI_MSG_QUEUE_NAME);
            _gui_mq_desc = (mqd_t)-1;
        }
    }

  private:
    mqd_t _gui_mq_desc;
};

} // namespace Monique