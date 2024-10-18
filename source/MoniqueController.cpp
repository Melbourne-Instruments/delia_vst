/**
 * @file NinaController.cpp
 * @brief Nina Edit Controller implementation.
 *
 * @copyright Copyright (c) 2022-2024 Melbourne Instruments, Australia
 */
#include "MoniqueController.h"
#include "base/source/fstreamer.h"
#include "base/source/fstring.h"
#include "common.h"
#include "monique_synth_parameters.h"
#include "pluginterfaces/base/futils.h"
#include "pluginterfaces/base/ustring.h"
#include "pluginterfaces/vst/ivstmidicontrollers.h"
#include <sstream>
#include <string>

namespace Steinberg {
namespace Vst {
namespace MoniqueSynth {

/**
 * @note	Define GUID for controller
 *
 */
FUID Controller::uid(0xF02697F4, 0x2DA144E5, 0x90E72DDE, 0x0E69CC22);

tresult PLUGIN_API Controller::initialize(FUnknown *context) {

    // Initialise the base edit controller
    tresult result = EditController::initialize(context);
    if (result == kResultTrue) {
        Parameter *param;
        auto state_def = Monique::normalised_state_param_init();
        auto common_def = Monique::common_param_init();
        auto preset_common = Monique::preset_common_param_init();
        auto layer_def = Monique::layer_param_init();
        auto global_def = Monique::global_param_init();
        auto preset_common_def = Monique::preset_common_param_init();

        // Add global params
        Monique::ParamNameGenerator name_gen;
        for (uint i = 0; i < Monique::GlobalParams::NUM_GLOBAL_PARAMS; i++) {
            auto pid = Monique::gen_param_id((Monique::GlobalParams)i);
            auto name = name_gen.get_parameter_name(pid);
            param = new RangeParameter(USTRING(name.c_str()), pid, USTRING(""), 0., 1., global_def.at(i));
            parameters.addParameter(param);
        }

        // Add preset common params
        for (uint i = 0; i < Monique::NUM_PRESET_COMMON_PARAMS; i++) {
            auto pid = Monique::gen_param_id((Monique::PresetCommonParameters)i);
            auto name = name_gen.get_parameter_name(pid);
            param = new RangeParameter(USTRING(name.c_str()), pid, USTRING(""), 0., 1., preset_common_def.at(i));
            parameters.addParameter(param);
        }

        // generate params unique to each layer
        for (int layer_n = 0; layer_n < Monique::NUM_LAYERS; layer_n++)

        // generate layer params
        {
            for (uint i = 0; i < Monique::LayerParameters::NUM_LAYER_PARAMS; i++) {
                auto pid = Monique::gen_param_id((Monique::LayerParameters)i, layer_n);
                auto name = name_gen.get_parameter_name(pid);
                param = new RangeParameter(USTRING(name.c_str()), pid, USTRING(""), 0., 1., layer_def.at(i));
                parameters.addParameter(param);
            }

            // add layer common parameters
            for (uint i = 0; i < Monique::CommonParameters::NUM_COMMON_PARAMS; i++) {
                auto pid = Monique::gen_param_id((Monique::CommonParameters)i, layer_n);
                auto name = name_gen.get_parameter_name(pid);
                param = new RangeParameter(USTRING(name.c_str()), pid, USTRING(""), 0., 1., common_def.at(i));
                parameters.addParameter(param);
            }

            // add parameters for each state
            for (int state_n = 0; state_n < (int)Monique::State::NUM_STATES; state_n++) {
                Monique::State state = (Monique::State)state_n;
                for (uint i = 0; i < (int)Monique::StateParameters::NUM_STATE_PARAMS; i++) {
                    auto pid = Monique::gen_param_id((Monique::StateParameters)i, state, layer_n);
                    auto name = name_gen.get_parameter_name(pid);
                    param = new RangeParameter(USTRING(name.c_str()), pid, USTRING(""), 0., 1., state_def.at(i));
                    parameters.addParameter(param);
                }

                // add mod matrix parameters
                for (int src = 0; src < (int)Monique::ModMatrixSrc::NUM_SRCS; src++) {
                    for (int dst = 0; dst < (int)Monique::ModMatrixDst::NUM_DSTS; dst++) {
                        const int id = (int)Monique::get_state_id_from_mod_matrix((Monique::ModMatrixSrc)src, (Monique::ModMatrixDst)dst);
                        auto pid = Monique::gen_param_id((Monique::ModMatrixSrc)src, (Monique::ModMatrixDst)dst, state, layer_n);
                        auto name = name_gen.get_parameter_name(pid);
                        param = new RangeParameter(USTRING(name.c_str()), pid, USTRING(""), 0., 1., (state_def.at(id)));
                        parameters.addParameter(param);
                    }
                }
            }
        }
    }
    return kResultTrue;
}

tresult PLUGIN_API Controller::terminate() {
    return EditController::terminate();
}

//------------------------------------------------------------------------
ParamID PLUGIN_API Controller::firstLayerParamId() {
    return Monique::NUM_GLOBAL_PARAMS + Monique::NUM_PRESET_COMMON_PARAMS;
}

//------------------------------------------------------------------------
ParamID PLUGIN_API Controller::firstLayerStateParamId() {
    return Monique::FIRST_STATE_PARAM;
}

//------------------------------------------------------------------------
int32 PLUGIN_API Controller::getNumLayerParams() {
    return Monique::NUM_LAYER_PARAMS + Monique::NUM_COMMON_PARAMS + (2 * Monique::TOTAL_STATE_PARAMS);
}

//------------------------------------------------------------------------
int32 PLUGIN_API Controller::getNumLayerStateParams() {
    return Monique::TOTAL_STATE_PARAMS;
}

tresult PLUGIN_API Controller::setParamNormalizedFromFile(ParamID tag,
    ParamValue value) {

    Parameter *pParam = EditController::getParameterObject(tag);

    if (!pParam)
        return kResultFalse;

    return setParamNormalized(tag, pParam->toNormalized(value));
}

tresult PLUGIN_API Controller::setComponentState(IBStream *fileStream) {

    return kResultTrue;
}

/**
 * @note will be queried 129 times for control messages
 *
 */
tresult PLUGIN_API Controller::getMidiControllerAssignment(
    int32 busIndex, int16 channel, CtrlNumber midiControllerNumber,
    ParamID &id /*out*/) {

    // TODO set ID's for returned params
    if ((midiControllerNumber == kPitchBend)) {
        return kResultTrue;
    }
    if (midiControllerNumber == kCtrlModWheel) {
        return kResultTrue;
    }
    if (midiControllerNumber == kAfterTouch) {
        return kResultTrue;
    }
    if (midiControllerNumber == kCtrlAllNotesOff) {
        return kResultTrue;
    }

    return kResultFalse;
}

ParamValue PLUGIN_API
Controller::plainParamToNormalized(ParamID tag, ParamValue plainValue) {
    return EditController::plainParamToNormalized(tag, plainValue);
}

ParamValue PLUGIN_API
Controller::normalizedParamToPlain(ParamID tag, ParamValue valueNormalized) {
    return EditController::normalizedParamToPlain(tag, valueNormalized);
}

} // namespace MoniqueSynth
} // namespace Vst
} // namespace Steinberg
