#pragma once
#include "monique_synth_parameters.h"
#include "pluginterfaces/base/ustring.h"
#include "pluginterfaces/vst/ivstmidicontrollers.h"
#include "public.sdk/source/vst/vsteditcontroller.h"
#include "public.sdk/source/vst/vstparameters.h"
#include <stdio.h>

namespace Steinberg {
namespace Vst {
namespace MoniqueSynth {

//-----------------------------------------------------------------------------
class Controller : public EditController, public IMidiMapping {
  public:
    // --- EditController Overrides
    Controller() { printf("\ncreate Monique Controller\n"); }

    ~Controller(){};
    tresult PLUGIN_API initialize(FUnknown *context);
    tresult PLUGIN_API terminate();

    ParamID PLUGIN_API firstLayerParamId();
    ParamID PLUGIN_API firstLayerStateParamId();
    int32 PLUGIN_API getNumLayerParams();
    int32 PLUGIN_API getNumLayerStateParams();

    // --- serialize-read from file
    tresult PLUGIN_API setComponentState(IBStream *fileStream);

    // --- IMidiMapping
    virtual tresult PLUGIN_API getMidiControllerAssignment(
        int32 busIndex, int16 channel, CtrlNumber midiControllerNumber,
        ParamID &id /*out*/);

    // override this function to mask off the layer state bits of the param ID, this is so we return the actual parameter
    virtual Parameter *getParameterObject(ParamID tag) override {
        return parameters.getParameter(tag);
    }

    // --- oridinarily not needed; see documentation on Automation for using these
    virtual ParamValue PLUGIN_API
    normalizedParamToPlain(ParamID id, ParamValue valueNormalized);
    virtual ParamValue PLUGIN_API plainParamToNormalized(ParamID id,
        ParamValue plainValue);

    // --- Our COM Creating Method
    static FUnknown *createInstance(void *) {
        return (IEditController *)new Controller;
    }

    // --- our globally unique ID value
    static FUID uid;

    // --- helper function for serialization
    tresult PLUGIN_API setParamNormalizedFromFile(ParamID tag, ParamValue value);

    // --- define the controller and interface
    OBJ_METHODS(Controller, EditController)
    DEFINE_INTERFACES
    DEF_INTERFACE(IMidiMapping)
    END_DEFINE_INTERFACES(EditController)
    REFCOUNT_METHODS(EditController)
};

} // namespace MoniqueSynth
} // namespace Vst
} // namespace Steinberg
