/**
 * @file factory.cpp
 * @brief
 * @version 0.1
 * @date 2021-09-30
 *
 * @copyright Copyright (c) 2023 Melbourne Instruments
 *
 */

#include "FactorySamplerController.h"
#include "FactorySamplerProcessor.h"
#include "FactoryTestController.h"
#include "FactoryTestProcessor.h"
#include "MoniqueController.h"
#include "MoniqueProcessor.h"
#include "public.sdk/source/main/pluginfactoryvst3.h"
#include "version.h"

//-----------------------------------------------------------------------------
bool InitModule() { return true; }

bool DeinitModule() { return true; }

//-----------------------------------------------------------------------------
#define kVersionString FULL_VERSION_STR

BEGIN_FACTORY_DEF(stringCompanyName, "", "")

//-----------------------------------------------------------------------------
// -- Monique VST
DEF_CLASS2(INLINE_UID_FROM_FUID(Vst::MoniqueSynth::Processor::uid),
    PClassInfo::kManyInstances, kVstAudioEffectClass, "Delia Vst",
    Vst::kDistributable, Vst::PlugType::kInstrumentSynth, kVersionString,
    kVstVersionString, Vst::MoniqueSynth::Processor::createInstance)

DEF_CLASS2(INLINE_UID_FROM_FUID(Vst::MoniqueSynth::Controller::uid),
    PClassInfo::kManyInstances, kVstComponentControllerClass,
    "Delia VST", Vst::kDistributable, "", kVersionString,
    kVstVersionString, Vst::MoniqueSynth::Controller::createInstance)

// -- Factory Test VST's
DEF_CLASS2(INLINE_UID_FROM_FUID(Vst::MoniqueFactory::FactoryTestProcessor::uid),
    PClassInfo::kManyInstances, kVstAudioEffectClass, "factory test",
    Vst::kDistributable, Vst::PlugType::kInstrumentSynth, kVersionString,
    kVstVersionString, Vst::MoniqueFactory::FactoryTestProcessor::createInstance)

DEF_CLASS2(INLINE_UID_FROM_FUID(Vst::MoniqueFactory::FactoryTestController::uid),
    PClassInfo::kManyInstances, kVstComponentControllerClass,
    "factory test", Vst::kDistributable, "", kVersionString,
    kVstVersionString, Vst::MoniqueFactory::FactoryTestController::createInstance)
DEF_CLASS2(INLINE_UID_FROM_FUID(Vst::MoniqueFactory::FactorySamplerProcessor::uid),
    PClassInfo::kManyInstances, kVstAudioEffectClass, "factory sampler",
    Vst::kDistributable, Vst::PlugType::kAnalyzer, kVersionString,
    kVstVersionString, Vst::MoniqueFactory::FactorySamplerProcessor::createInstance)
DEF_CLASS2(INLINE_UID_FROM_FUID(Vst::MoniqueFactory::FactorySamplerController::uid),
    PClassInfo::kManyInstances, kVstComponentControllerClass,
    "factory sampler", Vst::kDistributable, "", kVersionString,
    kVstVersionString, Vst::MoniqueFactory::FactorySamplerController::createInstance) //-----------------------------------------------------------------------------
END_FACTORY
