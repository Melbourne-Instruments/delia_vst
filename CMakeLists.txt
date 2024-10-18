project(delia_vst)
cmake_minimum_required(VERSION 3.12)

# Add c++ wave reading lib
add_subdirectory("${PROJECT_SOURCE_DIR}/AudioFile" SYSTEM)

# Xenomai
set(XENOMAI_BASE_DIR "/usr/xenomai" CACHE STRING "xenomai base dir path")

# Check if a cross-compile build
message("CMAKE_SYSROOT is " ${CMAKE_SYSROOT})
if(NOT "${CMAKE_SYSROOT}" STREQUAL "")
    set(CMAKE_SYSROOT "${CMAKE_SYSROOT}")
    set(CMAKE_FIND_ROOT_PATH ${CMAKE_SYSROOT})
    set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
    set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
    set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
    set(XENOMAI_BASE_DIR "${CMAKE_SYSROOT}/usr/xenomai")
endif()
message("XENOMAI_BASE_DIR is " ${XENOMAI_BASE_DIR})

option(WITH_XENOMAI "Enable Xenomai support" ON)

function(add_xenomai_to_target target)
    # from `xeno-config --skin=posix --cflags`
    set(XENOMAI_C_FLAGS "-D_REENTRANT -D__COBALT__ -D__COBALT_WRAP__")

    set(XENOMAI_INCLUDE_DIRS
        ${XENOMAI_BASE_DIR}/include
        ${XENOMAI_BASE_DIR}/include/cobalt
    )

    find_library(COBALT_LIB cobalt HINTS ${XENOMAI_BASE_DIR}/lib)
    target_compile_options(${target} PRIVATE ${XENOMAI_C_FLAGS})
    target_include_directories(${target} PRIVATE ${XENOMAI_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${COBALT_LIB} rt m)
endfunction()

set(VST3_SDK_ROOT "${PROJECT_SOURCE_DIR}/vst3sdk" CACHE PATH "Location of VST3 SDK")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(AddVST3SDK)

smtg_get_linux_architecture_name() # Sets var LINUX_ARCHITECTURE_NAME

set(monique-vst3_sources
factory.cpp
FastApprox.h
MoniqueCommon.h
MoniqueController.cpp
MoniqueController.h
MoniqueProcessor.cpp
MoniqueProcessor.h
MoniqueDriveCompensator.h
MoniqueEnvelope.h
MoniqueSlaveEnvelope.h
MoniqueLayer.h
MoniqueDriftSource.h
MoniqueLfo.h
MoniqueMatrix.h
MoniqueModel.h
MoniqueOscMix.h
MoniqueOutputPanner.h
MoniqueVoice.h
MoniqueWavetable.h
MoniqueWavetable.cpp
MoniqueLayerManager.h
MoniqueAnalogModel.h
MoniqueFilters.h
NoiseOscillator.h
SynthMath.h
MoniqueUnitTests.h
MoniqueCalibrator.h
GuiMsgQueue.h
MoniqueEffectsEngine.h
MoniqueExpander.h
MoniqueReverb.h
MoniqueReverb.cpp
Chorus.h
ChorusEngine.h
OnePoleLP.h
MoniqueDelay.h
MoniqueSampler.h
freeverb/allpass.cpp
freeverb/allpass.hpp
freeverb/allpass_t.hpp
freeverb/biquad.cpp
freeverb/biquad.hpp
freeverb/biquad_t.hpp
freeverb/comb.cpp
freeverb/comb.hpp
freeverb/comb_t.hpp
freeverb/delay.cpp
freeverb/delay.hpp
freeverb/delayline.cpp
freeverb/delayline.cpp
freeverb/delayline.hpp
freeverb/delayline_t.hpp
freeverb/earlyref.cpp
freeverb/earlyref.hpp
freeverb/earlyref_t.hpp
freeverb/efilter.cpp
freeverb/efilter.hpp
freeverb/efilter_t.hpp
freeverb/nrev.cpp
freeverb/nrev.hpp
freeverb/nrev_t.hpp
freeverb/nrevb.cpp
freeverb/nrevb.hpp
freeverb/nrevb_t.hpp
freeverb/progenitor.cpp
freeverb/progenitor.hpp
freeverb/progenitor_t.hpp
freeverb/progenitor2.cpp
freeverb/progenitor2.hpp
freeverb/progenitor2_t.hpp
freeverb/revbase.cpp
freeverb/revbase.hpp
freeverb/revbase_t.hpp
freeverb/slot.cpp
freeverb/slot.hpp
freeverb/slot_t.hpp
freeverb/strev.cpp
freeverb/strev.hpp
freeverb/strev_t.hpp
freeverb/utils.cpp
freeverb/utils.hpp
freeverb/utils_t.hpp
freeverb/zrev.cpp
freeverb/zrev.hpp
freeverb/zrev_t.hpp
freeverb/zrev2.cpp
freeverb/zrev2.hpp
freeverb/zrev2_t.hpp
test.cpp
FactorySampleMsgQueue.cpp
FactorySampleMsgQueue.h
FactorySamplerController.cpp
FactorySamplerController.h
FactorySamplerProcessor.cpp
FactorySamplerProcessor.h
FactoryTestController.cpp
FactoryTestController.h
FactoryTestProcessor.cpp
FactoryTestProcessor.h
version.h
)
FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)
PREPEND( monique-vst3_sources "source" "${monique-vst3_sources}")
set(INCLUDE_DIRS
    "${PROJECT_SOURCE_DIR}/source"
    "${PROJECT_SOURCE_DIR}/source/freeverb"
    "${PROJECT_SOURCE_DIR}/delia_common/include"
    "${PROJECT_SOURCE_DIR}/submodules/AudioFile")

message("make unit tests?")
message(STATUS "${LINUX_ARCHITECTURE_NAME}")

execute_process(
    COMMAND git log -1 --format=%H
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(SUBSTRING "${GIT_COMMIT_HASH}" 0 8 GIT_SHORT_HASH)
message("${GIT_SHORT_HASH}")
set(VST_X_VER "1")
set(VST_Y_VER "0")
set(VST_Z_VER "3")
configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.in
    ${CMAKE_SOURCE_DIR}/source/version.h
)

configure_file(
    ${CMAKE_SOURCE_DIR}/version.txt.in
    ${CMAKE_BINARY_DIR}/VST3/Debug/delia_vst.vst3/Contents/version.txt
)
if(NOT(LINUX_ARCHITECTURE_NAME STREQUAL "aarch64-linux"))
    message("make and run native tests")
    add_executable(moniqueUnitTesting ${monique-vst3_sources})
    target_compile_definitions(moniqueUnitTesting PRIVATE MONIQUE_NATIVE=1)

    target_compile_options(moniqueUnitTesting PRIVATE -std=c++20 -O0 -Wall -Wextra -Wno-psabi)
    target_include_directories(moniqueUnitTesting PRIVATE ${INCLUDE_DIRS})

    target_link_libraries(moniqueUnitTesting PRIVATE base -Bstatic sdk AudioFile fftw3 rt -std=c++20 -O0 )
    add_custom_target(run_unit_tests
        COMMAND moniqueUnitTesting
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}

        COMMENT "run generated tests in ${CMAKE_CURRENT_BINARY_DIR}"
    )
endif()

set(target delia_vst)
set(VST_INC_DIRS "${VST3_SDK_ROOT}" "${INCLUDE_DIRS}")
set(VST3_OUTPUT_DIR ${CMAKE_BINARY_DIR}/VST3)
smtg_add_vst3plugin(${target} ${VST3_SDK_ROOT} ${monique-vst3_sources})
target_include_directories(${target} PRIVATE ${VST_INC_DIRS})
if (${WITH_XENOMAI})
    message("Building with Xenomai support")
    add_xenomai_to_target(${target})
else()
    message("Building WITHOUT Xenomai support")
    target_compile_definitions(${target} PRIVATE NO_XENOMAI)
endif()
target_link_libraries(${target} PRIVATE base sdk AudioFile pthread)
message("here ${LINUX_ARCHITECTURE_NAME}")

if((LINUX_ARCHITECTURE_NAME STREQUAL "aarch64-linux"))
    message("set arm compile options")
    target_compile_options(delia_vst PRIVATE -mtune=cortex-a72 -funsafe-math-optimizations -std=c++20 -ftree-vectorize -flto -Ofast -ffast-math)

else()
    message("set 86 compile options")
    target_compile_definitions(moniqueUnitTesting PRIVATE MONIQUE_NATIVE=1)
    target_compile_options(delia_vst PRIVATE -Wall -Wextra -Wno-psabi -std=c++20 -O0)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math -ftree-vectorize")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math -fno-rtti -ftree-vectorize -std=c++20")

file(GLOB_RECURSE ALL_CXX_SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/source/*.cpp
    ${PROJECT_SOURCE_DIR}/source/*.h
    ${PROJECT_SOURCE_DIR}/source/*.h
    ${PROJECT_SOURCE_DIR}/source/freeverb/*.h
    ${PROJECT_SOURCE_DIR}/source/freeverb/*.hpp
    ${PROJECT_SOURCE_DIR}/source/freeverb/*.cpp    
)

# generate formatting target
add_custom_target(
    format
    COMMAND /usr/bin/clang-format-14
    -i
    -style=file
    ${ALL_CXX_SOURCE_FILES}
)

add_custom_command(TARGET delia_vst POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/source/RunCalibration.py ${CMAKE_BINARY_DIR}/VST3/Debug/delia_vst.vst3/Contents/)
add_custom_command(TARGET delia_vst POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/source/FilterTuningRoutine.py ${CMAKE_BINARY_DIR}/VST3/Debug/delia_vst.vst3/Contents/)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_custom_target(lint
    COMMAND sh -c "clang-tidy ${ALL_CXX_SOURCE_FILES}"
    VERBATIM
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)