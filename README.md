# DELIA VST #

DELIA VST app for the Raspberry Pi (64-bit).

### Building the DELIA VST ###

To cross compile for the Raspberry Pi, the Melbourne Instruments SDK *must* be used.
It is recommended to install the SDK in the /opt folder on your Ubuntu PC.

Once this has been done, source the environment script to set-up the build environment, for example:

$ source /opt/monique/1.0.0/environment-setup-cortexa72-elk-linux

create build folder and cmake script:
$ mkdir build
$ cd build
$ cmake ../ -DCMAKE_BUILD_TYPE=Debug #for debug build

build
$ make

format before pushing
$ cd build && make format

### Dependencies ###


---
Copyright 2023-2024 Melbourne Instruments, Australia.
