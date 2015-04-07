#!/usr/bin/env bash
# For setting up and compiling Release version on kongull

source environment.sh

rm -rf Release
mkdir Release
cd Release

CXX=icpc CC=icc FC=ifort cmake -DCMAKE_BUILD_TYPE=Release ..

make
