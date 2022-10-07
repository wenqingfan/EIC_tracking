#!/bin/bash

print_env.sh

## To run the reconstruction, we need the following global variables:
## - JUGGLER_INSTALL_PREFIX:   Install prefix for Juggler (simu/recon)
## - JUGGLER_DETECTOR:         the detector package we want to use for this benchmark
## - JUGGLER_DETECTOR_VERSION: the detector package we want to use for this benchmark
## - DETECTOR_PATH:            full path to the detector definitions
##
## You can ready options/env.sh for more in-depth explanations of the variables
## and how they can be controlled.
export DETECTOR_PATH=${DETECTOR_PATH}

if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=10
fi

export JUGGLER_FILE_NAME_TAG="single_pions"
export JUGGLER_GEN_FILE="${JUGGLER_FILE_NAME_TAG}.hepmc"

export JUGGLER_SIM_FILE="sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="rec_${JUGGLER_FILE_NAME_TAG}.root"

echo "JUGGLER_N_EVENTS = ${JUGGLER_N_EVENTS}"
echo "JUGGLER_DETECTOR = ${JUGGLER_DETECTOR}"

export RUN_DIR=/global/project/projectdirs/m3763/wenqing/eic/reconstruction_benchmarks/
export WORK_DIR=/global/project/projectdirs/m3763/wenqing/eic/

rootls -t ${JUGGLER_SIM_FILE}

if [[ -z "${ANALYSIS_ONLY}" ]] ;
then
  # Need to figure out how to pass file name to juggler from the commandline
  gaudirun.py ${RUN_DIR}/benchmarks/tracking/options/track_reconstruction_ePIC.py
  # gaudirun.py ${RUN_DIR}/benchmarks/tracking/options/track_reconstruction_YS.py
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running juggler"
    exit 1
  fi
fi

