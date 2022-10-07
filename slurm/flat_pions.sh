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

# uniform random generation
PMIN=("0" "2" "5" "30" "60" "0")
PMAX=("2" "5" "30" "60" "100" "60")
ETAMIN=("-3.5" "3.5" "-4.0")
ETAMAX=("3.5" "4.0" "-3.5")

## generate the input events
PDGID=211
echo random single pion id ${PDGID} event in p range from ${PMIN[$1]} to ${PMAX[$1]} and eta range from ${ETAMIN[$2]} to ${ETAMAX[$2]}
root -b -q "${WORK_DIR}/event_gen/single_part/gen_random_single_particle.cxx(${JUGGLER_N_EVENTS}, ${PDGID}, ${PMIN[$1]}, ${PMAX[$1]}, ${ETAMIN[$2]}, ${ETAMAX[$2]}, \"${JUGGLER_FILE_NAME_TAG}.hepmc\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script"
  exit 1
fi

echo "Running geant4 simulation"
## run geant4 simulations
ddsim --runType batch \
  --part.minimalKineticEnergy 1000*GeV  \
  --filter.tracker edep0 \
  -v WARNING \
  --numberOfEvents ${JUGGLER_N_EVENTS} \
  --compactFile ${DETECTOR_PATH}/${JUGGLER_DETECTOR_CONFIG}.xml \
  --inputFiles  ${JUGGLER_FILE_NAME_TAG}.hepmc \
  --outputFile  ${JUGGLER_SIM_FILE}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running script"
  exit 1
fi

rootls -t ${JUGGLER_SIM_FILE}

if [[ -z "${ANALYSIS_ONLY}" ]] ;
then
  # Need to figure out how to pass file name to juggler from the commandline
  # gaudirun.py ${RUN_DIR}/benchmarks/tracking/options/track_reconstruction.py
  # gaudirun.py ${RUN_DIR}/benchmarks/tracking/options/track_reconstruction_YS.py
  gaudirun.py /opt/benchmarks/reconstruction_benchmarks/benchmarks/tracking/options/track_reconstruction.py
  #gaudirun.py ${RUN_DIR}/benchmarks/tracking/options/track_reconstruction_ePIC.py
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running juggler"
    exit 1
  fi
fi

