#!/bin/bash

mkdir -p output_rec_hi
pushd output_rec_hi # create sub dir under outdir to avoid conflict if running multiple setting the same time

# create different running directory for different jobs
INPUT=$(( $1 + 0 ))
JOB=`printf "%05d" $INPUT`
echo $JOB
mkdir -p $JOB
pushd $JOB

cp /global/project/projectdirs/m3763/wenqing/eic/reconstruction_benchmarks/benchmarks/tracking/flat_pions_rec_only.sh .

export JUGGLER_N_EVENTS=1000

cp /global/project/projectdirs/m3763/wenqing/eic/output_sim_hi/sim_single_pions_${INPUT}.root sim_single_pions.edm4hep.root
bash flat_pions_rec_only.sh
echo check current directory
pwd
echo list files under current directory
ls -lhtr

mv rec_single_pions.root ../rec_single_pions_${INPUT}.root

popd
rm -rf $JOB
popd
