#!/bin/bash

mkdir -p output_sim
pushd output_sim # create sub dir under outdir to avoid conflict if running multiple setting the same time

# create different running directory for different jobs
INPUT=$(( $1 + 0 ))
JOB=`printf "%05d" $INPUT`
echo $JOB
mkdir -p $JOB
pushd $JOB

cp /global/project/projectdirs/m3763/wenqing/eic/reconstruction_benchmarks/benchmarks/tracking/flat_pions.sh .

export JUGGLER_N_EVENTS=1000

# $2 -- mom range, $3 -- eta range
bash flat_pions.sh $2 $3
echo check current directory
pwd
echo list files under current directory
ls -lhtr

mv single_pions.hepmc ../single_pions_${INPUT}.hepmc
mv sim_single_pions.edm4hep.root ../sim_single_pions_${INPUT}.root
mv rec_single_pions.root ../rec_single_pions_${INPUT}.root

popd
rm -rf $JOB
popd
