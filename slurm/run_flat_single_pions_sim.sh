#!/bin/bash

eicdir=/global/project/projectdirs/m3763/wenqing/eic/

echo start jobs in ${eicdir}
cd ${eicdir}

source eic-shell -n <<< "pwd;
    source mysetup.sh;
    echo juggler version is ${JUGGLER_INSTALL_PREFIX};
    bash flat_single_pions_sim.sh $1 $2 $3;
    popd
    "
