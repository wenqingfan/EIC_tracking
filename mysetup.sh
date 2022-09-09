#!/bin/sh

# accessing ATHENA data
export S3_ACCESS_KEY=eicS3read
export S3_SECRET_KEY=eicS3read

export EIC_SHELL_PREFIX=/global/project/projectdirs/m3763/$USER/eic/local

#source /opt/detector/setup.sh
source /opt/detector/epic-nightly/setup.sh
export DETECTOR_PATH=$EIC_SHELL_PREFIX/../epic
export JUGGLER_INSTALL_PREFIX=$EIC_SHELL_PREFIX
export LD_LIBRARY_PATH=$EIC_SHELL_PREFIX/lib:$LD_LIBRARY_PATH
echo detector path is ${DETECTOR_PATH}
echo juggler version is ${JUGGLER_INSTALL_PREFIX}

#source /opt/detector/setup.sh
#
#export ATHENA_PREFIX=/global/project/projectdirs/m3763/wenqing/eic/local
#
## echo using my own version of juggler by default
#export JUGGLER_INSTALL_PREFIX=${ATHENA_PREFIX}
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ATHENA_PREFIX}/lib
#echo juggler version is ${JUGGLER_INSTALL_PREFIX}

# alias
alias mywork='cd /global/project/projectdirs/alice/$USER/'
alias eic='cd /global/project/projectdirs/m3763/$USER/'
alias sq='squeue -u wenqing'

# useful stuff
alias ll='ls -lhrt'
alias c='clear'
alias rm='rm -i'
