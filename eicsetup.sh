#!/bin/sh

# accessing ATHENA data
export S3_ACCESS_KEY=eicS3read
export S3_SECRET_KEY=eicS3read

#source /opt/detector/setup.sh
source /opt/detector/athena-nightly/setup.sh

# alias
alias mywork='cd /global/project/projectdirs/alice/$USER/'
alias eic='cd /global/project/projectdirs/m3763/$USER/'
alias sq='squeue -u wenqing'

# useful stuff
alias ll='ls -lhrt'
alias c='clear'
alias rm='rm -i'
