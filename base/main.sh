#!/bin/bash 

#NEEDED
export HOME=$(pwd)
export PROOFANADIR=$(pwd)
#ROOT STUFF
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase

# run ALRB but no output
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
# set up ROOT with python-2.7
localSetupROOT 5.34.18-x86_64-slc6-gcc4.7 --skipConfirm

# extract to $HOME because that's where it looks
tar -zxf local.tar.gz -C $HOME

# export the .local directory under $HOME and prepend to $PYTHONPATH
#   by default, python should already search the directory, we're being explicit
export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages:$PYTHONPATH

# export the proxy file so that we have X509 access
#   note that this particular environment is set in the `config`
#   and that we need to export it within the current working directory
export X509_USER_PROXY=$HOME/$X509_USER_PROXY_FILENAME
echo $X509_USER_PROXY

# print directory listing just for sanity
ls -lavh

printf "Start time: "; /bin/date 
printf "Job is running on node: "; /bin/hostname 
printf "Job running as user: "; /usr/bin/id 
printf "Job is running in directory: "; /bin/pwd 

# You would normally change the --seedEt, --towerThresh, --noiseFilter, --digitization to suit your needs. Digitization is in MeV, 0 == no digitization.

python main.py --processNum ${1} --file ${2} --start ${3} --numEvents ${4} --stepSize 100 --seedEt=15 --towerThresh=6 --noiseFilter=0 --digitization=256

# Can run other scripts after processing with one set of configs

#python main.py --processNum ${1} --file ${2} --start ${3} --numEvents ${4} --stepSize 100 --seedEt=15 --towerThresh=6 --noiseFilter=0 --digitization=125
#python main.py --processNum ${1} --file ${2} --start ${3} --numEvents ${4} --stepSize 100 --seedEt=15 --towerThresh=6 --noiseFilter=0 --digitization=256
#python main.py --processNum ${1} --file ${2} --start ${3} --numEvents ${4} --stepSize 100 --seedEt=15 --towerThresh=6 --noiseFilter=0 --digitization=512



