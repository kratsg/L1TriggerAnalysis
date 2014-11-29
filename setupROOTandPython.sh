#!/bin/bash
echo "Setting up environment with ROOT=5.34.18 and Python=2.7 and exporting packages at .local"
localSetupROOT 5.34.18-x86_64-slc6-gcc4.7
export PYTHONPATH=$(pwd)/.local/lib/python2.7/site-packages:$PYTHONPATH

