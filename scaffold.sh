#!/bin/bash
if [ $# -ne 1 ]; then
  # Control will enter here if num params != 2
  printf "Usage:\n\t$0 <directory>\n" # not one parameter
  exit 2
fi
if [ -d "$1" ]; then
  # Control will enter here if $1 exists
  echo "'$1' already exists!"
  exit 2
fi

# one parameter, directory name, does not exist
mkdir $1 && cd $1
# symlink files to set it up
echo "Symlinking main.py..."
ln -s ../base/main.py main.py
echo "Symlinking main.sh..."
ln -s ../base/main.sh main.sh
echo "Symlinking make_config.rb..."
ln -s ../base/make_config.rb make_config.rb
echo "Symlinking make_plots.py..."
ln -s ../base/make_plots.py make_plots.py
echo "Symlinking merge.py..."
ln -s ../base/merge.py merge.py
echo "Symlinking plots_wrapper.py..."
ln -s ../base/plots_wrapper.py plots_wrapper.py
echo "Adding example dataset configuration file for make_config.rb..."
cp ../base/datasets.yml datasets.yml
echo "Adding example plot configuration file for plot_configs.rb..."
cp ../base/plot_configs.py plot_configs.py
echo "Creating output directories for condor jobs..."
mkdir data plots out
