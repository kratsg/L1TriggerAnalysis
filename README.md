L1TriggerAnalysis
=================

Level 1 Trigger Analysis for gFEX


This uses the [`atlas_jets`](https://github.com/kratsg/atlas_jets) package that I've been writing as well, to act as a wrapper for the incoming trigger data loaded from the files.

## File Structure
At the top level, there are multiple folders for each of the physics analysis that is to be done. The same code should apply to all of them, apart from configuration files `datasets` and `plots_config.py`. The former defines the list of files for your jobs to run over, and the latter defines extra configs for `make_plots.py`. As of right now, the only configuration for plotting is the `dataSetStr` which is placed on every plot made.

To generate a Condor config file, you must run
```
ruby make_config.rb
```
which will build a config file for you. Run `ruby make_config.rb -h` for available options. Afterwards, you can merge as many files as you'd like for the plotting
```
python merge.py --*kwargs
```
where you can get a list of keyword arguments from `python merge.py -h`. Similarly, you can make your plots via
```
python make_plots.py --*kwargs
```
and the list of keyword arguments from `python make_plots.py -h`. Notice that a base set of keyword arguments are the same to help with consistency and copy-paste if needed.

### Contact

File an issue or contact [Giordon Stark](https://github.com/kratsg) with any questions.
