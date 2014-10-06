L1TriggerAnalysis
=================

Level 1 Trigger Analysis for gFEX


This uses the [`atlas_jets`](https://github.com/kratsg/atlas_jets) package that I've been writing as well, to act as a wrapper for the incoming trigger data loaded from the files.

## Scaffolding
One can quickly scaffold a new physics analysis by running
```bash
./scaffold.sh newAnalysis
```
which will set up a new directory `newAnalysis/` with the correct symlinked files and example configuration files copied over for you to edit.

## File Structure
At the top level, there are multiple folders for each of the physics analysis that is to be done. The same code should apply to all of them, apart from configuration files [datasets.yml](base/datasets.yml) and [plot_configs.py](base/plot_configs.py). The former defines the list of files for your jobs to run over, and the latter defines extra configs for [make_plots.py](base/make_plots.py). As of right now, the only configuration for plotting is the `dataSetStr` which is placed on every plot made.

To generate a Condor config file, you must run
```
ruby make_config.rb
```
which will build a config file for you. Run `ruby make_config.rb -h` for available options. Once you make your `config` file, you can submit to Condor. Every time you make a change to [datasets.yml](base/datasets.yml), you need to re-make the config file. After the jobs are done, you can merge as many files as you'd like for the plotting
```
python merge.py --*kwargs
```
where you can get a list of keyword arguments from `python merge.py -h`. Once this is done, you can just make some plots via
```
python make_plots.py --*kwargs
```
and the list of keyword arguments from `python make_plots.py -h`. Notice that a base set of keyword arguments are the same to help with consistency and copy-paste if needed.

### Submitting a job
Simply, all you need to do is edit [main.sh](base/main.sh) from any physics analysis directory (this file is symlinked). Then you can just resubmit a condor job. The job submission process is on the idea that you want to generate a job with specific configurations and then just simply `cd` into the other directories and re-submit via `condor_submit config` without any extra work. I've set this up to do exactly that. Each directory is separate because the datasets are different and the plots need to be distinct.

### Obtaining the python packages
If you run `ruby make_config.rb` and it states that it is missing `local.tar.gz` which is not included in this repository (62MB), you can grab a copy here: [faxbox::/user/kratsg/L1TriggerAnalysis/local.tar.gz](http://faxbox.usatlas.org/user/kratsg/L1TriggerAnalysis/local.tar.gz) or by XRootD
```
xrdcp root://faxbox.usatlas.org//user/kratsg/L1TriggerAnalysis/local.tar.gz local.tar.gz
```

### Setting up Dataset configurations for Condor config
In particular, I provide a [datasets.yml](base/datasets.yml) file which [make_config.rb](base/make_config.rb) expects by default (you can pass in a different file using command line flags). This is a [YAML (Yet Another Markup Language)](http://www.yaml.org/YAML_for_ruby.html) file which you can read about. There are global and inline parameters which can be set using inheritance so that inline parameters override globally set parameters. For example, if all the files are located within the same directory, set a prefix
```yaml
prefix: root://faxbox.usatlas.org//user/kratsg/LArStudies/
```
which gets preprended to all files listed below. Similarly, if most or all of your files have the same number of events, feel free to set
```yaml
numEvents: 10000
```
globally. Then for specific files, you can set `numEvents` to a different value if that particularly file has a different number.

You also specify how many jobs are run per file, which is defining how to automatically partition the file up. For example, if we set `numEvents: 10000` and then
```yaml
numJobs: 10
```
Each job will run over 1000 events in the file. Like `numEvents`, you can also set this to a different value on a file-by-file basis.

### Setting up Plot configurations
I also provide a `plot_configs.py` which only contains one variable right now, the `dataSetStr`. This is just a LaTeX string that gets rendered when making plots (such as for TTbar or ZH->nu nu bbar). Matplotlib can parse it correctly with standard LaTeX, but if you need help making this string, contact [Giordon Stark](https://github.com/kratsg) with questions, or [file an issue](https://github.com/kratsg/L1TriggerAnalysis/issues/new).

## Contact
[File an issue](https://github.com/kratsg/L1TriggerAnalysis/issues/new) or contact [Giordon Stark](https://github.com/kratsg) with any questions.
