import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import argparse
import time

import sys
from plots_wrapper import PlotHelpers

# locally load the config file
sys.path.insert(0, '.')
import plot_configs as plotConfigs
sys.path.remove('.')

try:
  import cPickle as pickle
except:
  import pickle

parser = argparse.ArgumentParser(description='Process pickled data to make turn-on curves')
parser.add_argument('--seedEt',       type=float, required=True, dest='seedEt_thresh',   help='seed Et Threshold (GeV)')
parser.add_argument('--noiseFilter',  type=float, required=True, dest='noise_filter',    help='gTower noise filter (GeV)')
parser.add_argument('--towerThresh',  type=float, required=True, dest='tower_thresh',    help='gTower rho Et thresh (GeV)')
parser.add_argument('--digitization', type=float, required=True, dest='digitization',    help='digitization for gTower Et (MeV)')

# parse the arguments, throw errors if missing any
args = parser.parse_args()

startTime_wall      = time.time()
startTime_processor = time.clock()

filename_id = "seed{:0.0f}_noise{:0.0f}_signal{:0.0f}_digitization{:0.0f}".format(args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed{:0.0f}/leading_jets_{}.pkl".format(args.seedEt_thresh, filename_id)
data = pickle.load(file(filename))

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished reading in data:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

dataSetStr  = plotConfigs.dataSetStr
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

helpers = PlotHelpers(dataSetStr=dataSetStr, seedCutStr=seedCutStr, noiseCutStr=noiseCutStr, towerThrStr=towerThrStr)

startTime_wall = time.time()
startTime_processor = time.clock()

# buildin plots for pt
print "pt combined"
fig, ax = pl.subplots(figsize=helpers.figsize)
n, bins, weightedPatches = ax.hist(data['oJet.pt'], weights=data['weight'], bins=np.arange(0, 500, 2), stacked=True, fill=False, histtype='step', color='b', label=r'weighted', linewidth=helpers.linewidth, alpha=0.75)
axt = ax.twinx()
n, bins, unweightedPatches = axt.hist(data['oJet.pt'], bins=np.arange(0, 500, 2), stacked=True, fill=False, histtype='step', color='r', label=r'unweighted', linewidth=helpers.linewidth, alpha=0.75)

# http://matplotlib.org/examples/api/two_scales.html
for tl in ax.get_yticklabels():
  tl.set_color('b')

for tl in axt.get_yticklabels():
  tl.set_color('r')

# http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
patches = [weightedPatches[0], unweightedPatches[0]]
# http://matplotlib.org/examples/pylab_examples/legend_auto.html
labels = [p.get_label() for p in patches]

ax.set_yscale('log', nonposy='clip')
axt.set_yscale('log', nonposy='clip')

legend = ax.legend(patches, labels, fancybox=True, framealpha=0.75, fontsize=helpers.labelsize)
legend.get_frame().set_facecolor(helpers.light_grey)
legend.get_frame().set_linewidth(0.0)

helpers.add_labels(fig, ax, xlabel=r'$p_T$ [GeV]', ylabel=r'weighted count')
helpers.add_labels(fig, axt, ylabel=r'unweighted count')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.towerThrStr])
helpers.to_file(fig, ax, "plots/weighting/{}_pt.png".format(filename_id))
pl.close(fig)

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: {:0.2f} s \n\t Clock Time: {:0.2f} s".format((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
