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
filename = "plots/subjets/{}.pkl".format(filename_id)

signal = pickle.load(file('../ttbar/{}'.format(filename)))
background = pickle.load(file('../dijets/{}'.format(filename)))

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished reading in data:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

dataSetStr  = plotConfigs.dataSetStr
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

helpers = PlotHelpers(dataSetStr=dataSetStr, seedCutStr=seedCutStr, noiseCutStr=noiseCutStr, towerThrStr=towerThrStr, labelsize=28)

PtCutLabel = "$250\ \mathrm{GeV}\ \leq P_T^\mathrm{oJet} \leq\ 500\ \mathrm{GeV}$"
MassCutLabel = "$100\ \mathrm{GeV}\ \leq m^\mathrm{oJet} \leq\ 220\ \mathrm{GeV}$"

dataStorage = {}


def rescaleY(y, option="fake rate"):
  validOptions = ["fake rate", "rejection", "purity"]
  if option == "fake rate":
    return y
  elif option == "rejection":
    return 1./y
  elif option == "purity":
    return 1.-y
  else:
    raise ValueError("option does not have a valid value; {}".format(validOptions))

# we want to apply offline cuts to the signal, but not dijets -- both PtCut and
# MassCut

print "Making rates for subjet studies"
# each file will be based on a gJetEt_cut
for gJetEt_cut in [200., 250., 300., 350.]:
  print "\t"*1, "gJetEt_cut = {}".format(gJetEt_cut)
  gJetEtCutLabel = '$E_T^\mathrm{{gFEX\ jet}} \geq {:0.2f}\ \mathrm{{GeV}}$'.format(gJetEt_cut)

  offlineLabel = "\nOffline Selections:\n\t{}\n\t{}".format(PtCutLabel, MassCutLabel)
  filenameEnd = '_250oJetPt500_100oJetM220'
  for gTowerEt_cut in [5., 10., 15., 20., 25.]:
    print "\t"*2, "gTowerEt_cut = {}".format(gTowerEt_cut)
    # make a label for gTowerEt cut
    gTowerEtThrLabel = r'$N\left(E_T^\mathrm{{gFEX\ tower}} >\ {:0.2f}\ \mathrm{{GeV}}\right) \geq X$'.format(gTowerEt_cut)

    sigData = np.zeros((4, 4))

    bkgData = background['gJetEtCut_{:0.0f}'.format(gJetEt_cut)]['']['gTowerEtCut_{:0.0f}'.format(gTowerEt_cut)]['oJetnsjCut_{:d}'.format(0)]['vals']

    for gJetnsj_cut, i in zip([1, 2, 3, 4], range(4)):
      # add storage by oJetnsj_cut
      sigData[i] = signal['gJetEtCut_{:0.0f}'.format(gJetEt_cut)][filenameEnd]['gTowerEtCut_{:0.0f}'.format(gTowerEt_cut)]['gJetnsjCut_{:d}'.format(gJetnsj_cut)]['vals']

    fig, ax = pl.subplots(figsize=helpers.figsize)
    for sig, color, oJetnsj_cut in zip(sigData.T, helpers.colors, range(1, 5)):
      print "\t"*3, "oJetnsj_cut = {}".format(oJetnsj_cut)
      # make a label for oJet nsj cut
      oJetnsjCutLabel = r'$N(P_T^\mathrm{{oJet\ subjet}} >\ 20 \ \mathrm{{GeV}}) \geq {:d}$'.format(oJetnsj_cut)
      ax.plot(sig, bkgData, linestyle='steps-mid', alpha=0.75, color=color, marker='x', ms=20, mew=10, linewidth=0, label=oJetnsjCutLabel)

    ax.set_xlim((0.0, 1.1))
    ax.set_ylim((0.0, 2.0))
    ax.set_xticks(np.linspace(0., 1., 6))
    ax.set_yticks(np.linspace(0., 1., 6))
    helpers.add_legend(fig, ax, numpoints=1)
    helpers.add_labels(fig, ax, xlabel='signal', ylabel='background')
    helpers.add_description(fig, ax, align='tl', strings=[helpers.dataSetStr, 'iso., $\Delta R(\mathrm{gJet},\mathrm{oJet})\leq 1$', helpers.seedCutStr, helpers.noiseCutStr, gJetEtCutLabel, gTowerEtThrLabel, offlineLabel])
    helpers.add_grid(fig, ax)
    helpers.to_file(fig, ax, 'rates/{}_gJetEt{:0.0f}_gTowerEt{:0.0f}.png'.format(filename_id, gJetEt_cut, gTowerEt_cut))
    pl.close(fig)

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: {:0.2f} s \n\t Clock Time: {:0.2f} s".format((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
