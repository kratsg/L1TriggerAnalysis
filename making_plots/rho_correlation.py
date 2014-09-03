from atlas_jets import *
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import os
import argparse
import time

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

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

#this is a wrapper around file strings to ensure the directory exists
#       could use decorators...
def write_file(f):
  ensure_dir(f)
  return f

filename_id = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "ZH_leading_jets_%s.pkl" % (filename_id)

print "loading"

data = pickle.load(file(filename))

print "loaded"

bins_rho = np.arange(0.,100.,1.)
bins_vertices = np.arange(0.,100.,1.)

figsize = (16, 12)
labelsize = 30
titlesize = 36
linewidth = 4
light_grey = np.array([float(200)/float(255)]*3)
filled=False
cmap = 'jet'

dataSetStr  = '$ZH \\to \\nu\\nu b\\bar{b}$\n$\sqrt{s}=14 \ \mathrm{TeV}\ \langle\mu\\rangle=80$'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

label_oRho = r'Offline $\rho$ [GeV/area]'
label_gRho = r'gFEX $\rho$ [GeV/area]'

region_legend = {'all': r'towers: $\eta \in (-4.9,4.9)$',\
                1: r'towers: $\eta \in [-1.6,0.0)$',\
                2: r'towers: $\eta \in [0.0,1.6)$',\
                3: r'towers: $\eta \in (-4.9,-1.6)$',\
                4: r'towers: $\eta \in [1.6,4.9)$',\
                '3a': r'towers: $\eta \in (-2.5,-1.6)$',\
                '3b': r'towers: $\eta \in (-3.2,-2.5)$',\
                '3c': r'towers: $\eta \in (-4.9,-3.2)$',\
                '4a': r'towers: $\eta \in (1.6,2.5)$',\
                '4b': r'towers: $\eta \in (2.5,3.2)$',\
                '4c': r'towers: $\eta \in (3.2,4.9)$'}

from matplotlib.ticker import FuncFormatter

def latex_float(f, pos=0):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:
        return r"${}$".format(float_str)


formatter = FuncFormatter(latex_float)

startTime_wall      = time.time()
startTime_processor = time.clock()

corr = np.corrcoef(data['offline_rho'], data['gFEX_rho_1'])[0,1]
fig = pl.figure(figsize=figsize)
hist, _, _, _ = pl.hist2d(data['offline_rho'], data['gFEX_rho_1'], norm=LogNorm(), bins=(bins_rho, bins_rho) , alpha=0.75, cmap = cmap)

cbar = pl.colorbar(ticks=np.logspace(0, np.log10(np.max(hist)), 10), format=formatter)#'%.2e')
cbar.set_label('number density', fontsize=labelsize)
cbar.ax.tick_params(labelsize=labelsize-6)

pl.xlabel(label_oRho, fontsize=labelsize)
pl.ylabel('%s, %s' % (label_gRho, region_legend[1]), fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.text(0.95, 0.05, '%s\n%s\n$\mathrm{Corr} = %0.4f$' % (dataSetStr, towerThrStr, corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)

#add atlas simulation
#    internal
pl.text(0.05, 0.95, 'ATLAS', fontsize=42, style='italic', fontweight='bold', verticalalignment='top', horizontalalignment='left', transform=fig.gca().transAxes)
pl.text(0.27, 0.95, 'Preliminary', verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)
pl.text(0.05, 0.90, 'Simulation', verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)

pl.tick_params(axis='both', which='major', labelsize=labelsize-6)
pl.savefig( write_file('ZH/%s_%s_correlation.png' % (filename_id, 'gFEX_rho_1')), bbox_inches='tight')
pl.savefig( write_file('ZH/%s_%s_correlation.eps' % (filename_id, 'gFEX_rho_1')), bbox_inches='tight')
pl.savefig( write_file('ZH/%s_%s_correlation.pdf' % (filename_id, 'gFEX_rho_1')), bbox_inches='tight')
pl.close()



endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
