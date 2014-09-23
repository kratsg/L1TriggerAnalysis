import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib.patches import Rectangle
import os
import argparse

from scipy.interpolate import UnivariateSpline

try:
  import cPickle as pickle
except:
  import pickle

parser = argparse.ArgumentParser(description='Process pickled data to make turn-on curves')
parser.add_argument('-o', '--noiseFilter', type=float, required=True, dest='noise_filter',    help='gTower noise filter')
parser.add_argument('-d', '--seedEt', type=float, required=True, dest='seedEt_thresh', help='seed Et Threshold')

tower_thresholds = [3,4,6]

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

def FWHM(bins, vals):
  spline = UnivariateSpline(bins[:-1]+np.diff(bins)/2., vals-np.max(vals)/2., s=0)
  r1, r2 = spline.roots() # find the roots
  return np.abs(r1-r2)

bins_rho = np.arange(0.,100.,1.)
bins_vertices = np.arange(0.,100.,1.)

figsize = (16, 12)
labelsize = 28
titlesize = 36
linewidth = 1
light_grey = np.array([float(200)/float(255)]*3)
filled=False

basicTextStr = 'TTbar 14TeV $\langle\mu\\rangle=80$\n$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$\n$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % (args.seedEt_thresh, args.noise_filter)
textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

artists = []
labels = []
fwhmVals = []

fig = pl.figure(figsize=figsize)
for tower_thresh, color in zip(tower_thresholds,['b','r','c','m','y','k']):
  filename = "plots/pickle/seed%d_noise%d_signal%d_resolution_PtOffline_projection_noPileup.pkl" % (args.seedEt_thresh, args.noise_filter, tower_thresh)
  data = pickle.load(file(filename))
  hist, bins = np.histogram(data['170to180'], bins=100, density=True)
  fwhmVals.append('FWHM: %0.4f' % FWHM(bins, hist) )
  pl.plot(bins[:-1], hist, linestyle='steps-post', marker='o', alpha=0.75, markersize=10, color=color, linewidth=linewidth)
  hist, bins = np.histogram(data['200to220'], bins=100, density=True)
  fwhmVals.append('FWHM: %0.4f' % FWHM(bins, hist) )
  pl.plot(bins[:-1]+1.0, hist, linestyle='steps-post', marker='s', alpha=0.75, markersize=10, color=color, linewidth=linewidth)
  hist, bins = np.histogram(data['300to350'], bins=100, density=True)
  fwhmVals.append('FWHM: %0.4f' % FWHM(bins, hist) )
  pl.plot(bins[:-1]+2.0, hist, linestyle='steps-post', marker='v', alpha=0.75, markersize=10, color=color, linewidth=linewidth)
  # each of these returns a list of lines, but there's only one LOL
  artists.append(Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='none', linewidth=0))
  labels.append(r'$\rho(%d\ \mathrm{GeV})$' % tower_thresh)

pl.axvline(0., linewidth=4, color='k', alpha=0.75, zorder=1)
pl.axvline(1., linewidth=4, color='k', alpha=0.75, zorder=1)
pl.axvline(2., linewidth=4, color='k', alpha=0.75, zorder=1)

def dummyArtist(marker):
  return pl.Line2D((0,1),(0,0), color='k', marker=marker, linestyle='none', markersize=15 )

def blankArtist(color):
  return Rectangle((0,0), 1, 1, facecolor=color, edgecolor='none', linewidth=0)

dummyArtists = [dummyArtist('o'), dummyArtist('s'), dummyArtist('v')]
dummyLabels = [r'$170\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 180\ \mathrm{GeV}$', r'$200\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 220\ \mathrm{GeV}$', r'$300\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 350\ \mathrm{GeV}$']
artists = np.array(dummyArtists + artists + [blankArtist(color) for color in ['b','b','b','r','r','r','c','c','c']])
labels = np.array(dummyLabels + labels + fwhmVals)
#fwhmVals = np.array(fwhmVals).reshape(3,-1).T.flatten()

#reshape based on number of rows you want
legend = pl.legend(artists.reshape(5,-1).T.flatten(), labels.reshape(5,-1).T.flatten(), fancybox=True, framealpha=0.75, fontsize=labelsize*2/3, numpoints=1, ncol=3, mode='expand')
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$, $p_T$ bins offset in $x$', fontsize=labelsize)
pl.ylabel('normalized counts', fontsize=labelsize)
pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
pl.title('Y-Axis Projections of Resolution, with subtraction', fontsize=titlesize)
#pl.xlim((0.0, 1000.0))
pl.ylim((0.0, pl.ylim()[1]+0.5))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl.savefig( write_file('plots/2dhistograms/seed%d_noise%d_signalMerged_resolution_PtOffline_projection_noPileup.png' % (args.seedEt_thresh, args.noise_filter)) )
pl.close()

