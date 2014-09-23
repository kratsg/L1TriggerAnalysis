from atlas_jets import *
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import os
import argparse
import time

#for getting FWHM
from scipy.interpolate import UnivariateSpline

#for getting error function fitting in diff/int curves
from scipy.optimize import curve_fit
from scipy.special import erf

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

def profile(xbins, xvals, yvals):
  # this finds the difference between successive bins and then halves it, and
  #     then adds it back to the bins to get xpoints
  xpoints = xbins[:-1] + np.array([v-xbins[i-1] for i,v in enumerate(xbins)][1:])/2.
  # this digitizes our samples by figuring out which xbin each xval belongs to
  digitized = np.digitize(xvals, xbins)
  ymean = [np.mean(yvals[np.where(digitized == i)]) for i in np.arange(1,len(xbins)+1)][1:]
  yerr  = [np.std(yvals[np.where(digitized == i)]) for i in np.arange(1,len(xbins)+1)][1:]
  return xpoints, ymean, yerr

def profile_y(xbins, xvals, yvals):
  return profile(xbins, xvals, yvals)

def profile_x(ybins, yvals, xvals):
  return profile(ybins, yvals, xvals)

filename_id = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed%d/leading_jets_%s.pkl" % (args.seedEt_thresh, filename_id)

data = pickle.load(file(filename))
numEvents = data.size

tJet_exists = [tJet.Et > 0 for tJet in data['matched_trigger_jet'] ]

gTower_Et = np.zeros((4,np.sum(tJet_exists)))

offline_jet_nsj, trigger_jet_Et, seed_Et, gTower_Et[0], gTower_Et[1], gTower_Et[2], gTower_Et[3] = np.array([ [oJet.nsj, tJet.Et, tJet.seed.Et, tJet.towers_around[0].Et, tJet.towers_around[1].Et, tJet.towers_around[2].Et, tJet.towers_around[3].Et]  for oJet,tJet in zip(data['leading_offline_jet'], data['matched_trigger_jet']) if tJet.Et > 0.0 ]).T

seed_Et = np.array([seed_Et]).reshape((1,-1))
gTower_Et = np.array(gTower_Et).reshape((4,-1))
all_gTower_Et = np.sort(np.vstack((seed_Et,gTower_Et)), axis=0)[::-1,:]

figsize = (16, 12)
labelsize = 28
titlesize = 36
linewidth = 4
light_grey = np.array([float(200)/float(255)]*3)
filled=False
markersize=160

dataSetStr  = 'TTbar 14TeV $\langle\mu\\rangle=80$'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

basicTextStr = '%s\n%d Events\n%s\n%s\n%s' % (dataSetStr, numEvents, seedCutStr, noiseCutStr, towerThrStr)

textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

startTime_wall      = time.time()
startTime_processor = time.clock()

#basic distribution plots
fig = pl.figure(figsize=figsize)

colors = ['b','r','c','m','y','k']
labels = ['$E_T^0$','$E_T^1$','$E_T^2$','$E_T^3$', '$E_T^4$']
markers = ['o','v','^','s','D','8']

pl.hist(seed_Et.T, bins=np.arange(0.0,100.0,2.5), stacked=True, fill=False, histtype='step', color='b', label='$E_T^\mathrm{seed}$', linewidth=linewidth)
for datapts, color, label in zip( all_gTower_Et, colors[1:], labels ):
  pl.hist(datapts, bins=np.arange(0.0,100.0,2.5), stacked=True, fill=False, histtype='step', color=color, label=label, linewidth=linewidth)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel('trigger jet\'s $E_T^\mathrm{gTower}$ [GeV]', fontsize=labelsize)
pl.ylabel(r'count/2.5 GeV bin', fontsize=labelsize)
pl.yscale('log')
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.text(0.05, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='left', bbox=textprops)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file("plots/gTowers/%s_Et.png" % (filename_id)) )
pl.close()

profiles = {}
profiles_scaled = {}

for datapts, label in zip(all_gTower_Et, ['gTower 0','gTower 1','gTower 2','gTower 3', 'gTower 4']):
    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(datapts, offline_jet_nsj, bins=(np.arange(0.,100.,2.5), np.arange(1,7,1)), norm=LogNorm())
    points_y, mean_x, err_x = profile_x(edges_y, offline_jet_nsj, datapts)
    profiles[label.replace(' ','_')] = {'x': mean_x, 'y': points_y, 'e': err_x}
    pl.scatter(mean_x, points_y, s=markersize, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'$E_T^\mathrm{%s}$' % label, fontsize=labelsize)
    pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
    pl.title('%s' % label, fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 100.0))
    pl.ylim((1.0, 6.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/gTowers/%s_%s.png' % (filename_id, label.replace(' ','_') )) )
    pl.close()

    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(datapts/trigger_jet_Et, offline_jet_nsj, bins=(np.arange(0.,1.,0.05), np.arange(1,7,1)), norm=LogNorm())
    points_y, mean_x, err_x = profile_x(edges_y, offline_jet_nsj, datapts/trigger_jet_Et)
    profiles_scaled[label.replace(' ','_')] = {'x': mean_x, 'y': points_y, 'e': err_x}
    pl.scatter(mean_x, points_y, s=markersize, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'$E_T^\mathrm{%s} / E_T^\mathrm{tJet}$' % label, fontsize=labelsize)
    pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
    pl.title('%s' % label, fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 1.0))
    pl.ylim((1.0, 6.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/gTowers/%s_%s_over_tJet_Et.png' % (filename_id, label.replace(' ','_') )) )
    pl.close()

#profile_x's merged
fig = pl.figure(figsize=figsize)
for l,c,m,yshift in zip(['gTower 0','gTower 1','gTower 2','gTower 3', 'gTower 4'], colors, markers, np.arange(-0.4,0.6,0.2)):
  datapts = profiles[l.replace(' ','_')]
  pl.errorbar(datapts['x'], datapts['y']+yshift, ms=20, mfc=c, mec='k', marker=m, linestyle='--', linewidth=3, capsize=10, alpha=0.75, label=l, xerr=datapts['e'], ecolor='k')
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize, numpoints=1)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
pl.xlabel(r'$E_T^\mathrm{gTower}$', fontsize=labelsize)
pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
pl.xlim((0.0, 100.0))
pl.ylim((1.0, 6.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/gTowers/%s_profiles.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
for l,c,m,yshift in zip(['gTower 0','gTower 1','gTower 2','gTower 3', 'gTower 4'], colors, markers, np.arange(-0.4,0.6,0.2)):
  datapts = profiles_scaled[l.replace(' ','_')]
  pl.errorbar(datapts['x'], datapts['y']+yshift, ms=20, mfc=c, mec='k', marker=m, linestyle='--', linewidth=3, capsize=10, alpha=0.75, label=l, xerr=datapts['e'], ecolor='k')
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize, numpoints=1)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
pl.xlabel(r'$E_T^\mathrm{gTower} / E_T^\mathrm{tJet}$', fontsize=labelsize)
pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
pl.xlim((0.0, 0.5))
pl.ylim((1.0, 6.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/gTowers/%s_profiles_over_tJet_Et.png' % (filename_id)) )
pl.close()


endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
