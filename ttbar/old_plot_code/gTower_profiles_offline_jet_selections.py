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
  return xpoints, ymean

def profile_y(xbins, xvals, yvals):
  return profile(xbins, xvals, yvals)

def profile_x(ybins, yvals, xvals):
  return profile(ybins, yvals, xvals)
 
filename_id = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed%d/leading_jets_%s.pkl" % (args.seedEt_thresh, filename_id)

data = pickle.load(file(filename))
numEvents = data.size

data_dtype = [\
              ('oJet_Pt','float64'),\
              ('oJet_m','float64'),\
              ('oJet_n','uint8'),\
              ('oJet_subjetsPt','object'),\
              ('tJet_Et','float64'),\
              ('gTowers_Et','5float64')\
            ]

data = np.array([ (oJet.Pt, oJet.m, oJet.nsj, oJet.subjetsPt, tJet.Et, sorted([tJet.seed.Et] + [gTower.Et for gTower in tJet.towers_around]))  for oJet,tJet in zip(data['leading_offline_jet'], data['matched_trigger_jet']) if tJet.Et > 0.0 ], dtype=data_dtype)

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

def custom_filter(listOfObjects, cut, numLeft):
  ''' this can be called in two ways
      custom_filter(subjets, subjetPt_cut, nsj_cut)
      custom_filter(gTowers/tJet.Et, f_cut, gTower_N_cut)
  '''
  return (np.sum(np.array(listOfObjects) > cut) >= numLeft)

#basic distribution plots
fig = pl.figure(figsize=figsize)

colors = ['b','r','c','m','y','k']
labels = ['$E_T^\mathrm{seed}$','$E_T^1$','$E_T^2$','$E_T^3$', '$E_T^4$']
markers = ['o','v','^','s','D','8']

yvals = data['oJet_n'] #doesnt' change
tJet_Et = data['tJet_Et'] #doesn't change

num_gTowers = data['gTowers_Et'].shape[1]

profiles = [{}]*num_gTowers
profiles_scaled = [{}]*num_gTowers

bins_gTower_Et = np.arange(0.0,100.0,2.5)

for i in range(num_gTowers):
    xvals = data['gTowers_Et'][:,i]

    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(xvals, yvals, bins=(bins_gTower_Et, np.arange(1,7,1)), norm=LogNorm())
    points_y, mean_x = profile_x(edges_y, yvals, xvals)
    profiles[i] = {'x': mean_x, 'y': points_y}
    pl.scatter(mean_x, points_y, s=markersize, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'$E_T^\mathrm{gTower\ %s}$' % (i+1), fontsize=labelsize)
    pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
    pl.title('gTower %s' % (i+1), fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 100.0))
    pl.ylim((1.0, 6.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/gTowers/%s_gTower_%s.png' % (filename_id, (i+1) )) )
    pl.close()

    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(xvals/tJet_Et, yvals, bins=(np.arange(0.,1.,0.05), np.arange(1,7,1)), norm=LogNorm())
    points_y, mean_x = profile_x(edges_y, yvals, xvals/tJet_Et)
    profiles_scaled[i] = {'x': mean_x, 'y': points_y}
    pl.scatter(mean_x, points_y, s=markersize, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'$E_T^\mathrm{gTower\ %s} / E_T^\mathrm{tJet}$' % (i+1), fontsize=labelsize)
    pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
    pl.title('gTower %s' % (i+1), fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 1.0))
    pl.ylim((1.0, 6.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/gTowers/%s_gTower_%s_over_tJet_Et.png' % (filename_id, (i+1) )) )
    pl.close()

#profile_x's merged
fig = pl.figure(figsize=figsize)
for i in range(num_gTowers):
  c = colors[i]
  m = markers[i]
  datapts = profiles[i]
  pl.scatter(datapts['x'], datapts['y'], s=markersize, facecolor=c, edgecolor='k', marker=m, linewidth=1, alpha=0.75, label='gTower %s' % i)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize, scatterpoints=1)
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
for i in range(num_gTowers):
  c = colors[i]
  m = markers[i]
  datapts = profiles_scaled[i]
  pl.scatter(datapts['x'], datapts['y'], s=markersize, facecolor=c, edgecolor='k', marker=m, linewidth=1, alpha=0.75, label='gTower %s' % i)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize, scatterpoints=1)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
pl.xlabel(r'$E_T^\mathrm{gTower} / E_T^\mathrm{tJet}$', fontsize=labelsize)
pl.ylabel('offline jet\'s number of subjets (cut at 6)', fontsize=labelsize)
pl.xlim((0.0, 0.4))
pl.ylim((1.0, 6.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/gTowers/%s_profiles_over_tJet_Et.png' % (filename_id)) )
pl.close()


def custom_filter(listOfObjects, cut, numLeft):
  ''' this can be called in two ways
      custom_filter(subjets, subjetPt_cut, nsj_cut)
      custom_filter(gTowers/tJet.Et, f_cut, gTower_N_cut)
  '''
  if numLeft == 0:
    return True
  return (np.sum(np.array(listOfObjects) > cut) >= numLeft)

# building the offline jet selections

# subjet done with custom_filter()
offline_subjet_selection = [0, 1, 2, 3]
offline_subjetPt_selection = [20.,20.,20.,20.]
offline_jet_Pt_cut = (200 <= data['oJet_Pt'])&(data['oJet_Pt'] <= 300.)

# this just store the array of motherfucking truth values
offline_jet_m_selection = [data['oJet_m'] > 0.0, data['oJet_m'] > 0.0, (50 <= data['oJet_m'])&(data['oJet_m'] <= 110.), (100. <= data['oJet_m'])&(data['oJet_m'] <= 220.)]

for oJet_nsj_req, oJet_subjetPt_req, oJet_m_cut in zip(offline_subjet_selection, offline_subjetPt_selection, offline_jet_m_selection):
  #build the offline subjet selection
  oJet_subjet_cut = np.array(map(lambda x: custom_filter(x, oJet_subjetPt_req, oJet_nsj_req), data['oJet_subjetsPt']))

  #array[rows][col][sub-set if applicable]
  overall_cut = np.where( oJet_subjet_cut & oJet_m_cut & offline_jet_Pt_cut )

  tJet_Et = data[overall_cut]['tJet_Et']

  basicTextStr = '%s\n%d Events\n%s\n%s\n%s' % (dataSetStr, tJet_Et.size, seedCutStr, noiseCutStr, towerThrStr)

  fig = pl.figure(figsize=figsize)
  for i in range(num_gTowers):
    datapts = data[overall_cut]['gTowers_Et'][:,i]
    pl.hist(datapts, bins=bins_gTower_Et, color=colors[i], linewidth=3, alpha=0.75, histtype='step', fill=False, label='gTower %s' % i)
  legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
  legend.get_frame().set_facecolor(light_grey)
  legend.get_frame().set_linewidth(0.0)
  pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
  pl.xlabel(r'$E_T^\mathrm{gTower}$ [GeV]', fontsize=labelsize)
  pl.ylabel('count', fontsize=labelsize)
  pl.title('Leading oJet, number of subjets = %s' % oJet_nsj_req, fontsize=titlesize)
  pl.xlim((0.0, 100.0))
  pl.yscale('log', nonposy='clip')
  pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
  pl.tick_params(axis='both', which='major', labelsize=labelsize)
  pl.savefig( write_file('plots/gTowers_oJet_selection/%s_profiles_nsj%s.png' % (filename_id, oJet_nsj_req)) )
  pl.close()

  fig = pl.figure(figsize=figsize)
  for i in range(num_gTowers):
    datapts = data[overall_cut]['gTowers_Et'][:,i] / tJet_Et
    pl.hist(datapts, bins=np.arange(0.,1.,0.05), color=colors[i], linewidth=3, alpha=0.75, histtype='step', fill=False, label='gTower %s' % i)
  legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
  legend.get_frame().set_facecolor(light_grey)
  legend.get_frame().set_linewidth(0.0)
  pl.text(0.95, 0.05, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
  pl.xlabel(r'$E_T^\mathrm{gTower} / E_T^\mathrm{tJet}$', fontsize=labelsize)
  pl.ylabel('count', fontsize=labelsize)
  pl.title('Leading oJet, number of subjets = %s' % oJet_nsj_req, fontsize=titlesize)
  pl.xlim((0.0, 1.0))
  pl.yscale('log', nonposy='clip')
  pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
  pl.tick_params(axis='both', which='major', labelsize=labelsize)
  pl.savefig( write_file('plots/gTowers_oJet_selection/%s_profiles_over_tJet_Et_nsj%s.png' % (filename_id, oJet_nsj_req)) )
  pl.close()


endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
