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

def distance(etas, phis):
  return ((etas[0] - etas[1])**2. + (phis[0] - phis[1])**2.)**0.5

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
 
def FWHM(bins, vals):
  spline = UnivariateSpline(bins[:-1]+np.diff(bins)/2., vals-np.max(vals)/2., s=0)
  roots = spline.roots() # find the roots
  r1, r2 = roots[0], roots[-1]
  return np.abs(r1-r2)

filename_id = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed%d/leading_jets_%s.pkl" % (args.seedEt_thresh, filename_id)

data = pickle.load(file(filename))

offline_jet_Pt, offline_jet_eta, offline_jet_phi, offline_jet_m = np.array([[oJet.Pt, oJet.eta, oJet.phi, oJet.m] for oJet in data['leading_offline_jet']]).T
trigger_jet_Et, trigger_jet_area, trigger_jet_region = np.array([[tJet.Et, tJet.area, tJet.region()] for tJet in data['matched_trigger_jet']]).T
weights = data['weight']

figsize = (16, 12)
labelsize = 28
titlesize = 36
linewidth = 4
light_grey = np.array([float(200)/float(255)]*3)
filled=False

dataSetStr  = 'Dijets 14 TeV, JZ0W.147910'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

basicTextStr = '%s\n%d Events\n%s\n%s\n%s' % (dataSetStr, data.size, seedCutStr, noiseCutStr, towerThrStr)

textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

startTime_wall      = time.time()
startTime_processor = time.clock()

#standard kinematic plots
fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_Pt, bins=100, density=False)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$p_T^\mathrm{oJet}$ [GeV]', fontsize=labelsize)
pl.ylabel('unweighted counts', fontsize=labelsize)
pl.yscale('log', nonposy='clip')
pl.text(0.95, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_pt_hist = {'Pt': offline_jet_Pt, 'weights': weights}
pickle.dump(pl_pt_hist, file( write_file('plots/pickle/%s_oJet_Pt_unweighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_Pt_unweighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_Pt, bins=100, weights=weights)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$p_T^\mathrm{oJet}$ [GeV]', fontsize=labelsize)
pl.ylabel('weighted counts', fontsize=labelsize)
pl.yscale('log', nonposy='clip')
pl.text(0.95, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_pt_hist = {'Pt': offline_jet_Pt, 'weights': weights}
pickle.dump(pl_pt_hist, file( write_file('plots/pickle/%s_oJet_Pt_weighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_Pt_weighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_eta, bins=100, density=False)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$\eta^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel('unweighted counts', fontsize=labelsize)
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_eta_hist = {'eta': offline_jet_eta, 'weights': weights}
pickle.dump(pl_eta_hist, file( write_file('plots/pickle/%s_oJet_eta_unweighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_eta_unweighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_eta, bins=100, weights=weights)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$\eta^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel('weighted counts', fontsize=labelsize)
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.rc('font', **{'size':labelsize})
fig.tight_layout()
pl_eta_hist = {'eta': offline_jet_eta, 'weights': weights}
pickle.dump(pl_eta_hist, file( write_file('plots/pickle/%s_oJet_eta_weighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_eta_weighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_phi, bins=100, density=False)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$\phi^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel('unweighted counts', fontsize=labelsize)
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_phi_hist = {'phi': offline_jet_phi, 'weights': weights}
pickle.dump(pl_phi_hist, file( write_file('plots/pickle/%s_oJet_phi_unweighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_phi_unweighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_phi, bins=100, weights=weights)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$\phi^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel('weighted counts', fontsize=labelsize)
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig.tight_layout()
pl_phi_hist = {'phi': offline_jet_phi, 'weights': weights}
pickle.dump(pl_phi_hist, file( write_file('plots/pickle/%s_oJet_phi_weighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_phi_weighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_m, bins=100, density=False)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$m^\mathrm{oJet}$ [GeV]', fontsize=labelsize)
pl.ylabel('unweighted counts', fontsize=labelsize)
pl.yscale('log', nonposy='clip')
pl.text(0.95, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_m_hist = {'m': offline_jet_m, 'weights': weights}
pickle.dump(pl_m_hist, file( write_file('plots/pickle/%s_oJet_m_unweighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_m_unweighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(offline_jet_m, bins=100, weights=weights)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$m^\mathrm{oJet}$ [GeV]', fontsize=labelsize)
pl.ylabel('weighted counts', fontsize=labelsize)
pl.yscale('log', nonposy='clip')
pl.text(0.95, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_m_hist = {'m': offline_jet_m, 'weights': weights}
pickle.dump(pl_m_hist, file( write_file('plots/pickle/%s_oJet_m_weighted.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_oJet_m_weighted.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
hist, bins = np.histogram(weights, bins=100, density=False)
pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, linewidth=linewidth)
pl.xlabel(r'$\omega_{event} * \sigma * \varepsilon_{filter} / N_{events}$', fontsize=labelsize)
pl.ylabel('counts', fontsize=labelsize)
pl.text(0.95, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.rc('font', **{'size':labelsize})
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl_weight_hist = {'weights': weights}
pickle.dump(pl_weight_hist, file( write_file('plots/pickle/%s_weights.pkl' % (filename_id)), 'w+') )
pl.savefig( write_file('plots/offline_jet_kinematics/%s_weights.png' % (filename_id)) )
pl.close()


fig = pl.figure(figsize=figsize)
counts, edges_x, edges_y, H = pl.hist2d(offline_jet_eta, offline_jet_phi, bins=(np.arange(-2.7,2.7,0.2), np.arange(-3.2,3.2,0.2)), norm=LogNorm(), alpha=0.75, weights=weights)
points_x, mean_y = profile_y(edges_x, offline_jet_eta, offline_jet_phi)
points_y, mean_x = profile_x(edges_y, offline_jet_phi, offline_jet_eta)
cbar = pl.colorbar()
pl.scatter(points_x, mean_y, s=80, facecolor='none', edgecolor='k', marker='o', linewidth=2, label='profile y')
pl.scatter(mean_x, points_y, s=80, facecolor='none', edgecolor='k', marker='x', linewidth=2, label='profile x')
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$\eta^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel(r'$\phi^\mathrm{oJet}$', fontsize=labelsize)
cbar.set_label('number density', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/offline_jet_kinematics/%s_angular_positions.png' % (filename_id)) )
pl.close()

#y projection slices
fig = pl.figure(figsize=figsize)
selections = {'left':   np.where((-1.9 < offline_jet_eta)&(offline_jet_eta < -1.5)),\
              'middle': np.where((-0.2 < offline_jet_eta)&(offline_jet_eta < 0.2)),\
              'right':  np.where((1.5 < offline_jet_eta)&(offline_jet_eta < 1.9))}
pl.hist(offline_jet_phi[selections['left']], bins=100, weights=weights[selections['left']], stacked=True, fill=False, histtype='step', linewidth=linewidth, alpha=0.75, color='b', label='$-1.9\ < \eta^\mathrm{oJet} <\ -1.5$')
pl.hist(offline_jet_phi[selections['middle']], bins=100, weights=weights[selections['middle']], stacked=True, fill=False, histtype='step', linewidth=linewidth, alpha=0.75, color='c', label='$-0.2\ < \eta^\mathrm{oJet} <\ 0.2$')
pl.hist(offline_jet_phi[selections['right']], bins=100, weights=weights[selections['right']], stacked=True, fill=False, histtype='step', linewidth=linewidth, alpha=0.75, color='r', label='$1.5\ < \eta^\mathrm{oJet} <\ 1.9$')
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$\eta^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel('weighted counts', fontsize=labelsize)
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
pl.title('Y-Axis Projections of $\phi^\mathrm{oJet}$', fontsize=titlesize)
pl.yscale('log', nonposy='clip')
#pl.xlim((-1.0,1.0))
#pl.ylim((0.0, 1000.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl.savefig( write_file('plots/offline_jet_kinematics/%s_offline_jet_phi_projection_y.png' % (filename_id)) )
pl.close()


fig = pl.figure(figsize=figsize)
counts, edges_x, edges_y, H = pl.hist2d(offline_jet_eta, offline_jet_Pt, bins=(np.arange(-2.7,2.9,0.2), np.arange(0.,375.,10.)), norm=LogNorm(), alpha=0.75, weights=weights)
points_x, mean_y = profile_y(edges_x, offline_jet_eta, offline_jet_Pt)
points_y, mean_x = profile_x(edges_y, offline_jet_Pt, offline_jet_eta)
cbar = pl.colorbar()
pl.scatter(points_x, mean_y, s=80, facecolor='none', edgecolor='k', marker='o', linewidth=2, label='profile y')
pl.scatter(mean_x, points_y, s=80, facecolor='none', edgecolor='k', marker='x', linewidth=2, label='profile x')
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$\eta^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel(r'$p_T^\mathrm{oJet}$', fontsize=labelsize)
cbar.set_label('number density', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/offline_jet_kinematics/%s_eta_pt_corr.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt, offline_jet_phi, bins=(np.arange(0.,350.,10.), np.arange(-3.2,3.2,0.2)), norm=LogNorm(), alpha=0.75, weights=weights)
points_x, mean_y = profile_y(edges_x, offline_jet_Pt, offline_jet_phi)
points_y, mean_x = profile_x(edges_y, offline_jet_phi, offline_jet_Pt)
cbar = pl.colorbar()
pl.scatter(points_x, mean_y, s=80, facecolor='none', edgecolor='k', marker='o', linewidth=2, label='profile y')
pl.scatter(mean_x, points_y, s=80, facecolor='none', edgecolor='k', marker='x', linewidth=2, label='profile x')
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$p_T^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel(r'$\phi^\mathrm{oJet}$', fontsize=labelsize)
cbar.set_label('number density', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/offline_jet_kinematics/%s_pt_phi_corr.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
counts, edges_x, edges_y, H = pl.hist2d(offline_jet_eta, offline_jet_m, bins=(np.arange(-2.7,2.9,0.2), np.arange(0.,175.,10.)), norm=LogNorm(), alpha=0.75, weights=weights)
points_x, mean_y = profile_y(edges_x, offline_jet_eta, offline_jet_m)
points_y, mean_x = profile_x(edges_y, offline_jet_m, offline_jet_eta)
cbar = pl.colorbar()
pl.scatter(points_x, mean_y, s=80, facecolor='none', edgecolor='k', marker='o', linewidth=2, label='profile y')
pl.scatter(mean_x, points_y, s=80, facecolor='none', edgecolor='k', marker='x', linewidth=2, label='profile x')
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$\eta^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel(r'$m^\mathrm{oJet}$', fontsize=labelsize)
cbar.set_label('number density', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/offline_jet_kinematics/%s_eta_m_eta.png' % (filename_id)) )
pl.close()

fig = pl.figure(figsize=figsize)
counts, edges_x, edges_y, H = pl.hist2d(offline_jet_m, offline_jet_phi, bins=(np.arange(0.,175.,10.), np.arange(-3.2,3.2,0.2)), norm=LogNorm(), alpha=0.75, weights=weights)
points_x, mean_y = profile_y(edges_x, offline_jet_m, offline_jet_phi)
points_y, mean_x = profile_x(edges_y, offline_jet_phi, offline_jet_m)
cbar = pl.colorbar()
pl.scatter(points_x, mean_y, s=80, facecolor='none', edgecolor='k', marker='o', linewidth=2, label='profile y')
pl.scatter(mean_x, points_y, s=80, facecolor='none', edgecolor='k', marker='x', linewidth=2, label='profile x')
pl.text(0.05, 0.95, basicTextStr, transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$m^\mathrm{oJet}$', fontsize=labelsize)
pl.ylabel(r'$\phi^\mathrm{oJet}$', fontsize=labelsize)
cbar.set_label('number density', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/offline_jet_kinematics/%s_m_phi_corr.png' % (filename_id)) )
pl.close()


endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

