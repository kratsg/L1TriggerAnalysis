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

def binomial_errors(hist_ratio, hist_one, hist_two):
  errors = []
  for w, num, den in zip(hist_ratio, hist_one, hist_two):
    # root.cern.ch/root/html/src/TH1.cxx.html#l5.yxD
    # formula cited (for histograms [num, den] with no errors) is:
    #     w = num/den
    #     if w = 1:
    #             sigma = 0
    #     else:
    #             sigma = abs( (1 - 2*w + w**2) / den**2 )
    if w == 1.0:
      errors.append(0.0)
    else:
      errors.append( (np.abs( (1.-2.*w + w**2.)/den**2.))**0.5/2. )
  return errors

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

offline_jet_Pt = {}
trigger_jet_Et = {}
trigger_jet_area = {}

offline_jet_Pt['all'] = np.array([oJet.Pt for oJet in data['leading_offline_jet']])
trigger_jet_Et['all'], trigger_jet_area['all'], trigger_jet_region, trigger_jet_eta = np.array([[tJet.Et, tJet.area, tJet.region(), tJet.eta] for tJet in data['matched_trigger_jet']]).T

#we added subregions to the forward region
sub_region_a = (1.6 < np.fabs(trigger_jet_eta))&(np.fabs(trigger_jet_eta) < 2.5)
sub_region_b = (2.5 < np.fabs(trigger_jet_eta))&(np.fabs(trigger_jet_eta) < 3.2)
sub_region_c = (3.2 < np.fabs(trigger_jet_eta))&(np.fabs(trigger_jet_eta) < 4.9)

sub_region = {}
sub_region['3a'] = np.where( sub_region_a&(trigger_jet_region == 3) )
sub_region['3b'] = np.where( sub_region_b&(trigger_jet_region == 3) )
sub_region['3c'] = np.where( sub_region_c&(trigger_jet_region == 3) )
sub_region['4a'] = np.where( sub_region_a&(trigger_jet_region == 4) )
sub_region['4b'] = np.where( sub_region_b&(trigger_jet_region == 4) )
sub_region['4c'] = np.where( sub_region_c&(trigger_jet_region == 4) )

trigger_jet_Et_correction = {}
trigger_jet_Et_noPileup = {}
numEvents = {}
trigger_jet_exists = {}
trigger_jet_exists_noPileup = {}
# the new regions aren't working because there are not enough jets in them!!!
regions = [1,2,3,4,'3a','3b','3c','4a','4b','4c']

for i in regions:
  try:
    int(i)
    region = np.where(trigger_jet_region == i)
  except:
    region = sub_region[i]

  regionStr = 'region_%s' % i
  offline_jet_Pt[regionStr] = offline_jet_Pt['all'][region]
  trigger_jet_Et[regionStr] = trigger_jet_Et['all'][region]
  trigger_jet_area[regionStr] = trigger_jet_area['all'][region]

  #handle 3x and 4x
  try:
    ii = int(i)
  except:
    ii = int(i[0])

  trigger_jet_Et_correction[regionStr] = trigger_jet_area[regionStr]*data['gFEX_rho_%d' % ii][region]
  trigger_jet_Et_noPileup[regionStr] = trigger_jet_Et[regionStr] - trigger_jet_Et_correction[regionStr]
  numEvents[regionStr] = offline_jet_Pt[regionStr].size
  trigger_jet_exists[regionStr] = np.where(trigger_jet_Et[regionStr] > 0.)
  trigger_jet_exists_noPileup[regionStr] = np.where(trigger_jet_Et_noPileup[regionStr] > 0.)

del trigger_jet_region
#del offline_jet_Pt['all']
del trigger_jet_Et['all']
del trigger_jet_area['all']

bins_efficiency = np.arange(0.,2000.,10.)
width_efficiency = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])
bins_multiplicity = np.arange(0.0,50.0,1.0)

bins_rho = np.arange(0.,100.,1.)
bins_vertices = np.arange(0.,100.,1.)

figsize = (16, 12)
labelsize = 28
titlesize = 36
linewidth = 4
light_grey = np.array([float(200)/float(255)]*3)
filled=False

dataSetStr  = 'TTbar 14TeV $\langle\mu\\rangle=80$'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

basicTextStr = {}
for i in regions:
  basicTextStr['region_%s' % i] = '%s\n%d Events \n%s\n%s\n%s' % (dataSetStr, numEvents['region_%s' % i], seedCutStr, noiseCutStr, towerThrStr)

basicTextStr['default'] = '%s\n%d Events\n%s\n%s\n%s' % (dataSetStr, (numEvents['region_1'] + numEvents['region_2'] + numEvents['region_3'] + numEvents['region_4']), seedCutStr, noiseCutStr, towerThrStr)

textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

label_oRho = r'Offline $\rho$ [GeV]'
label_gRho = r'gFEX $\rho$ [GeV]'
label_nVer = r'Number of primary vertices (vxp_nTracks $\geq$ 2)'

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

startTime_wall      = time.time()
startTime_processor = time.clock()

#multiplicity on gTowers
fig = pl.figure(figsize=figsize)
# b r c m y k
where = np.where(offline_jet_Pt['all'] > 0.)
pl.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size, linestyle='steps-post', alpha=0.75, color='b', label='no cuts\n%d events' % where[0].size, linewidth=linewidth)
where = np.where((offline_jet_Pt['all'] > 100.)&(offline_jet_Pt['all'] < 150.))
pl.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size, linestyle='steps-post', alpha=0.75, color='r', label='$100 < p_T^\mathrm{oJet} < 150$\n%d events' % where[0].size, linewidth=linewidth)
where = np.where((offline_jet_Pt['all'] > 150.)&(offline_jet_Pt['all'] < 200.))
pl.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size, linestyle='steps-post', alpha=0.75, color='c', label='$150 < p_T^\mathrm{oJet} < 200$\n%d events' % where[0].size, linewidth=linewidth)
where = np.where((offline_jet_Pt['all'] > 200.)&(offline_jet_Pt['all'] < 250.))
pl.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size, linestyle='steps-post', alpha=0.75, color='m', label='$200 < p_T^\mathrm{oJet} < 250$\n%d events' % where[0].size, linewidth=linewidth)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$E_T^\mathrm{gTower}$ [GeV]', fontsize=labelsize)
pl.ylabel('gTower multiplicity / event', fontsize=labelsize)
pl.yscale('log', nonposy='clip')
pl.text(0.05, 0.95, basicTextStr['default'], transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='left', bbox=textprops)
pl.ylim((0.0, 1284.0))
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
fig.tight_layout()
pl.savefig( write_file('plots/multiplicity/%s.png' % (filename_id)) )
pl.close()

#for col,legend in zip(['gFEX_rho_all','gFEX_rho_1','gFEX_rho_2','gFEX_rho_3','gFEX_rho_4'],['all',1,2,3,4]):
for i in regions:
  region = 'region_%s' % i

  try:
    #make a correlation of corrected Et versus trigger jet Et
    corr = np.corrcoef(trigger_jet_Et[region], trigger_jet_Et_correction[region])[0,1]
    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(trigger_jet_Et[region], trigger_jet_Et_correction[region], bins=np.arange(0.,1500.,10.), norm=LogNorm(), alpha=0.75)
    points_x, mean_y = profile_y(edges_x, trigger_jet_Et[region], trigger_jet_Et_correction[region])
    cbar = pl.colorbar()
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    pl.text(0.95, 0.95, '%s\nCorr = %0.4f' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'$E_T^{\mathrm{gJet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel(r'$\rho*A^\mathrm{gJet}$ [GeV]', fontsize=labelsize)
    pl.xlim((0.,1000.))
    pl.ylim((0.,200.))
    cbar.set_label('number density', fontsize=labelsize)
    #pl.xlim((0.0, 1000.0))
    #pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/trigger_jet_Et_correction/%s_trigger_jet_Et_correction_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not make correlation of corrected Et versus trigger jet Et" % region
    pass
  
  try:
    corr = np.corrcoef(offline_jet_Pt[region][trigger_jet_exists[region]], trigger_jet_Et[region][trigger_jet_exists[region]])[0,1]
    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt[region][trigger_jet_exists[region]], trigger_jet_Et[region][trigger_jet_exists[region]], bins=np.arange(0.,1500.,10.), norm=LogNorm())
    points_x, mean_y = profile_y(edges_x, offline_jet_Pt[region][trigger_jet_exists[region]], trigger_jet_Et[region][trigger_jet_exists[region]])
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.05, '%s\nCorr = %0.4f' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel(r'trigger $E_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.title('no subtraction', fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 1000.0))
    pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/jet_energy_correlation/%s_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not make jet energy correlation" % region
    pass

  try:
    corr = np.corrcoef(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], trigger_jet_Et_noPileup[region][trigger_jet_exists_noPileup[region]])[0,1]
    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], trigger_jet_Et_noPileup[region][trigger_jet_exists_noPileup[region]], bins=np.arange(0.,1500.,10.), norm=LogNorm())
    points_x, mean_y = profile_y(edges_x, offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], trigger_jet_Et_noPileup[region][trigger_jet_exists_noPileup[region]])
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.05, '%s\nCorr = %0.4f' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel(r'trigger $E_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.title('with subtraction', fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 1000.0))
    pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/jet_energy_correlation/%s_noPileup_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not make jet energy correlation (no Pileup)" % region
    pass

  #we want to make resolution plots now
  #   note that we ignore trigger jets that are negative or less than zero
  resolution = (trigger_jet_Et[region] - offline_jet_Pt[region])/offline_jet_Pt[region]
  resolution_noPileup = (trigger_jet_Et_noPileup[region] - offline_jet_Pt[region])/offline_jet_Pt[region]
  #eta = np.array([oJet.eta for oJet in data['leading_offline_jet']])
  #rho = data['offline_rho']
  #vxp_n = data['vxp_n']

  try:
    corr = np.corrcoef(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution[trigger_jet_exists_noPileup[region]])[0,1]
    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution[trigger_jet_exists_noPileup[region]], bins=100, norm=LogNorm())
    points_x, mean_y = profile_y(edges_x, offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution[trigger_jet_exists_noPileup[region]])
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.95, '%s\nCorr = %0.4f' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel(r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', fontsize=labelsize)
    pl.title('no subtraction', fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 1000.0))
    #pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/resolution/%s_resolution_PtOffline_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not make resolution correlation" % region
    pass

  try:
    corr = np.corrcoef(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution_noPileup[trigger_jet_exists_noPileup[region]])[0,1]
    fig = pl.figure(figsize=figsize)
    counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution_noPileup[trigger_jet_exists_noPileup[region]], bins=100, norm=LogNorm())
    points_x, mean_y = profile_y(edges_x, offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution_noPileup[trigger_jet_exists_noPileup[region]])
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.95, '%s\nCorr = %0.4f' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
    pl.xlabel(r'offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel(r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', fontsize=labelsize)
    pl.title('with subtraction', fontsize=titlesize)
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlim((0.0, 1000.0))
    #pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/resolution/%s_resolution_PtOffline_noPileup_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not make resolution correlation (no Pileup)" % region
    pass

  try:
    #y projection slices
    fig = pl.figure(figsize=figsize)
    hist, bins = np.histogram(resolution[np.where((170. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 180.)&(trigger_jet_Et_noPileup[region] > 0.))], bins=100, density=True)
    fwhm = FWHM(bins, hist)
    pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='b', label='$170\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 180\ \mathrm{GeV}$\nFWHM = %0.4f' % fwhm, linewidth=linewidth)
    hist, bins = np.histogram(resolution[np.where((200. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 220.)&(trigger_jet_Et_noPileup[region] > 0.))], bins=100, density=True)
    fwhm = FWHM(bins, hist)
    pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='r', label='$200\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 220\ \mathrm{GeV}$\nFWHM = %0.4f' % fwhm, linewidth=linewidth)
    hist, bins = np.histogram(resolution[np.where((300. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 350.)&(trigger_jet_Et_noPileup[region] > 0.))], bins=100, density=True)
    fwhm = FWHM(bins, hist)
    pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='c', label='$300\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 350\ \mathrm{GeV}$\nFWHM = %0.4f' % fwhm, linewidth=linewidth)
    legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
    legend.get_frame().set_facecolor(light_grey)
    legend.get_frame().set_linewidth(0.0)
    pl.xlabel(r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', fontsize=labelsize)
    pl.ylabel('normalized counts', fontsize=labelsize)
    pl.text(0.95, 0.05, basicTextStr[region], transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.title('Y-Axis Projections of Resolution', fontsize=titlesize)
    pl.xlim((-1.0,1.0))
    #pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    fig.tight_layout()
    pl_res_proj = {'170to180': resolution[np.where((170. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 180.)&(trigger_jet_Et_noPileup[region] > 0.))],\
                   '200to220': resolution[np.where((200. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 220.)&(trigger_jet_Et_noPileup[region] > 0.))],\
                   '300to350': resolution[np.where((300. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 350.)&(trigger_jet_Et_noPileup[region] > 0.))]}
    pickle.dump(pl_res_proj, file( write_file('plots/pickle/%s_resolution_PtOffline_projection_region%s.pkl' % (filename_id, i)), 'w+') )
    pl.savefig( write_file('plots/resolution/%s_resolution_PtOffline_projection_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not resolution projection" % region
    pass

  try:
    fig = pl.figure(figsize=figsize)
    hist, bins = np.histogram(resolution_noPileup[np.where((170. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 180.)&(trigger_jet_Et_noPileup[region] > 0.))], bins=100, density=True)
    fwhm = FWHM(bins, hist)
    pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='b', label='$170\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 180\ \mathrm{GeV}$\nFWHM = %0.4f' % fwhm, linewidth=linewidth)
    hist, bins = np.histogram(resolution_noPileup[np.where((200. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 220.)&(trigger_jet_Et_noPileup[region] > 0.))], bins=100, density=True)
    fwhm = FWHM(bins, hist)
    pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='r', label='$200\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 220\ \mathrm{GeV}$\nFWHM = %0.4f' % fwhm, linewidth=linewidth)
    hist, bins = np.histogram(resolution_noPileup[np.where((300. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 350.)&(trigger_jet_Et_noPileup[region] > 0.))], bins=100, density=True)
    fwhm = FWHM(bins, hist)
    pl.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='c', label='$300\ \mathrm{GeV} < p_T^\mathrm{oJet} <\ 350\ \mathrm{GeV}$\nFWHM = %0.4f' % fwhm, linewidth=linewidth)
    legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
    legend.get_frame().set_facecolor(light_grey)
    legend.get_frame().set_linewidth(0.0)
    pl.xlabel(r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', fontsize=labelsize)
    pl.ylabel('normalized counts', fontsize=labelsize)
    pl.text(0.95, 0.05, basicTextStr[region], transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.title('Y-Axis Projections of Resolution, with subtraction', fontsize=titlesize)
    pl.xlim((-1.0,1.0))
    #pl.ylim((0.0, 1000.0))
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    fig.tight_layout()
    pl_res_proj = {'170to180': resolution_noPileup[np.where((170. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 180.)&(trigger_jet_Et_noPileup[region] > 0.))],\
                   '200to220': resolution_noPileup[np.where((200. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 220.)&(trigger_jet_Et_noPileup[region] > 0.))],\
                   '300to350': resolution_noPileup[np.where((300. < offline_jet_Pt[region])&(offline_jet_Pt[region] < 350.)&(trigger_jet_Et_noPileup[region] > 0.))]}
    pickle.dump(pl_res_proj, file( write_file('plots/pickle/%s_resolution_PtOffline_projection_noPileup_region%s.pkl' % (filename_id, i)), 'w+') )
    pl.savefig( write_file('plots/resolution/%s_resolution_PtOffline_projection_noPileup_region%s.png' % (filename_id, i)) )
    pl.close()
  except:
    print "Error for %s: could not make resolution projection (no Pileup)" % region
    pass

  for triggerEt_thresh in [0., 100., 110., 130., 150., 200., 220.]:
    '''ADD IN np.interp(triggerEt_thresh, trigger_Et, rho*A) and shift subtracted eff curves'''
    eff_curve_shift = np.interp(triggerEt_thresh, trigger_jet_Et[region], trigger_jet_Et_correction[region])

    textstr = '%s\n$E_T^\mathrm{gFEX\ jet} >\ %d\ \mathrm{GeV}$' % (basicTextStr[region], triggerEt_thresh)

    #*_noPileup means the data with the pileup subtracted

    hist_efficiency_den, _          = np.histogram(offline_jet_Pt[region], bins=bins_efficiency)
    hist_efficiency_num, _          = np.histogram(offline_jet_Pt[region][np.where(trigger_jet_Et[region] > triggerEt_thresh + eff_curve_shift)], bins=bins_efficiency)
    hist_efficiency_num_noShift, _  = np.histogram(offline_jet_Pt[region][np.where(trigger_jet_Et[region] > triggerEt_thresh)], bins=bins_efficiency)
    hist_efficiency_num_noPileup, _ = np.histogram(offline_jet_Pt[region][np.where(trigger_jet_Et_noPileup[region] > triggerEt_thresh)], bins=bins_efficiency)

    nonzero_bins = np.where(hist_efficiency_den != 0)
    #compute integral and differential curves
    hist_eff_curve_diff          = np.true_divide(hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    hist_eff_curve_diff_noShift  = np.true_divide(hist_efficiency_num_noShift[nonzero_bins], hist_efficiency_den[nonzero_bins])
    hist_eff_curve_diff_noPileup = np.true_divide(hist_efficiency_num_noPileup[nonzero_bins], hist_efficiency_den[nonzero_bins])

    hist_eff_curve_int          = np.true_divide(np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
    hist_eff_curve_int_noPileup = np.true_divide(np.cumsum(hist_efficiency_num_noPileup[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])

    #get halfway in between really
    xpoints_efficiency = bins_efficiency[:-1] + width_efficiency/2.

    #binomial errors s^2 = n * p * q
    errors_eff_diff          = binomial_errors(hist_eff_curve_diff, hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    errors_eff_diff_noShift  = binomial_errors(hist_eff_curve_diff_noShift, hist_efficiency_num_noShift[nonzero_bins], hist_efficiency_den[nonzero_bins])
    errors_eff_diff_noPileup = binomial_errors(hist_eff_curve_diff_noPileup, hist_efficiency_num_noPileup[nonzero_bins], hist_efficiency_den[nonzero_bins])

    errors_eff_int          = binomial_errors(hist_eff_curve_int, np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
    errors_eff_int_noPileup = binomial_errors(hist_eff_curve_int_noPileup, np.cumsum(hist_efficiency_num_noPileup[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])

    xlim_efficiency = (0.0, 450.0) #GeV
    ylim_efficiency = (0.0,1.1)

    #define erfx used for error fitting
    def func(x, a, b, c, d):
      # note that b == sigma here, see wiki for more info
      return a*erf( (x-c)/b ) + d

    def fit_func(x, y):
      try:
        popt, pcov = curve_fit(func, x, y, p0=(1., 50., triggerEt_thresh, 0.5))
      except RuntimeError:
        return -1.
      return popt[1] #return the width (b)

    fig = pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff, yerr=errors_eff_diff, ecolor='black', label='no subtraction\nWidth = %0.4f' % w, linewidth=linewidth)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noShift)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noShift, yerr=errors_eff_diff_noShift, ecolor='black', label='no sub. and shift\nWidth = %0.4f' % w, linewidth=linewidth)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noPileup)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noPileup, yerr=errors_eff_diff_noPileup, ecolor='black', label='with subtraction\nWidth = %0.4f' % w, linewidth=linewidth)
    pl.text(0.95, 0.05, '%s\nShift: %0.4f' % (textstr, eff_curve_shift), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlim(xlim_efficiency)
    pl.ylim(ylim_efficiency)
    legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
    legend.get_frame().set_facecolor(light_grey)
    legend.get_frame().set_linewidth(0.0)
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl_eff_diff = {'xdata'          : xpoints_efficiency,\
                   'ydata'          : hist_eff_curve_diff,\
                   'ydata_noPileup' : hist_eff_curve_diff_noPileup,\
                   'ydata_noShift'  : hist_eff_curve_diff_noShift,\
                   'xerr'           : 1.0,\
                   'yerr'           : errors_eff_diff,\
                   'yerr_noPileup'  : errors_eff_diff_noPileup,\
                   'yerr_noShift'   : errors_eff_diff_noShift,\
                   'num'            : hist_efficiency_num,\
                   'num_noPileup'   : hist_efficiency_num_noPileup,\
                   'num_noShift'    : hist_efficiency_num_noShift,\
                   'den'            : hist_efficiency_den,\
                   'bins'           : bins_efficiency,\
                   'nonzero_bins'   : nonzero_bins,\
                   'eff_curve_shift': eff_curve_shift}
    pickle.dump(pl_eff_diff, file( write_file('plots/pickle/matched_differential_%s_trigger%d_region%s.pkl' % (filename_id, triggerEt_thresh, i)), 'w+') )
    pl.savefig( write_file('plots/differential/%s_trigger%d_region%s.png' % (filename_id, triggerEt_thresh, i)) )
    pl.close()

    fig = pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Integral', fontsize=labelsize)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_int)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_int, yerr=errors_eff_int, ecolor='black', label='no subtraction\nWidth = %0.4f' % w, linewidth=linewidth)
    w = fit_func(xpoints_efficiency[nonzero_bins]-eff_curve_shift, hist_eff_curve_int_noPileup)
    pl.errorbar(xpoints_efficiency[nonzero_bins]-eff_curve_shift, hist_eff_curve_int_noPileup, yerr=errors_eff_int_noPileup, ecolor='black', label='with subtraction\nWidth = %0.4f' % w, linewidth=linewidth)
    pl.text(0.95, 0.05, '%s\nShift: %0.4f' % (textstr, eff_curve_shift), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
    pl.xlim(xlim_efficiency)
    pl.ylim(ylim_efficiency)
    pl.legend()
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl_eff_int = {'xdata'          : xpoints_efficiency,\
                  'ydata'          : hist_eff_curve_int,\
                  'ydata_noPileup' : hist_eff_curve_int_noPileup,\
                  'xerr'           : 1.0,\
                  'yerr'           : errors_eff_int,\
                  'yerr_noPileup'  : errors_eff_int_noPileup,\
                  'num'            : hist_efficiency_num,\
                  'num_noPileup'   : hist_efficiency_num_noPileup,\
                  'den'            : hist_efficiency_den,\
                  'bins'           : bins_efficiency,\
                  'nonzero_bins'   : nonzero_bins,\
                  'eff_curve_shift': eff_curve_shift}
    pickle.dump(pl_eff_diff, file( write_file('plots/pickle/matched_integral_%s_trigger%d_region%s.pkl' % (filename_id, triggerEt_thresh, i)), 'w+') )
    pl.savefig( write_file('plots/integral/%s_trigger%d_region%s.png' % (filename_id, triggerEt_thresh, i)) )
    pl.close()


#buildin plots for rho
fig = pl.figure(figsize=figsize)
pl.hist(data['offline_rho'], bins=bins_rho, stacked=True, fill=filled, histtype='step', color='b', label=r'offline Kt4LCTopo', linewidth=linewidth)
pl.hist(data['gFEX_rho_all'], bins=bins_rho, stacked=True, fill=filled, histtype='step', color='r', label=region_legend['all'], linewidth=linewidth)
pl.hist(data['gFEX_rho_1'],   bins=bins_rho, stacked=True, fill=filled, histtype='step', color='c', label=region_legend[1], linewidth=linewidth)
pl.hist(data['gFEX_rho_2'],   bins=bins_rho, stacked=True, fill=filled, histtype='step', color='m', label=region_legend[2], linewidth=linewidth)
pl.hist(data['gFEX_rho_3'],   bins=bins_rho, stacked=True, fill=filled, histtype='step', color='y', label=region_legend[3], linewidth=linewidth)
pl.hist(data['gFEX_rho_4'],   bins=bins_rho, stacked=True, fill=filled, histtype='step', color='k', label=region_legend[4], linewidth=linewidth)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$\rho$ [GeV]', fontsize=labelsize)
pl.ylabel(r'count', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.text(0.95, 0.05, basicTextStr['default'], transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file("plots/pileup/seed%d_noise%d_signal%d_rho.png" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh)) )
pl.close()

fig=pl.figure(figsize=figsize)
diff21 = data['gFEX_rho_2'] - data['gFEX_rho_1']
diff43 = data['gFEX_rho_4'] - data['gFEX_rho_3']
xlimit = np.ceil(np.max([np.max(diff21),np.max(diff43)])/10)*10
hist, bins = np.histogram(data['gFEX_rho_2'] - data['gFEX_rho_1'], density=True, bins=np.arange(-xlimit, xlimit+1., 1.))
fwhm = FWHM(bins, hist)
pl.plot(bins[:-1], hist, linestyle='steps-post', color='b', alpha=0.75, label='region 2 - region 1\nFWHM: %0.4f' % fwhm, linewidth=linewidth)
hist, bins = np.histogram(data['gFEX_rho_4'] - data['gFEX_rho_3'], density=True, bins=np.arange(-xlimit, xlimit+1., 1.))
fwhm = FWHM(bins, hist)
pl.plot(bins[:-1], hist, linestyle='steps-post', color='r', alpha=0.75, label='region 4 - region 3\nFWHM: %0.4f' % fwhm, linewidth=linewidth)
legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize)
legend.get_frame().set_facecolor(light_grey)
legend.get_frame().set_linewidth(0.0)
pl.xlabel(r'$\Delta\rho$ [GeV]', fontsize=labelsize)
pl.ylabel(r'normalized counts', fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.text(0.95, 0.05, basicTextStr['default'], transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file("plots/pileup/%s_deltaRho.png" % (filename_id)) )
pl.yscale('log', nonposy='clip')
pl.savefig( write_file("plots/pileup/%s_deltaRho_logy.png" % (filename_id)) )
pl.close()

corr = np.corrcoef(data['vxp_n'], data['offline_rho'])[0,1]
fig = pl.figure(figsize=figsize)
pl.hist2d(data['vxp_n'], data['offline_rho'], norm=LogNorm(), bins=(bins_vertices, bins_rho) )
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel(label_nVer, fontsize=labelsize)
pl.ylabel(label_oRho, fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.text(0.95, 0.95, '%s\nCorr = %0.4f' % (basicTextStr['default'], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/pileup/%s_offlineRho.png' % (filename_id)) )
pl.close()

for col,legend in zip(['gFEX_rho_all','gFEX_rho_1','gFEX_rho_2','gFEX_rho_3','gFEX_rho_4'],['all',1,2,3,4]):

  try:
    corr = np.corrcoef(data['vxp_n'], data[col])[0,1]
    fig = pl.figure(figsize=figsize)
    pl.hist2d(data['vxp_n'], data[col], norm=LogNorm(), bins=(bins_vertices, bins_rho) )
    cbar = pl.colorbar()
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlabel(label_nVer, fontsize=labelsize)
    pl.ylabel('%s, %s' % (label_gRho, region_legend[legend]), fontsize=labelsize)
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.text(0.95, 0.95, '%s\nCorr = %0.4f' % (basicTextStr['default'], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/pileup/%s_%s.png' % (filename_id, col)) )
    pl.close()
  except:
    print "Error for %s: could not make correlation between gFEX rho and primary vertices" % col
    pass

  try:
    corr = np.corrcoef(data['offline_rho'], data[col])[0,1]
    fig = pl.figure(figsize=figsize)
    pl.hist2d(data['offline_rho'], data[col], norm=LogNorm(), bins=(bins_rho, bins_rho) )
    cbar = pl.colorbar()
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlabel(label_oRho, fontsize=labelsize)
    pl.ylabel('%s, %s' % (label_gRho, region_legend[legend]), fontsize=labelsize)
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.text(0.95, 0.95, '%s\nCorr = %0.4f' % (basicTextStr['default'], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/pileup/%s_%s_correlation.png' % (filename_id, col)) )
    pl.close()
  except:
    print "Error for %s: could not make correlation between offline rho and gFEX rho" % col
    pass

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
