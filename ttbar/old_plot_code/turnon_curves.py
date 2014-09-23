from atlas_jets import *
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import os
import argparse
import time

import itertools #for building up the products

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

def FWHM(bins, vals):
  spline = UnivariateSpline(bins[:-1]+np.diff(bins)/2., vals-np.max(vals)/2., s=0)
  roots = spline.roots() # find the roots
  r1, r2 = roots[0], roots[-1]
  return np.abs(r1-r2)

filename_id = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed%d/leading_jets_%s.pkl" % (args.seedEt_thresh, filename_id)

data = pickle.load(file(filename))

oJet_dtype = [\
              ('Pt','float64'),\
              ('nsj','uint64'),\
              ('subjets','object')]
tJet_dtype = [\
              ('Et','float64'),\
              ('area','float64'),\
              ('gTowers','object'),\
              ('region','uint8'),\
              ('eta','float64')]

oJets = []
tJets = []

for oJet,tJet in zip(data['leading_offline_jet'], data['matched_trigger_jet']):
  oJets.append( (oJet.Pt, oJet.nsj, oJet.subjetsPt) )
  tJets.append( (tJet.Et, tJet.area, [tJet.seed.Et] + [gTower.Et for gTower in tJet.towers_around], tJet.region(), tJet.eta) )

oJets = np.array(oJets, dtype=oJet_dtype)
tJets = np.array(tJets, dtype=tJet_dtype)

#region definition
region_legend = {'all': r'towers: $\eta \in (-4.9,4.9)$',\
                '1': r'towers: $\eta \in [-1.6,0.0)$',\
                '2': r'towers: $\eta \in [0.0,1.6)$',\
                '3': r'towers: $\eta \in (-4.9,-1.6)$',\
                '4': r'towers: $\eta \in [1.6,4.9)$',\
                '3a': r'towers: $\eta \in (-2.5,-1.6)$',\
                '3b': r'towers: $\eta \in (-3.2,-2.5)$',\
                '3c': r'towers: $\eta \in (-4.9,-3.2)$',\
                '4a': r'towers: $\eta \in (1.6,2.5)$',\
                '4b': r'towers: $\eta \in (2.5,3.2)$',\
                '4c': r'towers: $\eta \in (3.2,4.9)$'}

# the new regions aren't working because there are not enough jets in them???
all_regions = ['1','2','3','4','3a','3b','3c','4a','4b','4c']

#we added subregions to the forward region
region_1 = tJets['region'] == 1
region_2 = tJets['region'] == 2
region_3 = tJets['region'] == 3
region_4 = tJets['region'] == 4
sub_region_a = (1.6 < np.fabs(tJets['eta']))&(np.fabs(tJets['eta']) < 2.5)
sub_region_b = (2.5 < np.fabs(tJets['eta']))&(np.fabs(tJets['eta']) < 3.2)
sub_region_c = (3.2 < np.fabs(tJets['eta']))&(np.fabs(tJets['eta']) < 4.9)

region = {}
region['1']  = np.where( region_1 )
region['2']  = np.where( region_2 )
region['3']  = np.where( region_3 )
region['4']  = np.where( region_4 )
region['3a'] = np.where( sub_region_a&region_3 )
region['3b'] = np.where( sub_region_b&region_3 )
region['3c'] = np.where( sub_region_c&region_3 )
region['4a'] = np.where( sub_region_a&region_4 )
region['4b'] = np.where( sub_region_b&region_4 )
region['4c'] = np.where( sub_region_c&region_4 )

#styling options for the plots
figsize = (16, 12)
labelsize = 28
titlesize = 36
linewidth = 2
light_grey = np.array([float(200)/float(255)]*3)
filled=False
textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

dataSetStr  = 'TTbar 14TeV $\langle\mu\\rangle=80$'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

basicTextStr = {}
basicTextStr['default'] = '%s\n%d Events\n%s\n%s\n%s' % (dataSetStr, oJets.size, seedCutStr, noiseCutStr, towerThrStr)
for i in all_regions:
  basicTextStr['%s' % i] = '%s\n%d Events \n%s\n%s\n%s' % (dataSetStr, oJets[region[i]].size, seedCutStr, noiseCutStr, towerThrStr)

region_legend = {'all': r'towers: $\eta \in (-4.9,4.9)$',\
                '1': r'towers: $\eta \in [-1.6,0.0)$',\
                '2': r'towers: $\eta \in [0.0,1.6)$',\
                '3': r'towers: $\eta \in (-4.9,-1.6)$',\
                '4': r'towers: $\eta \in [1.6,4.9)$',\
                '3a': r'towers: $\eta \in (-2.5,-1.6)$',\
                '3b': r'towers: $\eta \in (-3.2,-2.5)$',\
                '3c': r'towers: $\eta \in (-4.9,-3.2)$',\
                '4a': r'towers: $\eta \in (1.6,2.5)$',\
                '4b': r'towers: $\eta \in (2.5,3.2)$',\
                '4c': r'towers: $\eta \in (3.2,4.9)$'}

bins_efficiency = np.arange(0.,1000.,10.)
width_efficiency = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])
xpoints_efficiency = bins_efficiency[:-1] + width_efficiency/2.

startTime_wall      = time.time()
startTime_processor = time.clock()


#time to build the huge fucking loops
''' each plot:
      - is a different trigger jet Et requirement
      - is a different offline subjet requirement
      - is a different region
      - contains the different combinations of gTower requirements
'''
triggerEt_thresholds = [0., 100., 110., 130., 150., 200., 220.] #GeV
offline_subjetPt_thresholds = [15., 20., 25., 30.] #GeV
offline_nsj_thresholds = [2, 3] #applied after subjetPt cut
tJet_gTower_f_thresholds = [0.0, 0.05, 0.10, 0.15] # E_T^i / E_T^{tJet}
tJet_gTower_N_thresholds = [0, 2, 3] # number of gTowers >= f_threshold

markers = {'0': '<', '2': 'v', '3': '>'} #markers for tJet_gTower_f_thresholds
colors = {'0.0': 'k', '0.05': 'r', '0.1': 'b', '0.15': 'g'} #colors for tJet_gTower_N_thresholds

def custom_filter(listOfObjects, cut, numLeft):
  ''' this can be called in two ways
      custom_filter(subjets, subjetPt_cut, nsj_cut)
      custom_filter(gTowers/tJet.Et, f_cut, gTower_N_cut)
  '''
  return (np.sum(np.array(listOfObjects) > cut) >= numLeft)

xlim_efficiency = (0.0, 750.0) #GeV
ylim_efficiency = (0.0,1.1)

#define erfx used for error fitting
def func(x, a, b, c, d):
  # note that b == sigma here, see wiki for more info
  return a*erf( (x-c)/b ) + d

for regionStr, tJet_Et_cut, oJet_subjetPt_cut, oJet_nsj_cut in itertools.product(*[all_regions, triggerEt_thresholds, offline_subjetPt_thresholds, offline_nsj_thresholds]):

  print regionStr, tJet_Et_cut, oJet_subjetPt_cut, oJet_nsj_cut

  textstr = '%s\n$E_T^\mathrm{gFEX\ jet} >\ %d\ \mathrm{GeV}$' % (basicTextStr[regionStr], tJet_Et_cut)

  def fit_func(x, y):
    try:
      popt, pcov = curve_fit(func, x, y, p0=(1., 50., tJet_Et_cut, 0.5))
    except RuntimeError:
      return -1.
    return popt[1] #return the width (b)

  # build up the relevant data for offline Pt and trigger Et
  offline_jet_Pt = oJets['Pt'][region[regionStr]]
  trigger_jet_Et = tJets['Et'][region[regionStr]]

  tJet_Et_correction = tJets['area']*data['gFEX_rho_%s' %  regionStr[0]]
  tJet_Et_correction = tJet_Et_correction[region[regionStr]]
  eff_curve_shift = np.interp(tJet_Et_cut, trigger_jet_Et, tJet_Et_correction)

  # build up the same relevant data for offline Pt and trigger Et with pileup subtracted
  #     also filter out cases in which we have a negative trigger jet
  trigger_jet_Et_noPileup = trigger_jet_Et - tJet_Et_correction
  offline_jet_Pt_noPileup = offline_jet_Pt[np.where(trigger_jet_Et_noPileup > 0.)]
  trigger_jet_Et_noPileup = trigger_jet_Et_noPileup[np.where(trigger_jet_Et_noPileup > 0.)]

  # apply the subjetPt cut
  oJet_subjet_cut = np.array(map(lambda x: custom_filter(x, oJet_subjetPt_cut, oJet_nsj_cut), oJets['subjets'][region[regionStr]]))
  oJet_subjet_cut_noPileup = oJet_subjet_cut[np.where(trigger_jet_Et_noPileup > 0.)]

  offline_jet_Pt = offline_jet_Pt[np.where(oJet_subjet_cut)]
  trigger_jet_Et = trigger_jet_Et[np.where(oJet_subjet_cut)]

  offline_jet_Pt_noPileup = offline_jet_Pt_noPileup[np.where(oJet_subjet_cut_noPileup)]
  trigger_jet_Et_noPileup = trigger_jet_Et_noPileup[np.where(oJet_subjet_cut_noPileup)]

  # build the denominator
  hist_efficiency_den, _            = np.histogram(offline_jet_Pt, bins=bins_efficiency)
  hist_efficiency_den_noPileup, _   = np.histogram(offline_jet_Pt_noPileup, bins=bins_efficiency)

  # start making the plot here, finish up after the loop
  fig = pl.figure(figsize=(figsize[0]+4, figsize[1]))

  pl_eff_diff = []

  for tJet_f_cut, tJet_N_cut in itertools.product(*[tJet_gTower_f_thresholds, tJet_gTower_N_thresholds]):
    print "\t", tJet_f_cut, tJet_N_cut

    # build the numerator, use other cuts
    tJet_towers_cut = []
    for gTowers, tJet_Et in zip(tJets['gTowers'][region[regionStr]][np.where(oJet_subjet_cut)], trigger_jet_Et):
      tJet_towers_cut.append( custom_filter(gTowers/tJet_Et, tJet_f_cut, tJet_N_cut) )
    tJet_towers_cut = np.array(tJet_towers_cut)
    tJet_towers_cut_noPileup = tJet_towers_cut[np.where(trigger_jet_Et_noPileup > 0.)]

    print "\t\t", offline_jet_Pt.size, offline_jet_Pt_noPileup.size 
    print "\t\t", offline_jet_Pt[np.where( (trigger_jet_Et > tJet_Et_cut + eff_curve_shift)&(tJet_towers_cut) )].size, offline_jet_Pt_noPileup[np.where( (trigger_jet_Et_noPileup > tJet_Et_cut)&(tJet_towers_cut_noPileup) )].size
    hist_efficiency_num, _          = np.histogram(offline_jet_Pt[np.where( (trigger_jet_Et > tJet_Et_cut + eff_curve_shift)&(tJet_towers_cut) )], bins=bins_efficiency)
    hist_efficiency_num_noPileup, _ = np.histogram(offline_jet_Pt_noPileup[np.where( (trigger_jet_Et_noPileup > tJet_Et_cut)&(tJet_towers_cut_noPileup) )], bins=bins_efficiency)

    nonzero_bins = np.where(hist_efficiency_den != 0)
    nonzero_bins_noPileup = np.where(hist_efficiency_den_noPileup != 0)
    #compute differential curves
    hist_eff_curve_diff             = np.true_divide(hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    hist_eff_curve_diff_noPileup    = np.true_divide(hist_efficiency_num_noPileup[nonzero_bins_noPileup], hist_efficiency_den_noPileup[nonzero_bins_noPileup])

    #hist_eff_curve_int              = np.true_divide(np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
    #hist_eff_curve_int_noPileup     = np.true_divide(np.cumsum(hist_efficiency_num_noPileup[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])

    #binomial errors s^2 = n * p * q
    errors_eff_diff                 = binomial_errors(hist_eff_curve_diff, hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    errors_eff_diff_noPileup        = binomial_errors(hist_eff_curve_diff_noPileup, hist_efficiency_num_noPileup[nonzero_bins_noPileup], hist_efficiency_den[nonzero_bins_noPileup])

    #errors_eff_int                  = binomial_errors(hist_eff_curve_int, np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
    #errors_eff_int_noPileup         = binomial_errors(hist_eff_curve_int_noPileup, np.cumsum(hist_efficiency_num_noPileup[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])

    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff, linestyle='--', yerr=errors_eff_diff, ecolor='black', marker=markers[str(tJet_N_cut)], c=colors[str(tJet_f_cut)], mfc=colors[str(tJet_f_cut)], ms=10, mec=colors[str(tJet_f_cut)], label='no sub.\nWidth = %0.4f\n$f_\mathrm{gTower} = %0.2f,\ n \geq %d$' % (w, tJet_f_cut, tJet_N_cut), linewidth=linewidth, alpha=0.75)
    w = fit_func(xpoints_efficiency[nonzero_bins_noPileup], hist_eff_curve_diff_noPileup)
    pl.errorbar(xpoints_efficiency[nonzero_bins_noPileup], hist_eff_curve_diff_noPileup, yerr=errors_eff_diff_noPileup, ecolor='black', marker=markers[str(tJet_N_cut)], c=colors[str(tJet_f_cut)], mfc=colors[str(tJet_f_cut)], ms=10, mec=colors[str(tJet_f_cut)], label='w/ sub.\nWidth = %0.4f\n$f_\mathrm{gTower} = %0.2f,\ n \geq %d$' % (w, tJet_f_cut, tJet_N_cut), linewidth=linewidth, alpha=0.75)

    pl_eff_diff.append({'xdata'     : xpoints_efficiency,\
                   'ydata'          : hist_eff_curve_diff,\
                   'ydata_noPileup' : hist_eff_curve_diff_noPileup,\
                   'xerr'           : 1.0,\
                   'yerr'           : errors_eff_diff,\
                   'yerr_noPileup'  : errors_eff_diff_noPileup,\
                   'num'            : hist_efficiency_num,\
                   'num_noPileup'   : hist_efficiency_num_noPileup,\
                   'den'            : hist_efficiency_den,\
                   'bins'           : bins_efficiency,\
                   'nonzero_bins'   : nonzero_bins,\
                   'nonzero_bins_noPileup': nonzero_bins_noPileup,\
                   'eff_curve_shift': eff_curve_shift})


  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
  pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
  pl.text(0.95, 0.05, '%s\nShift: %0.4f' % (textstr, eff_curve_shift), transform=fig.gca().transAxes, fontsize=labelsize/2, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
  pl.xlim(xlim_efficiency)
  pl.ylim(ylim_efficiency)
  legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize/2, ncol=6, loc=3, bbox_to_anchor=(-0.05, 1.02, 1., .102)) #expand top
  legend.get_frame().set_facecolor(light_grey)
  legend.get_frame().set_linewidth(0.0)
  pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
  pl.tick_params(axis='both', which='major', labelsize=labelsize)

  pickle.dump(pl_eff_diff, file( write_file('plots/pickle/differential_%s_region%s_tJetEtcut%d_oJetsubjetPt%d_oJetnsj%d.pkl' % (filename_id, regionStr, tJet_Et_cut, oJet_subjetPt_cut, oJet_nsj_cut)), 'w+') )
  pl.savefig( write_file('plots/differential/%s_region%s_tJetEtcut%d_oJetsubjetPt%d_oJetnsj%d.pdf' % (filename_id, regionStr, tJet_Et_cut, oJet_subjetPt_cut, oJet_nsj_cut)) , bbox_extra_artists=(legend,), bbox_inches='tight')
  pl.close()


endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
