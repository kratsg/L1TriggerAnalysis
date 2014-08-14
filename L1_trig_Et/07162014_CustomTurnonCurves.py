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

filename_id_noNoiseCut = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, 0, args.tower_thresh, args.digitization)
filename_id_wNoiseCut  = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename_noNoiseCut    = "data/seed%d/leading_jets_%s.pkl" % (15, filename_id_noNoiseCut)
filename_wNoiseCut     = "data/seed%d/leading_jets_%s.pkl" % (15, filename_id_wNoiseCut)

oJet_dtype = [\
              ('m','float64'),\
              ('Pt','float64'),\
              ('nsj','uint64'),\
              ('subjets','object')]
tJet_dtype = [\
              ('Et','float64'),\
              ('area','float64')]

data = pickle.load(file(filename_noNoiseCut))
oJets_noNoiseCut = []
tJets_noNoiseCut = []
rho_noNoiseCut   = []
for oJet,tJet,rho in zip(data['leading_offline_jet'], data['matched_trigger_jet'], data['gFEX_rho_all']):
  if abs(oJet.eta) > 1.6:
    continue
  oJets_noNoiseCut.append( (oJet.m, oJet.Pt, oJet.nsj, oJet.subjetsPt) )
  tJets_noNoiseCut.append( (tJet.Et, tJet.area) )
  rho_noNoiseCut.append( rho )
oJets_noNoiseCut = np.array(oJets_noNoiseCut, dtype=oJet_dtype)
tJets_noNoiseCut = np.array(tJets_noNoiseCut, dtype=tJet_dtype)
rho_noNoiseCut   = np.array(rho_noNoiseCut, dtype=np.float)
del data

data = pickle.load(file(filename_wNoiseCut))
oJets_wNoiseCut = []
tJets_wNoiseCut = []
rho_wNoiseCut   = []
for oJet,tJet,rho in zip(data['leading_offline_jet'], data['matched_trigger_jet'], data['gFEX_rho_all']):
  if abs(oJet.eta) > 1.6:
    continue
  oJets_wNoiseCut.append( (oJet.m, oJet.Pt, oJet.nsj, oJet.subjetsPt) )
  tJets_wNoiseCut.append( (tJet.Et, tJet.area) )
  rho_wNoiseCut.append( rho )
oJets_wNoiseCut = np.array(oJets_wNoiseCut, dtype=oJet_dtype)
tJets_wNoiseCut = np.array(tJets_wNoiseCut, dtype=tJet_dtype)
rho_wNoiseCut   = np.array(rho_noNoiseCut, dtype=np.float)
del data

#styling options for the plots
figsize = (16, 12)
labelsize = 28
titlesize = 36
linewidth = 2
light_grey = np.array([float(200)/float(255)]*3)
filled=False
textprops = dict(boxstyle='round', facecolor=light_grey, alpha=0.5, linewidth=0.0)

dataSetStr  = r'$t\bar{t}\ \sqrt{s}=14$ TeV $\langle\mu\rangle=80$'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

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
offline_subjetPt_thresholds = [20., 20., 20.] #GeV
offline_nsj_thresholds = [1, 2, 3] #applied after subjetPt cut

def custom_filter(listOfObjects, cut, numLeft):
  ''' this can be called in two ways
      custom_filter(subjets, subjetPt_cut, nsj_cut)
      custom_filter(gTowers/tJet.Et, f_cut, gTower_N_cut)
  '''
  return (np.sum(np.array(listOfObjects) > cut) >= numLeft)

xlim_efficiency = (0.0, 450.0) #GeV
ylim_efficiency = (0.0,1.1)

#define erfx used for error fitting
def func(x, a, b, c, d):
  # note that b == sigma here, see wiki for more info
  return a*erf( (x-c)/b ) + d

'''different plots required:

require |eta_offline| <= 1.6
require seed cut = 15 GeV
x3 for each digitization (125, 256, 512)
isolated offline jets dR > 2 and matched trigger jet dR < 1

  - R = 1.0 trimmed jets and >= 1 subjet with P_T^subjet > 20 GeV
  - R = 1.0 trimmed jets and >= 1 subjet with P_T^subjet > 20 GeV and 100 <= m_oJet <= 220 GeV
  - R = 1.0 trimmed jets and >= 2 subjet with P_T^subjet > 20 GeV and 100 <= m_oJet <= 220 GeV
  - R = 1.0 trimmed jets and >= 3 subjet with P_T^subjet > 20 GeV and 100 <= m_oJet <= 220 GeV


  numerator iterated over (for each above):
    - no subtraction and G140 and no shift
    - no subtraction and G140 and E_T^gTower < 3 GeV noise cut
    - no subtraction and G140 and shift
    - w/ subtraction and G140 (upper thresh of 6 GeV)

'''
subjets_cut_noNoiseCut = {}
subjets_cut_wNoiseCut = {}
for oJet_subjetPt_cut, oJet_nsj_cut in zip(offline_subjetPt_thresholds, offline_nsj_thresholds):
  subjets_cut_noNoiseCut[oJet_nsj_cut] = np.array(map(lambda x: custom_filter(x, oJet_subjetPt_cut, oJet_nsj_cut), oJets_noNoiseCut['subjets']))
  subjets_cut_wNoiseCut[oJet_nsj_cut] = np.array(map(lambda x: custom_filter(x, oJet_subjetPt_cut, oJet_nsj_cut), oJets_wNoiseCut['subjets']))

mass_cut_noNoiseCut = (100. <= oJets_noNoiseCut['m'])&(oJets_noNoiseCut['m'] <= 220.)
mass_cut_wNoiseCut = (100. <= oJets_wNoiseCut['m'])&(oJets_wNoiseCut['m'] <= 220.)

offline_cuts_noNoiseCut = [\
  subjets_cut_noNoiseCut[1],\
  subjets_cut_noNoiseCut[1]&mass_cut_noNoiseCut,\
  subjets_cut_noNoiseCut[2]&mass_cut_noNoiseCut,\
  subjets_cut_noNoiseCut[3]&mass_cut_noNoiseCut]

offline_cuts_wNoiseCut = [\
  subjets_cut_wNoiseCut[1],\
  subjets_cut_wNoiseCut[1]&mass_cut_wNoiseCut,\
  subjets_cut_wNoiseCut[2]&mass_cut_wNoiseCut,\
  subjets_cut_wNoiseCut[3]&mass_cut_wNoiseCut]

tJet_Et_cut = 140.

def fit_func(x, y):
  try:
    popt, pcov = curve_fit(func, x, y, p0=(1., 50., tJet_Et_cut, 0.5))
  except RuntimeError:
    return -1.
  return popt[1] #return the width (b)


plotStrings = [\
  '$|\eta^\mathrm{oJet}| \leq 1.6$\nanti-Kt R=1.0 (5% trimmed) jets\n$\geq 1 D=0.3$ subjet w/ $P_T^\mathrm{subjet} > 20$ GeV',\
  '$|\eta^\mathrm{oJet}| \leq 1.6$\n$100 \leq m^\mathrm{oJet} \leq 220$ GeV\nanti-Kt R=1.0 (5% trimmed) jets\n$\geq 1 D=0.3$ subjet w/ $P_T^\mathrm{subjet} > 20$ GeV',\
  '$|\eta^\mathrm{oJet}| \leq 1.6$\n$100 \leq m^\mathrm{oJet} \leq 220$ GeV\nanti-Kt R=1.0 (5% trimmed) jets\n$\geq 2 D=0.3$ subjet w/ $P_T^\mathrm{subjet} > 20$ GeV',\
  '$|\eta^\mathrm{oJet}| \leq 1.6$\n$100 \leq m^\mathrm{oJet} \leq 220$ GeV\nanti-Kt R=1.0 (5% trimmed) jets\n$\geq 3 D=0.3$ subjet w/ $P_T^\mathrm{subjet} > 20$ GeV']

nn = 0
for cuts_noNoiseCut, cuts_wNoiseCut, plotStr in zip(offline_cuts_noNoiseCut, offline_cuts_wNoiseCut, plotStrings):

  # prepare the cuts by np.where(Truth-arrays)
  cuts_noNoiseCut = np.where(cuts_noNoiseCut)
  cuts_wNoiseCut  = np.where(cuts_wNoiseCut)

  # build up the relevant data for offline Pt and trigger Et
  offline_jet_Pt_noNoiseCut = oJets_noNoiseCut['Pt'][cuts_noNoiseCut]
  offline_jet_Pt_wNoiseCut  = oJets_wNoiseCut['Pt'][cuts_wNoiseCut]

  trigger_jet_Et_noNoiseCut = tJets_noNoiseCut['Et'][cuts_noNoiseCut]
  trigger_jet_Et_wNoiseCut = tJets_wNoiseCut['Et'][cuts_wNoiseCut]

  tJet_Et_correction_noNoiseCut = (tJets_noNoiseCut['area']*rho_noNoiseCut)[cuts_noNoiseCut]
  #tJet_Et_correction_wNoiseCut = (tJets_wNoiseCut['area']*rho_wNoiseCut)[cuts_wNoiseCut]

  eff_curve_shift_noNoiseCut = np.interp(tJet_Et_cut, trigger_jet_Et_noNoiseCut, tJet_Et_correction_noNoiseCut)
  #eff_curve_shift_wNoiseCut = np.interp(tJet_Et_cut, trigger_jet_Et_wNoiseCut, tJet_Et_correction_wNoiseCut)

  trigger_jet_Et_noNoiseCut_noPileup = trigger_jet_Et_noNoiseCut - tJet_Et_correction_noNoiseCut
  offline_jet_Pt_noNoiseCut_noPileup = offline_jet_Pt_noNoiseCut[np.where(trigger_jet_Et_noNoiseCut_noPileup > 0.)]
  trigger_jet_Et_noNoiseCut_noPileup = trigger_jet_Et_noNoiseCut_noPileup[np.where(trigger_jet_Et_noNoiseCut_noPileup > 0.)]

  hist_efficiency_den = [\
    np.histogram(offline_jet_Pt_noNoiseCut, bins=bins_efficiency)[0],\
    np.histogram(offline_jet_Pt_wNoiseCut, bins=bins_efficiency)[0],\
    np.histogram(offline_jet_Pt_noNoiseCut, bins=bins_efficiency)[0],\
    np.histogram(offline_jet_Pt_noNoiseCut_noPileup, bins=bins_efficiency)[0] ]

  hist_efficiency_num = [\
    np.histogram(offline_jet_Pt_noNoiseCut[np.where( (trigger_jet_Et_noNoiseCut > tJet_Et_cut) )], bins=bins_efficiency)[0],\
    np.histogram(offline_jet_Pt_wNoiseCut[np.where( (trigger_jet_Et_wNoiseCut > tJet_Et_cut) )], bins=bins_efficiency)[0],\
    np.histogram(offline_jet_Pt_noNoiseCut[np.where( (trigger_jet_Et_noNoiseCut > tJet_Et_cut + eff_curve_shift_noNoiseCut) )], bins=bins_efficiency)[0],\
    np.histogram(offline_jet_Pt_noNoiseCut_noPileup[np.where( (trigger_jet_Et_noNoiseCut_noPileup > tJet_Et_cut) )], bins=bins_efficiency)[0] ]

  labels = [\
    'no subtraction',\
    'no subtraction\n$E_T^\mathrm{gTower} > 3$ GeV (noise cut)',\
    'no subtraction, shifted',\
    'w/ subtraction\n'+r'$\rho(E_T^\mathrm{gTower} < 6$ GeV$)$' ]

  numEvents = [\
    offline_jet_Pt_noNoiseCut.size,\
    offline_jet_Pt_wNoiseCut.size,\
    offline_jet_Pt_noNoiseCut.size,\
    offline_jet_Pt_noNoiseCut_noPileup.size ]

  markers = ['<','v','>','^']
  colors = ['k','r','b','g']


  # start making the plot here, finish up after the loop
  fig = pl.figure(figsize=(figsize[0]+4, figsize[1]))

  for den, num, label, numEvent, marker, color in zip(hist_efficiency_den, hist_efficiency_num, labels, numEvents, markers, colors):

    textstr = '%s\n%d Events\n%s\n$E_T^\mathrm{gFEX\ jet} >\ %d\ \mathrm{GeV}$\n%s' % (dataSetStr, numEvent, seedCutStr, tJet_Et_cut, plotStr)

    nonzero_bins = np.where(den != 0)
    #compute differential curves
    hist_eff_curve_diff             = np.true_divide(num[nonzero_bins], den[nonzero_bins])

    #binomial errors s^2 = n * p * q
    errors_eff_diff                 = binomial_errors(hist_eff_curve_diff, num[nonzero_bins], den[nonzero_bins])

    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff, linestyle='-', yerr=errors_eff_diff, ecolor='black', marker=marker, c=color, mfc=color, ms=10, mec=color, label='%s\nWidth = %0.4f' % (label, w), linewidth=linewidth, alpha=0.75)

  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
  pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
  pl.text(0.95, 0.05, '%s\nShift: %0.4f' % (textstr, eff_curve_shift_noNoiseCut), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)
  pl.xlim(xlim_efficiency)
  pl.ylim(ylim_efficiency)
  legend = pl.legend(fancybox=True, framealpha=0.75, fontsize=labelsize, ncol=2, loc=3, mode='expand', bbox_to_anchor=(0.0, 1.02, 1., .102)) #expand top
  legend.get_frame().set_facecolor(light_grey)
  legend.get_frame().set_linewidth(0.0)
  pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
  pl.tick_params(axis='both', which='major', labelsize=labelsize)

  pl.savefig( write_file('plots/differential_hardcode/%s_%d.pdf' % (filename_id_wNoiseCut, nn)) , bbox_extra_artists=(legend,), bbox_inches='tight')
  pl.close()
  nn += 1


endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
