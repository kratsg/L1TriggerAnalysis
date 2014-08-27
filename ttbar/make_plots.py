from atlas_jets import *
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import os
import argparse
import time

from plots_wrapper import *

#for getting FWHM
from scipy.interpolate import UnivariateSpline, interp1d

#for getting error function fitting in diff/int curves
from scipy.optimize import curve_fit
from scipy.special import erf, erfinv

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

filename_id = "seed%d_noise%d_signal%d_digitization%d" % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed%d/leading_jets_%s.pkl" % (args.seedEt_thresh, filename_id)
data = pickle.load(file(filename))

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished reading in data:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

dataSetStr  = '$t\\bar{t}$\n$\sqrt{s}=14 \ \mathrm{TeV}\ \langle\mu\\rangle=80$'
seedCutStr  = '$E_T^\mathrm{seed} >\ %d\ \mathrm{GeV}$' % args.seedEt_thresh
noiseCutStr = '$E_T^\mathrm{tower} >\ %d\ \mathrm{GeV}$' % args.noise_filter
towerThrStr = '$\\rho\left(E_T^\mathrm{tower} <\ %d\ \mathrm{GeV}\\right)$' % args.tower_thresh

helpers = PlotHelpers(dataSetStr=dataSetStr, seedCutStr=seedCutStr, noiseCutStr=noiseCutStr, towerThrStr=towerThrStr)

tJetEt_correction = np.zeros(data.size)
tJetEt_subtracted = np.zeros(data.size)
# the new regions aren't working because there are not enough jets in them!!!
regions = {1: '', 2: '', 3: '', 4: '', '3a': '', '3b': '', '3c': '', '4a': '', '4b': '', '4c': ''}

for region in regions.keys():
  region_cut = helpers.region_cut(data['tJet.eta'], region)
  regions[region] = region_cut

  if region not in [1,2,3,4]:
    region_parsed = int(region[0])
  else:
    region_parsed = region

  if region_parsed == 3:
    rho = 'gFEX_rho_1'
  elif region_parsed == 4:
    rho = 'gFEX_rho_2'
  else:
    rho = 'gFEX_rho_{:d}'.format(region_parsed)

  region_cut = np.where(region_cut)

  tJetEt_correction[region_cut] = data['tJet.area'][region_cut]*data[rho][region_cut]
  tJetEt_subtracted[region_cut] = data['tJet.et'][region_cut] - tJetEt_correction[region_cut]

tJet_exists = data['tJet.et'] > 0.
tJet_exists_subtracted = tJetEt_subtracted > 0.

bins_efficiency = np.arange(0.,2000.,10.)
width_efficiency = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])
bins_multiplicity = np.arange(0.0,100.0,2.0)

bins_rho = np.arange(0.,100.,1.)
bins_vertices = np.arange(0.,100.,1.)

startTime_wall      = time.time()
startTime_processor = time.clock()

try:
  #multiplicity on gTowers
  fig, ax = pl.subplots(figsize=helpers.figsize)

  where = np.where( helpers.btwn(data['oJet.pt'], 0., None) )
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,\
            linestyle='steps-post',\
            alpha=0.75,\
            color=helpers.colors[0],\
            label='no cuts\n{:d} events'.format(where[0].size),\
            linewidth=helpers.linewidth)

  where = np.where( helpers.btwn(data['oJet.pt'], 100., 150.) )
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,\
            linestyle='steps-post',\
            alpha=0.75,\
            color=helpers.colors[1],\
            label='$100 < p_T^\mathrm{{oJet}} < 150$\n{:d} events'.format(where[0].size),\
            linewidth=helpers.linewidth)

  where = np.where( helpers.btwn(data['oJet.pt'], 150., 200.) )
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,\
            linestyle='steps-post',\
            alpha=0.75,\
            color=helpers.colors[2],\
            label='$150 < p_T^\mathrm{{oJet}} < 200$\n{:d} events'.format(where[0].size),\
            linewidth=helpers.linewidth)

  where = np.where( helpers.btwn(data['oJet.pt'], 200., 250.) )
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,\
            linestyle='steps-post',\
            alpha=0.75,\
            color=helpers.colors[3],\
            label='$200 < p_T^\mathrm{{oJet}} < 250$\n{:d} events'.format(where[0].size),\
            linewidth=helpers.linewidth)

  helpers.add_legend(fig, ax)
  helpers.add_labels(fig, ax, xlabel='$E_T^\mathrm{gTower}$ [GeV]', ylabel='gTower multiplicity / event')
  helpers.add_grid(fig, ax)
  helpers.add_description(fig, ax, align='bl', strings=[helpers.dataSetStr])
  ax.set_yscale('log', nonposy='clip')
  ax.set_ylim((0.0, 1284.0))
  helpers.to_file(fig, ax, 'plots/multiplicity/{}.png'.format(filename_id) )
  pl.close(fig)
except:
  print "Could not make multiplicity plot"
  pass

for i, cut in regions.iteritems():
  region = 'region_%s' % i
  region_cut = np.where(cut)

  try:
    # make a correlation of corrected Et versus trigger jet Et
    x = data['tJet.et'][region_cut]
    y = tJetEt_correction[region_cut]
    bins_x = bins_y = np.arange(0., 1500., 10.)
    xlim = (0., 1000.)
    ylim = (0., 200.)
    label_x = r'$E_T^{\mathrm{gJet}}$ [GeV]'
    label_y = r'$\rho*A^\mathrm{gJet}$ [GeV]'

    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim, ylim, profile_y=True, align='tr',
                     strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/trigger_jet_Et_correction/{}_trigger_jet_Et_correction_region{}.png'.format(filename_id, i))

    pl.close(fig)
  except:
    print "Error for {}: could not make correlation of corrected Et versus trigger jet Et".format(region)
    pass

  try:
    x = data['oJet.pt'][np.where(tJet_exists & cut)]
    y = data['tJet.et'][np.where(tJet_exists & cut)]
    bins_x = bins_y = np.arange(0., 1500., 10.)
    xlim = ylim = (0., 1000.)
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'trigger $E_T^{\mathrm{jet}}$ [GeV]'
    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim, ylim, profile_y=True, align='tr', title=r'no $\rho$ subtraction',
                     strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/jet_energy_correlation/{}_region{}.png'.format(filename_id, i))

    pl.close(fig)
  except:
    print "Error for {}: could not make jet energy correlation".format(region)
    pass

  try:
    x = data['oJet.pt'][np.where(tJet_exists_subtracted & cut)]
    y = tJetEt_subtracted[np.where(tJet_exists_subtracted & cut)]
    bins_x = bins_y = np.arange(0., 1500., 10.)
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'trigger $E_T^{\mathrm{jet}}$ [GeV]'
    xlim = ylim = (0., 1000.)
    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim, ylim, profile_y=True, align='br', title=r'with $\rho$ subtraction',
                     strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/jet_energy_correlation/{}_noPileup_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "Error for {}: could not make jet energy correlation (no Pileup)".format(region)
    pass

  # ###### START HERE #######

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
    counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution[trigger_jet_exists_noPileup[region]], bins=100, norm=LogNorm(), alpha=0.75, cmap = cmap)
    points_x, mean_y = profile_y(edges_x, offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution[trigger_jet_exists_noPileup[region]])
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.95, '%s\n$\mathrm{Corr} = %0.4f$' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
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
    counts, edges_x, edges_y, H = pl.hist2d(offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution_noPileup[trigger_jet_exists_noPileup[region]], bins=100, norm=LogNorm(), alpha=0.75, cmap = cmap)
    points_x, mean_y = profile_y(edges_x, offline_jet_Pt[region][trigger_jet_exists_noPileup[region]], resolution_noPileup[trigger_jet_exists_noPileup[region]])
    pl.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    cbar = pl.colorbar(H)
    pl.text(0.95, 0.95, '%s\n$\mathrm{Corr} = %0.4f$' % (basicTextStr[region], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
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

  points_x, mean_y = profile_y(np.arange(0.,1500.,10.), trigger_jet_Et[region], trigger_jet_Et_correction[region])
  f = interp1d(points_x, mean_y,bounds_error=False,fill_value=0.,kind='cubic')
  print region

  for triggerEt_thresh in [0., 140.]:
    print "\t",triggerEt_thresh
    '''ADD IN np.interp(triggerEt_thresh, trigger_Et, rho*A) and shift subtracted eff curves'''
    eff_curve_shift = f([triggerEt_thresh])[0]
    print "\t",eff_curve_shift
    #eff_curve_shift = np.interp(triggerEt_thresh, trigger_jet_Et[region], trigger_jet_Et_correction[region])

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
      return popt #return (a,b,c,d); width = b

    def peak_point(w):
      return  w[1]*erfinv((0.95-w[3])/w[0])+w[2]

    def make_label(w):
      if np.all(w == -1) or ~np.isfinite(peak_point(w)):
        return ''
      else:
        return '$w = {0:0.1f},\ x_{{0.95}} = {1:0.1f}$'.format(w[1], peak_point(w))

    fig = pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff, yerr=errors_eff_diff, ecolor='black', label='no subtraction\n{0}'.format(make_label(w)), linewidth=linewidth)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noShift)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noShift, yerr=errors_eff_diff_noShift, ecolor='black', label='no sub. and shift\n{0}'.format(make_label(w)), linewidth=linewidth)
    w = fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noPileup)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_diff_noPileup, yerr=errors_eff_diff_noPileup, ecolor='black', label='with subtraction\n{0}'.format(make_label(w)), linewidth=linewidth)
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
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve_int, yerr=errors_eff_int, ecolor='black', label='no subtraction\n{0}'.format(make_label(w)), linewidth=linewidth)
    w = fit_func(xpoints_efficiency[nonzero_bins]-eff_curve_shift, hist_eff_curve_int_noPileup)
    pl.errorbar(xpoints_efficiency[nonzero_bins]-eff_curve_shift, hist_eff_curve_int_noPileup, yerr=errors_eff_int_noPileup, ecolor='black', label='with subtraction\n{0}'.format(make_label(w)), linewidth=linewidth)
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
pl.hist2d(data['vxp_n'], data['offline_rho'], norm=LogNorm(), bins=(bins_vertices, bins_rho) , alpha=0.75, cmap = cmap)
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel(label_nVer, fontsize=labelsize)
pl.ylabel(label_oRho, fontsize=labelsize)
pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
pl.text(0.95, 0.95, '%s\n$\mathrm{Corr} = %0.4f$' % (basicTextStr['default'], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
pl.tick_params(axis='both', which='major', labelsize=labelsize)
pl.savefig( write_file('plots/pileup/%s_offlineRho.png' % (filename_id)) )
pl.close()

for col,legend in zip(['gFEX_rho_all','gFEX_rho_1','gFEX_rho_2','gFEX_rho_3','gFEX_rho_4'],['all',1,2,3,4]):

  try:
    corr = np.corrcoef(data['vxp_n'], data[col])[0,1]
    fig = pl.figure(figsize=figsize)
    pl.hist2d(data['vxp_n'], data[col], norm=LogNorm(), bins=(bins_vertices, bins_rho) , alpha=0.75, cmap = cmap)
    cbar = pl.colorbar()
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlabel(label_nVer, fontsize=labelsize)
    pl.ylabel('%s, %s' % (label_gRho, region_legend[legend]), fontsize=labelsize)
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.text(0.95, 0.95, '%s\n$\mathrm{Corr} = %0.4f$' % (basicTextStr['default'], corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='top', horizontalalignment='right', bbox=textprops)
    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/pileup/%s_%s.png' % (filename_id, col)) )
    pl.close()
  except:
    print "Error for %s: could not make correlation between gFEX rho and primary vertices" % col
    pass

  try:
    corr = np.corrcoef(data['offline_rho'], data[col])[0,1]
    fig = pl.figure(figsize=figsize)
    pl.hist2d(data['offline_rho'], data[col], norm=LogNorm(), bins=(bins_rho, bins_rho) , alpha=0.75, cmap = cmap)
    cbar = pl.colorbar()
    cbar.set_label('number density', fontsize=labelsize)
    pl.xlabel(label_oRho, fontsize=labelsize)
    pl.ylabel('%s, %s' % (label_gRho, region_legend[legend]), fontsize=labelsize)
    pl.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)
    pl.text(0.95, 0.05, '%s\n%s\n$\mathrm{Corr} = %0.4f$' % (dataSetStr, towerThrStr, corr), transform=fig.gca().transAxes, fontsize=labelsize, verticalalignment='bottom', horizontalalignment='right', bbox=textprops)

    #add atlas simulation
    #    internal
    pl.text(0.05, 0.95, 'ATLAS', fontsize=42, style='italic', fontweight='bold', verticalalignment='top', horizontalalignment='left', transform=fig.gca().transAxes)
    pl.text(0.27, 0.95, 'Preliminary', verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)
    pl.text(0.05, 0.90, 'Simulation', verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)

    pl.tick_params(axis='both', which='major', labelsize=labelsize)
    pl.savefig( write_file('plots/pileup/%s_%s_correlation.png' % (filename_id, col)), bbox_inches='tight')
    pl.savefig( write_file('plots/pileup/%s_%s_correlation.eps' % (filename_id, col)), bbox_inches='tight')
    pl.savefig( write_file('plots/pileup/%s_%s_correlation.pdf' % (filename_id, col)), bbox_inches='tight')
    pl.close()
  except:
    print "Error for %s: could not make correlation between offline rho and gFEX rho" % col
    pass

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
