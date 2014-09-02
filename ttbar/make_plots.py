from atlas_jets import *
from plots_wrapper import *
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import argparse
import time

# for getting FWHM
from scipy.interpolate import interp1d

# for getting error function fitting in diff/int curves
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
print "Finished reading in data:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

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

  if region not in [1, 2, 3, 4]:
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

bins_efficiency = np.arange(0., 2000., 10.)
width_efficiency = np.array([x - bins_efficiency[i-1] for i, x in enumerate(bins_efficiency)][1:])
bins_multiplicity = np.arange(0.0, 100.0, 2.0)

bins_rho = np.arange(0., 100., 1.)
bins_vertices = np.arange(0., 100., 1.)

startTime_wall = time.time()
startTime_processor = time.clock()

try:
  # multiplicity on gTowers
  fig, ax = pl.subplots(figsize=helpers.figsize)

  where = np.where(helpers.btwn(data['oJet.pt'], 0., None))
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,
          linestyle='steps-post',
          alpha=0.75,
          color=helpers.colors[0],
          label='no cuts\n{:d} events'.format(where[0].size),
          linewidth=helpers.linewidth)

  where = np.where(helpers.btwn(data['oJet.pt'], 100., 150.))
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,
          linestyle='steps-post',
          alpha=0.75,
          color=helpers.colors[1],
          label='$100 < p_T^\mathrm{{oJet}} < 150$\n{:d} events'.format(where[0].size),
          linewidth=helpers.linewidth)

  where = np.where(helpers.btwn(data['oJet.pt'], 150., 200.))
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,
          linestyle='steps-post',
          alpha=0.75,
          color=helpers.colors[2],
          label='$150 < p_T^\mathrm{{oJet}} < 200$\n{:d} events'.format(where[0].size),
          linewidth=helpers.linewidth)

  where = np.where(helpers.btwn(data['oJet.pt'], 200., 250.))
  ax.plot(bins_multiplicity[:-1], np.cumsum(np.sum(data['gTower_distribution'][where]).astype(float)[::-1])[::-1]/where[0].size,
          linestyle='steps-post',
          alpha=0.75,
          color=helpers.colors[3],
          label='$200 < p_T^\mathrm{{oJet}} < 250$\n{:d} events'.format(where[0].size),
          linewidth=helpers.linewidth)

  helpers.add_legend(fig, ax)
  helpers.add_labels(fig, ax,
                     xlabel='$E_T^\mathrm{gTower}$ [GeV]',
                     ylabel='gTower multiplicity / event')
  helpers.add_grid(fig, ax)
  helpers.add_description(fig, ax,
                          align='bl',
                          strings=[helpers.dataSetStr])
  ax.set_yscale('log', nonposy='clip')
  ax.set_ylim((0.0, 1284.0))
  helpers.to_file(fig, ax, 'plots/multiplicity/{}.png'.format(filename_id))
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

    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_y=True, align='tr',
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
    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_y=True, align='tr', title=r'no $\rho$ subtraction',
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
    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_y=True, align='br', title=r'with $\rho$ subtraction',
                     strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/jet_energy_correlation/{}_noPileup_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "Error for {}: could not make jet energy correlation (no Pileup)".format(region)
    pass

  # we want to make resolution plots now
  # note that we ignore trigger jets that are negative or less than zero
  resolution = (data['tJet.et'][np.where(cut)]/data['oJet.pt'][np.where(cut)]) - 1.0
  resolution_subtracted = (tJetEt_subtracted[np.where(cut)]/data['oJet.pt'][np.where(cut)]) - 1.0
  # eta = np.array([oJet.eta for oJet in data['leading_offline_jet']])
  # rho = data['offline_rho']
  # vxp_n = data['vxp_n']

  try:
    x = data['oJet.pt'][np.where(tJet_exists_subtracted & cut)]
    y = resolution[np.where(tJet_exists_subtracted & cut)]
    bins_x = bins_y = 100
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$'
    xlim = (0., 1000.)
    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, profile_y=True, align='tr', title=r'no $\rho$ subtraction',
                     strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "Error for {}: could not make resolution correlation".format(region)
    pass

  try:
    x = data['oJet.pt'][np.where(tJet_exists_subtracted & cut)]
    y = resolution_subtracted[np.where(tJet_exists_subtracted & cut)]
    bins_x = bins_y = 100
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$'
    xlim = (0., 1000.)
    fig, ax = corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, profile_y=True, align='tr', title=r'with $\rho$ subtraction',
                     strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_withRhoSubtraction_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "Error for {}: could not make resolution correlation (with Rho subtraction)".format(region)
    pass

  try:
    # y projection slices
    pl_res_proj = {}
    fig, ax = pl.subplots(figsize=self.figsize)
    for oJetPt_cuts in [(170., 180.), (200., 220.), (300., 350.)]:
      oJetPt_cut = helpers.btwn(data['oJet.pt'], oJetPt_cuts[0], oJetPt_cuts[1])
      hist, bins = np.histogram(resolution[np.where(cut & oJetPt_cut & tJet_exists_subtracted)], bins=100, density=True)
      fwhm = FWHM(bins, hist)
      ax.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='b', label=r'${:d}\ \mathrm{{GeV}} < p_T^\mathrm{{oJet}} <\ {:d}\ \mathrm{{GeV}}$\nFWHM = {:0.4f}'.format(oJetPt_cuts[0], oJetPt_cuts[1], fwhm), linewidth=helpers.linewidth)

      pl_res_proj['{:d}to{:d}'.format(oJetPt_cuts[0], oJetPt_cuts[1])] = resolution[np.where(cut & oJetPt_cut & tJet_exists_subtracted)]

    helpers.add_legend(fig, ax)
    helpers.add_labels(fig, ax, xlabel=r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', ylabel='normalized counts', title='Y-Axis Projections of Resolution')
    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
    ax.set_xlim((-1.0, 1.0))
    helpers.add_grid(fig, ax)
    pickle.dump(pl_res_proj, file(helpers.write_file('plots/pickle/{}_resolution_PtOffline_projection_region{}.pkl'.format(filename_id, i)), 'w+'))
    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_projection_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "Error for {}: could not resolution projection".format(region)
    pass

  try:
    pl_res_proj = {}
    fig, ax = pl.subplots(figsize=self.figsize)
    for oJetPt_cuts in [(170., 180.), (200., 220.), (300., 350.)]:
      oJetPt_cut = helpers.btwn(data['oJet.pt'], oJetPt_cuts[0], oJetPt_cuts[1])
      hist, bins = np.histogram(resolution_subtracted[np.where(cut & oJetPt_cut & tJet_exists_subtracted)], bins=100, density=True)
      fwhm = FWHM(bins, hist)
      ax.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='b', label=r'${:d}\ \mathrm{{GeV}} < p_T^\mathrm{{oJet}} <\ {:d}\ \mathrm{{GeV}}$\nFWHM = {:0.4f}'.format(oJetPt_cuts[0], oJetPt_cuts[1], fwhm), linewidth=helpers.linewidth)

      pl_res_proj['{:d}to{:d}'.format(oJetPt_cuts[0], oJetPt_cuts[1])] = resolution_subtracted[np.where(cut & oJetPt_cut & tJet_exists_subtracted)]

    helpers.add_legend(fig, ax)
    helpers.add_labels(fig, ax, xlabel=r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', ylabel='normalized counts', title=r'Y-Axis Projections of Resolution, with $\rho$ subtraction')
    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
    ax.set_xlim((-1.0, 1.0))
    helpers.add_grid(fig, ax)
    pickle.dump(pl_res_proj, file(helpers.write_file('plots/pickle/{}_resolution_PtOffline_withRhoSubtraction_projection_region{}.pkl'.format(filename_id, i)), 'w+'))
    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_withRhoSubtraction_projection_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "Error for {}: could not make resolution projection (no Pileup)".format(region)
    pass

  # ###### START HERE ######

  points_x, mean_y = helpers.profile_y(np.arange(0., 1500., 10.), data['tJet.et'][np.where(cut)], tJetEt_correction[np.where(cut)])
  f = interp1d(points_x, mean_y, bounds_error=False, fill_value=0., kind='cubic')
  print region

  for tJetEt_cut in [0., 140.]:
    print "\t", tJetEt_cut
    # ADD IN np.interp(tJetEt_cut, trigger_Et, rho*A) and shift subtracted eff curves
    eff_curve_shift = f([tJetEt_cut])[0]
    print "\t", eff_curve_shift
    # eff_curve_shift = np.interp(tJetEt_cut, trigger_jet_Et[region], trigger_jet_Et_correction[region])
    tJetEtThrStr = r'$E_T^\mathrm{{gFEX\ jet}} >\ {:0.2f}\ \mathrm"{GeV}}$'.format(tJetEt_cut)

    pl_eff_diff = {'eff_curve_shift': eff_curve_shift}
    fig, ax = pl.subplots(figsize=helpers.figsize)

    # helpers.add_turnon() returns xpoints_efficiency, hist_eff_curve, errors_eff, nonzero_bins, w
    x, y, e, n, w = helpers.add_turnon(fig, ax, data=data['oJet.pt'], den_cut=cut, num_cut=(data['tJet.et'] > tJetEt_cut + eff_curve_shift), label='no subtraction', p0=(1., 50., tJetEt_cut, 0.5))
    pl_eff_diff['xdata'] = x
    pl_eff_diff['ydata'] = y
    pl_eff_diff['xerr'] = 1.0
    pl_eff_diff['yerr'] = e
    pl_eff_diff['nonzero_bins'] = n
    pl_eff_diff['fit'] = w
    x, y, e, n, w = helpers.add_turnon(fig, ax, data=data['oJet.pt'], den_cut=cut, num_cut=(data['tJet.et'] > tJetEt_cut), label='no subt., no shift', p0=(1., 50., tJetEt_cut, 0.5))
    pl_eff_diff['xdata_noShift'] = x
    pl_eff_diff['ydata_noShift'] = y
    pl_eff_diff['xerr_noShift'] = 1.0
    pl_eff_diff['yerr_noShift'] = e
    pl_eff_diff['nonzero_bins_noShift'] = n
    pl_eff_diff['fit_noShift'] = w
    x, y, e, n, w = helpers.add_turnon(fig, ax, data=data['oJet.pt'], den_cut=cut, num_cut=(tJetEt_subtracted > tJetEt_cut), label='with subtraction', p0=(1., 50., tJetEt_cut, 0.5))
    pl_eff_diff['xdata_subtracted'] = x
    pl_eff_diff['ydata_subtracted'] = y
    pl_eff_diff['xerr_subtracted'] = 1.0
    pl_eff_diff['yerr_subtracted'] = e
    pl_eff_diff['nonzero_bins_subtracted'] = n
    pl_eff_diff['fit_subtracted'] = w

    ax.set_xlim((0.0, 450.0))
    ax.set_ylim((0.0, 1.1))

    helpers.add_labels(fig, ax, xlabel='offline $p_T^{\mathrm{jet}}$ [GeV]', ylabel='Trigger Efficiency - Differential')
    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr, 'Shift: {:0.4f}'.format(eff_curve_shift)])
    helpers.add_legend(fig, ax)
    helpers.add_grid(fig, ax)
    pickle.dump(pl_eff_diff, file( write_file('plots/pickle/matched_differential_%s_trigger%d_region%s.pkl' % (filename_id, tJetEt_cut, i)), 'w+') )
    pl.savefig( write_file('plots/differential/%s_trigger%d_region%s.png' % (filename_id, tJetEt_cut, i)) )
    pl.close()

# ##### ACTUALLY START HERE #####

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
