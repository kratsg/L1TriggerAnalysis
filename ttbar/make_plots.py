import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import argparse
import time

import sys
from plots_wrapper import PlotHelpers

# locally load the config file
sys.path.insert(0, '.')
import plot_configs as plotConfigs
sys.path.remove('.')

# for getting FWHM
from scipy.interpolate import interp1d

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

filename_id = "seed{:0.0f}_noise{:0.0f}_signal{:0.0f}_digitization{:0.0f}".format(args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)
filename = "data/seed{:0.0f}/leading_jets_{}.pkl".format(args.seedEt_thresh, filename_id)
data = pickle.load(file(filename))

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished reading in data:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

dataSetStr  = plotConfigs.dataSetStr
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

bins_rho = np.arange(0., 70., 0.5)
bins_vertices = np.arange(0., 100., 1.)

startTime_wall = time.time()
startTime_processor = time.clock()

# this will be used to make resolution plots below line 200
resolution = (data['tJet.et']/data['oJet.pt']) - 1.0
resolution_subtracted = (tJetEt_subtracted/data['oJet.pt']) - 1.0

for i, cut in regions.iteritems():
  region = 'region_%s' % i
  region_cut = np.where(cut)
  print region
  try:
    print "\t", "correction Et versus tJetEt"
    # make a correlation of corrected Et versus trigger jet Et
    x = data['tJet.et'][region_cut]
    y = tJetEt_correction[region_cut]
    bins_x = bins_y = np.arange(0., 1500., 10.)
    xlim = (0., 1000.)
    ylim = (0., 200.)
    label_x = r'$E_T^{\mathrm{gJet}}$ [GeV]'
    label_y = r'$\rho*A^\mathrm{gJet}$ [GeV]'

    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_y=True, align='tr',
                             strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/trigger_jet_Et_correction/{}_trigger_jet_Et_correction_region{}.png'.format(filename_id, i))

    pl.close(fig)
  except:
    print "\t", "Error for {}: could not make correlation of corrected Et versus trigger jet Et".format(region)
    pl.close(fig)
    pass

  try:
    print "\t", "oJetPt versus tJetEt"
    x = data['oJet.pt'][np.where(tJet_exists & cut)]
    y = data['tJet.et'][np.where(tJet_exists & cut)]
    bins_x = bins_y = np.arange(0., 1500., 10.)
    xlim = ylim = (0., 1000.)
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'trigger $E_T^{\mathrm{jet}}$ [GeV]'
    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_y=True, align='br', title=r'no $\rho$ subtraction',
                             strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/jet_energy_correlation/{}_region{}.png'.format(filename_id, i))

    pl.close(fig)
  except:
    print "\t", "Error for {}: could not make jet energy correlation".format(region)
    pl.close(fig)
    pass

  try:
    print "\t", "oJetPt versus corrected tJetEt"
    x = data['oJet.pt'][np.where(tJet_exists_subtracted & cut)]
    y = tJetEt_subtracted[np.where(tJet_exists_subtracted & cut)]
    bins_x = bins_y = np.arange(0., 1500., 10.)
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'trigger $E_T^{\mathrm{jet}}$ [GeV]'
    xlim = ylim = (0., 1000.)
    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_y=True, align='br', title=r'with $\rho$ subtraction',
                             strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/jet_energy_correlation/{}_noPileup_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "\t", "Error for {}: could not make jet energy correlation (no Pileup)".format(region)
    pl.close(fig)
    pass

  # we want to make resolution plots now
  # note that we ignore trigger jets that are negative or less than zero

  print "\t", "making resolution plots"
  try:
    print "\t"*2, "oJetPt versus resolution"
    x = data['oJet.pt'][np.where(tJet_exists_subtracted & cut)]
    y = resolution[np.where(tJet_exists_subtracted & cut)]
    bins_x = bins_y = 100
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$'
    xlim = (0., 1000.)
    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, profile_y=True, align='tr', title=r'no $\rho$ subtraction',
                             strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "\t"*2, "Error for {}: could not make resolution correlation".format(region)
    pl.close(fig)
    pass

  try:
    print "\t"*2, "oJetPt versus corrected resolution"
    x = data['oJet.pt'][np.where(tJet_exists_subtracted & cut)]
    y = resolution_subtracted[np.where(tJet_exists_subtracted & cut)]
    bins_x = bins_y = 100
    label_x = r'offline $p_T^{\mathrm{jet}}$ [GeV]'
    label_y = r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$'
    xlim = (0., 1000.)
    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, profile_y=True, align='tr', title=r'with $\rho$ subtraction',
                             strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_withRhoSubtraction_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "\t", "Error for {}: could not make resolution correlation (with Rho subtraction)".format(region)
    pl.close(fig)
    pass

  print "\t", "making y-projection of resolution plots"
  try:
    print "\t"*2, "y-projection slices of resolution"
    # y projection slices
    pl_res_proj = {}
    fig, ax = pl.subplots(figsize=helpers.figsize)
    for oJetPt_cuts in [(170., 180.), (200., 220.), (300., 350.)]:
      oJetPt_cut = helpers.btwn(data['oJet.pt'], oJetPt_cuts[0], oJetPt_cuts[1])
      hist, bins = np.histogram(resolution[np.where(cut & oJetPt_cut & tJet_exists_subtracted)], bins=100, density=True)
      fwhm = helpers.FWHM(bins, hist)
      ax.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='b', label=r'${:0.0f}\ \mathrm{{GeV}} < p_T^\mathrm{{oJet}} <\ {:0.0f}\ \mathrm{{GeV}}$\nFWHM = {:0.4f}'.format(oJetPt_cuts[0], oJetPt_cuts[1], fwhm), linewidth=helpers.linewidth)

      pl_res_proj['{:0.0f}to{:0.0f}'.format(oJetPt_cuts[0], oJetPt_cuts[1])] = resolution[np.where(cut & oJetPt_cut & tJet_exists_subtracted)]

    helpers.add_legend(fig, ax)
    helpers.add_labels(fig, ax, xlabel=r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', ylabel='normalized counts', title='Y-Axis Projections of Resolution')
    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
    ax.set_xlim((-1.0, 1.0))
    helpers.add_grid(fig, ax)
    pickle.dump(pl_res_proj, file(helpers.write_file('plots/pickle/{}_resolution_PtOffline_projection_region{}.pkl'.format(filename_id, i)), 'w+'))
    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_projection_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "\t"*2, "Error for {}: could not make resolution projection".format(region)
    pl.close(fig)
    pass

  try:
    print "\t"*2, "y-projection slices of corrected resolution"
    pl_res_proj = {}
    fig, ax = pl.subplots(figsize=helpers.figsize)
    for oJetPt_cuts in [(170., 180.), (200., 220.), (300., 350.)]:
      oJetPt_cut = helpers.btwn(data['oJet.pt'], oJetPt_cuts[0], oJetPt_cuts[1])
      hist, bins = np.histogram(resolution_subtracted[np.where(cut & oJetPt_cut & tJet_exists_subtracted)], bins=100, density=True)
      fwhm = helpers.FWHM(bins, hist)
      ax.plot(bins[:-1], hist, linestyle='steps-post', alpha=0.75, color='b', label=r'${:0.0f}\ \mathrm{{GeV}} < p_T^\mathrm{{oJet}} <\ {:0.0f}\ \mathrm{{GeV}}$\nFWHM = {:0.4f}'.format(oJetPt_cuts[0], oJetPt_cuts[1], fwhm), linewidth=helpers.linewidth)

      pl_res_proj['{:0.0f}to{:0.0f}'.format(oJetPt_cuts[0], oJetPt_cuts[1])] = resolution_subtracted[np.where(cut & oJetPt_cut & tJet_exists_subtracted)]

    helpers.add_legend(fig, ax)
    helpers.add_labels(fig, ax, xlabel=r'resolution $\frac{E_T^\mathrm{gFEX} - p_T^\mathrm{offline}}{p_T^\mathrm{offline}}$', ylabel='normalized counts', title=r'Y-Axis Projections of Resolution, with $\rho$ subtraction')
    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
    ax.set_xlim((-1.0, 1.0))
    helpers.add_grid(fig, ax)
    pickle.dump(pl_res_proj, file(helpers.write_file('plots/pickle/{}_resolution_PtOffline_withRhoSubtraction_projection_region{}.pkl'.format(filename_id, i)), 'w+'))
    helpers.to_file(fig, ax, 'plots/resolution/{}_resolution_PtOffline_withRhoSubtraction_projection_region{}.png'.format(filename_id, i))
    pl.close(fig)
  except:
    print "\t"*2, "Error for {}: could not make resolution projection (no Pileup)".format(region)
    pl.close(fig)
    pass

  points_x, mean_y, err_y = helpers.profile_y(np.arange(0., 1500., 10.), data['tJet.et'][region_cut], tJetEt_correction[region_cut])
  f = interp1d(points_x, mean_y, bounds_error=False, fill_value=0., kind='cubic')

  print "\t", "efficiency curves using tJet Et cut in numerator"
  for tJetEt_cut in [0., 140.]:
    # ADD IN np.interp(tJetEt_cut, trigger_Et, rho*A) and shift subtracted eff curves
    eff_curve_shift = f([tJetEt_cut])[0]
    print "\t"*2, "tJetEt_cut = {}, eff_curve_shift = {}".format(tJetEt_cut, eff_curve_shift)
    # eff_curve_shift = np.interp(tJetEt_cut, trigger_jet_Et[region], trigger_jet_Et_correction[region])
    tJetEtThrStr = r'$E_T^\mathrm{{gFEX\ jet}} >\ {:0.2f}\ \mathrm{{GeV}}$'.format(tJetEt_cut)

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
    pickle.dump(pl_eff_diff, file(helpers.write_file('plots/pickle/matched_differential_%s_trigger%d_region%s.pkl' % (filename_id, tJetEt_cut, i)), 'w+'))
    helpers.to_file(fig, ax, 'plots/differential/{}_trigger{:0.0f}_region{}.png'.format(filename_id, tJetEt_cut, i))
    pl.close(fig)

# buildin plots for rho
print "rho histograms"
fig, ax = pl.subplots(figsize=helpers.figsize)
ax.hist(data['offline_rho'], bins=bins_rho, stacked=True, fill=False, histtype='step', color='b', label=r'offline Kt4LCTopo', linewidth=helpers.linewidth)
ax.hist(data['gFEX_rho_all'], bins=bins_rho, stacked=True, fill=False, histtype='step', color='r', label=helpers.region_legend('all'), linewidth=helpers.linewidth)
ax.hist(data['gFEX_rho_1'],   bins=bins_rho, stacked=True, fill=False, histtype='step', color='c', label=helpers.region_legend(1), linewidth=helpers.linewidth)
ax.hist(data['gFEX_rho_2'],   bins=bins_rho, stacked=True, fill=False, histtype='step', color='m', label=helpers.region_legend(2), linewidth=helpers.linewidth)
ax.hist(data['gFEX_rho_3'],   bins=bins_rho, stacked=True, fill=False, histtype='step', color='y', label=helpers.region_legend(3), linewidth=helpers.linewidth)
ax.hist(data['gFEX_rho_4'],   bins=bins_rho, stacked=True, fill=False, histtype='step', color='k', label=helpers.region_legend(4), linewidth=helpers.linewidth)
helpers.add_legend(fig, ax)
helpers.add_labels(fig, ax, xlabel=r'$\rho$ [GeV]', ylabel=r'count')
helpers.add_grid(fig, ax)
helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.towerThrStr])
helpers.to_file(fig, ax, "plots/pileup/{}_rho.png".format(filename_id))
pl.close(fig)

print "delta rho histograms"
fig, ax = pl.subplots(figsize=helpers.figsize)
diff21 = data['gFEX_rho_2'] - data['gFEX_rho_1']
diff43 = data['gFEX_rho_4'] - data['gFEX_rho_3']
xlimit = np.ceil(np.max([np.max(diff21), np.max(diff43)])/10)*10

hist, bins = np.histogram(data['gFEX_rho_2'] - data['gFEX_rho_1'], density=True, bins=np.arange(-xlimit, xlimit+1., 1.))
fwhm = helpers.FWHM(bins, hist)
ax.plot(bins[:-1], hist, linestyle='steps-post', color='b', alpha=0.75, label='region 2 - region 1\nFWHM: %0.4f' % fwhm, linewidth=helpers.linewidth)

hist, bins = np.histogram(data['gFEX_rho_4'] - data['gFEX_rho_3'], density=True, bins=np.arange(-xlimit, xlimit+1., 1.))
fwhm = helpers.FWHM(bins, hist)
ax.plot(bins[:-1], hist, linestyle='steps-post', color='r', alpha=0.75, label='region 4 - region 3\nFWHM: %0.4f' % fwhm, linewidth=helpers.linewidth)

helpers.add_legend(fig, ax)
helpers.add_labels(fig, ax, xlabel=r'$\Delta\rho$ [GeV]', ylabel=r'normalized counts')
helpers.add_grid(fig, ax)
helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.towerThrStr])

helpers.to_file(fig, ax, "plots/pileup/{}_deltaRho.png".format(filename_id))
ax.set_yscale('log', nonposy='clip')
helpers.to_file(fig, ax, "plots/pileup/{}_deltaRho_logy.png".format(filename_id))
pl.close(fig)

print "vxpN versus oRho"
x = data['vxp_n']
y = data['offline_rho']
bins_x = bins_vertices
bins_y = bins_rho
label_x = helpers.labels['vxp.number']
label_y = helpers.labels['rho.offline']
xlim = ylim = (0., 1000.)
fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, align='tr',
                         strings=[helpers.dataSetStr])

helpers.to_file(fig, ax, 'plots/pileup/{}_offlineRho.png'.format(filename_id))
pl.close(fig)

print "vxpN(1) / oRho(2) versus"
for col, legend in zip(['gFEX_rho_all', 'gFEX_rho_1', 'gFEX_rho_2', 'gFEX_rho_3', 'gFEX_rho_4'],
                       ['all', 1, 2, 3, 4]):

  try:
    print "(1)\t{}".format(col)
    x = data['vxp_n']
    y = data[col]
    bins_x = bins_vertices
    bins_y = bins_rho
    label_x = helpers.labels['vxp.number']
    label_y = '{}, {}'.format(helpers.labels['rho.offline'], helpers.region_legend(legend))
    xlim = ylim = (0., 1000.)
    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, align='tr',
                             strings=[helpers.dataSetStr, helpers.towerThrStr])

    helpers.to_file(fig, ax, 'plots/pileup/{}_{}.png'.format(filename_id, col))
    pl.close(fig)
  except:
    print "Error for {}: could not make correlation between gFEX rho and primary vertices".format(col)
    pl.close(fig)
    pass

  try:
    print "(2)\t{}".format(col)
    x = data['offline_rho']
    y = data[col]
    bins_x = bins_rho
    bins_y = bins_rho
    label_x = helpers.labels['rho.offline']
    label_y = '{}, {}'.format(helpers.labels['rho.gFEX'], helpers.region_legend(legend))
    xlim = ylim = (0., 1000.)
    fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, align='br',
                             strings=[helpers.dataSetStr, helpers.towerThrStr])

    helpers.add_atlas(fig, ax, level=1)

    helpers.to_file(fig, ax, 'plots/pileup/{}_{}_correlation.png'.format(filename_id, col))
    helpers.to_file(fig, ax, 'plots/pileup/{}_{}_correlation.pdf'.format(filename_id, col))
    helpers.to_file(fig, ax, 'plots/pileup/{}_{}_correlation.eps'.format(filename_id, col))
    pl.close(fig)
  except:
    print "Error for {}: could not make correlation between offline rho and gFEX rho".format(col)
    pl.close(fig)
    pass

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: {:0.2f} s \n\t Clock Time: {:0.2f} s".format((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
