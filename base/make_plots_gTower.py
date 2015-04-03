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

try:
  import cPickle as pickle
except:
  import pickle

from matplotlib.colors import LogNorm

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

helpers = PlotHelpers(dataSetStr=dataSetStr, seedCutStr=seedCutStr, noiseCutStr=noiseCutStr, towerThrStr=towerThrStr, labelsize=28)

tJet_exists = data['tJet.et'] > 0.

bins_efficiency = np.arange(0., 2000., 10.)
width_efficiency = np.array([x - bins_efficiency[i-1] for i, x in enumerate(bins_efficiency)][1:])
bins_multiplicity = np.arange(0.0, 1024.0, 32.0)

bins_rho = np.arange(0., 70., 0.5)
bins_vertices = np.arange(0., 100., 1.)

startTime_wall = time.time()
startTime_processor = time.clock()

try:
  print "gTower multiplicity"
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
  pl.close(fig)
  pass

valid_gJets = np.where(tJet_exists)

try:
  print "running out profiles for gTowers"
  fig, ax = pl.subplots(figsize=helpers.figsize)

  labels = ['$E_T^0$', '$E_T^1$', '$E_T^2$', '$E_T^3$']

  gTower_Et = np.zeros((3, np.sum(data[valid_gJets].size)))

  seed_Et, gTower_Et[0], gTower_Et[1], gTower_Et[2] = np.array([[seed.et, towers[0].et, towers[1].et, towers[2].et] for seed, towers in zip(data['tJet.seed'], data['tJet.towers'])])[valid_gJets].T
  all_gTower_Et = np.sort(np.vstack((seed_Et, gTower_Et)), axis=0)[::-1, :]

  ax.hist(seed_Et, bins=np.arange(0.0, 100.0, 2.5), stacked=True, fill=False, histtype='step', color=helpers.colors[0], label='$E_T^\mathrm{seed}$', linewidth=helpers.linewidth, weights=data['weight'][np.where(tJet_exists)])
  for datapts, color, label in zip(all_gTower_Et, helpers.colors[1:], labels):
    ax.hist(datapts, bins=np.arange(0.0, 100.0, 2.5), stacked=True, fill=False, histtype='step', color=color, label=label, linewidth=helpers.linewidth, weights=data['weight'][np.where(tJet_exists)])
  helpers.add_legend(fig, ax)
  helpers.add_grid(fig, ax)
  helpers.add_labels(fig, ax, xlabel='trigger jet\'s $E_T^\mathrm{gTower}$ [GeV]', ylabel=r'count/2.5 GeV bin')
  helpers.add_description(fig, ax, align='bl', strings=[helpers.dataSetStr, helpers.seedCutStr])
  ax.set_yscale('log')
  helpers.to_file(fig, ax, 'plots/gTowers/{}_Et.png'.format(filename_id))
  pl.close(fig)
except:
  print "Error for {}: could not make combined histograms of gTowers".format('gTowers')
  pl.close(fig)
  pass


profiles = {}
profiles_scaled = {}

for datapts, label in zip(all_gTower_Et, ['gTower 0', 'gTower 1', 'gTower 2', 'gTower 3']):
  try:
    fig, ax = pl.subplots(figsize=helpers.figsize)
    counts, edges_x, edges_y, im = ax.hist2d(datapts, data['oJet.nsj'][valid_gJets], bins=(np.arange(0., 100., 2.5), np.arange(1, 7, 1)), norm=LogNorm(), alpha=0.75, cmap=helpers.cmap, weights=data['weight'][np.where(tJet_exists)])
    points_y, mean_x, err_x = helpers.profile_x(edges_y, data['oJet.nsj'][valid_gJets], datapts)
    profiles[label.replace(' ', '_')] = {'x': mean_x, 'y': points_y, 'e': err_x}
    ax.scatter(mean_x, points_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)

    ticks = np.logspace(0, np.log10(np.max(counts)), 10)
    cbar = fig.colorbar(im, ticks=ticks, format=helpers.label_formatter)
    helpers.format_cbar(cbar)

    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr])
    helpers.add_labels(fig, ax, xlabel=r'$E_T^\mathrm{{{}}}$'.format(label), ylabel='offline jet\'s number of subjets (cut at 6)', title=label)

    ax.set_xlim((0.0, 100.0))
    ax.set_ylim((1.0, 6.0))
    helpers.add_grid(fig, ax)
    helpers.to_file(fig, ax, 'plots/gTowers/{}_{}.png'.format(filename_id, label.replace(' ', '_')))
    pl.close(fig)

    fig, ax = pl.subplots(figsize=helpers.figsize)
    counts, edges_x, edges_y, im = ax.hist2d(datapts/data['tJet.et'][valid_gJets], data['oJet.nsj'][valid_gJets], bins=(np.arange(0., 1., 0.05), np.arange(1, 7, 1)), norm=LogNorm(), alpha=0.75, cmap=helpers.cmap, weights=data['weight'][np.where(tJet_exists)])
    points_y, mean_x, err_x = helpers.profile_x(edges_y, data['oJet.nsj'][valid_gJets], datapts/data['tJet.et'][valid_gJets])
    profiles_scaled[label.replace(' ', '_')] = {'x': mean_x, 'y': points_y, 'e': err_x}
    ax.scatter(mean_x, points_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)

    ticks = np.logspace(0, np.log10(np.max(counts)), 10)
    cbar = fig.colorbar(im, ticks=ticks, format=helpers.label_formatter)
    helpers.format_cbar(cbar)

    helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr])
    helpers.add_labels(fig, ax, xlabel=r'$E_T^\mathrm{{{}}} / E_T^\mathrm{{tJet}}$'.format(label), ylabel='offline jet\'s number of subjets (cut at 6)', title=label)

    ax.set_xlim((0.0, 1.0))
    ax.set_ylim((1.0, 6.0))
    helpers.add_grid(fig, ax)
    helpers.to_file(fig, ax, 'plots/gTowers/{}_{}_over_tJetEt.png'.format(filename_id, label.replace(' ', '_')))
    pl.close(fig)
  except:
    print "Error for {}: could not make correlations with oJet.nsj".format(label)
    pl.close(fig)
    raise
    pass

# profile_x's merged
try:
  print "merging profiles for gTowers"
  fig, ax = pl.subplots(figsize=helpers.figsize)
  for l, c, m, yshift in zip(['gTower 0', 'gTower 1', 'gTower 2', 'gTower 3'], helpers.colors, helpers.markers, np.arange(-0.3, 0.5, 0.2)):
    datapts = profiles[l.replace(' ', '_')]
    ax.errorbar(datapts['x'], datapts['y']+yshift, ms=20, mfc=c, mec='k', marker=m, linestyle='--', linewidth=3, capsize=10, alpha=0.75, label=l, xerr=datapts['e'], ecolor='k')
  helpers.add_legend(fig, ax, numpoints=1)
  helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr])
  helpers.add_labels(fig, ax, xlabel=r'$E_T^\mathrm{gTower}$', ylabel='offline jet\'s number of subjets (cut at 6)')
  ax.set_xlim((0.0, 150.0))
  ax.set_ylim((1.0, 7.0))

  ax.set_yticklabels('')  # clear major tick labels on y-axis
  # add the labels on the minor ticks
  ax.set_yticks([1.5, 2.5, 3.5, 4.5, 5.5], minor=True)
  ax.set_yticklabels(['1', '2', '3', '4', '5'], minor=True)

  helpers.add_grid(fig, ax)
  helpers.to_file(fig, ax, 'plots/gTowers/{}_profiles.png'.format(filename_id))
  pl.close(fig)

  fig, ax = pl.subplots(figsize=helpers.figsize)
  for l, c, m, yshift in zip(['gTower 0', 'gTower 1', 'gTower 2', 'gTower 3'], helpers.colors, helpers.markers, np.arange(-0.3, 0.5, 0.2)):
    datapts = profiles_scaled[l.replace(' ', '_')]
    ax.errorbar(datapts['x'], datapts['y']+yshift, ms=20, mfc=c, mec='k', marker=m, linestyle='--', linewidth=3, capsize=10, alpha=0.75, label=l, xerr=datapts['e'], ecolor='k')
  helpers.add_legend(fig, ax, numpoints=1)
  helpers.add_description(fig, ax, align='br', strings=[helpers.dataSetStr, helpers.seedCutStr])
  helpers.add_labels(fig, ax, xlabel=r'$E_T^\mathrm{gTower} / E_T^\mathrm{tJet}$', ylabel='offline jet\'s number of subjets (cut at 6)')
  ax.set_xlim((0.0, 0.5))
  ax.set_ylim((1.0, 7.0))

  ax.set_yticklabels('')  # clear major tick labels on y-axis
  # add the labels on the minor ticks
  ax.set_yticks([1.5, 2.5, 3.5, 4.5, 5.5], minor=True)
  ax.set_yticklabels(['1', '2', '3', '4', '5'], minor=True)

  helpers.add_grid(fig, ax)
  helpers.to_file(fig, ax, 'plots/gTowers/{}_profiles_over_tJetEt.png'.format(filename_id))
  pl.close(fig)
except:
  print "Error for {}: could not merge gTower subjet profile_x's".format('gTowers')
  pl.close(fig)
  raise
  pass

# subjet studies comparing leading gTowers with offline subjet Pts
try:
  print "subjet studies comparing leading gTowers with leading offline subjet Pts"
  fig, ax = pl.subplots(figsize=helpers.figsize)
  # make a correlation of corrected Et versus trigger jet Et
  y = all_gTower_Et[0]
  x = np.array([el[0] for el in data['oJet.subjetsPt'][valid_gJets]])
  weights = data['weight'][np.where(tJet_exists)]
  bins_x = np.arange(0., 200., 4.)
  bins_y = np.arange(0., 100., 2.)
  xlim = (0., 200.)
  ylim = (0., 100.)
  label_y = r'leading $E_T^{\mathrm{gTower}}$/2 GeV [GeV]'
  label_x = r'leading $P_T^{\mathrm{oJet\ subjet}}$/4 GeV [GeV]'
  fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, xlim=xlim, ylim=ylim, profile_x=True, align='tl',
                           strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr], weights=weights)

  ax.plot((0., 200.), (args.seedEt_thresh, args.seedEt_thresh), 'k-', linewidth=helpers.linewidth, zorder=0)

  helpers.to_file(fig, ax, 'plots/subjets/{}_gTower_oJetSubjetPt.png'.format(filename_id))

  pl.close(fig)
except:
  print "Error for {}: could not make correlation plot between leading gTower and offline subjet Pt".format("first subjet")
  raise

# the new regions aren't working because there are not enough jets in them!!!
regions = {1: '', 2: '', 3: '', 4: '', '3a': '', '3b': '', '3c': '', '4a': '', '4b': '', '4c': ''}

for region in regions.keys():
  region_cut = helpers.region_cut(data['tJet.eta'], region)
  regions[region] = region_cut

NoCut = np.ones(np.sum(tJet_exists)).astype(bool)
PtCut = (250. <= data['oJet.pt'][np.where(tJet_exists)]) & (data['oJet.pt'][np.where(tJet_exists)] <= 500.)
MassCut = (100. <= data['oJet.m'][np.where(tJet_exists)]) & (data['oJet.m'][np.where(tJet_exists)] <= 220.)

PtCutLabel = "$250\ \mathrm{GeV}\ \leq P_T^\mathrm{oJet} \leq\ 500\ \mathrm{GeV}$"
MassCutLabel = "\n$100\ \mathrm{GeV}\ \leq m^\mathrm{oJet} \leq\ 220\ \mathrm{GeV}$"

offlineCuts = [NoCut, PtCut, PtCut & MassCut]

dataStorage = {}

print "Making subjet efficiency studies"
for gJetEt_cut in [200., 250., 300., 350.]:
  print "\t"*1, "gJetEt_cut = {}".format(gJetEt_cut)
  # add storage by gJetEt cut
  dataStorage['gJetEtCut_{:0.0f}'.format(gJetEt_cut)] = {}
  # create trigger cut
  triggerCut = (data['tJet.et'] >= gJetEt_cut)[np.where(tJet_exists)]
  gJetEtCutLabel = '$E_T^\mathrm{{gFEX\ jet}} \geq {:0.2f}\ \mathrm{{GeV}}$'.format(gJetEt_cut)

  fig, ax = pl.subplots(figsize=helpers.figsize)
  x, y, e, n, w = helpers.add_turnon(fig, ax, data=data['oJet.pt'][np.where(tJet_exists)], den_cut=MassCut, num_cut=(data['tJet.et'][np.where(tJet_exists)] > gJetEt_cut), label='no subt., no shift', p0=(1., 50., gJetEt_cut, 0.5))

  ax.set_xlim((0.0, 450.0))
  ax.set_ylim((0.0, 1.1))

  helpers.add_labels(fig, ax, xlabel='offline $p_T^{\mathrm{jet}}$ [GeV]', ylabel='Trigger Efficiency - Differential')
  helpers.add_description(fig, ax, align='tl', strings=[helpers.dataSetStr, 'iso., $\Delta R(\mathrm{gJet},\mathrm{oJet})\leq 1$', helpers.seedCutStr, helpers.noiseCutStr, MassCutLabel])
  helpers.add_legend(fig, ax, loc='lower left')
  helpers.add_grid(fig, ax)

  helpers.to_file(fig, ax, 'plots/subjets/{}_differential_gJetEt{:0.0f}{}.png'.format(filename_id, gJetEt_cut, '_100oJetM220'))
  pl.close(fig)

  for offlineCut, offlineLabel, filenameEnd in zip(offlineCuts, [r'', PtCutLabel, "{}{}".format(PtCutLabel, MassCutLabel)], ['', '_250oJetPt500', '_250oJetPt500_100oJetM220']):
    print "\t"*2, filenameEnd
    # add storage by offline cut
    dataStorage['gJetEtCut_{:0.0f}'.format(gJetEt_cut)][filenameEnd] = {}
    oJetsubjets_count = np.array(map(lambda x: np.sum(x > 20.0), data['oJet.subjetsPt'][np.where(tJet_exists)]))

    for gTowerEt_cut in [5., 10., 15., 20., 25.]:
      print "\t"*3, "gTowerEt_cut = {}".format(gTowerEt_cut)
      # add storage by gTowerEt cut
      dataStorage['gJetEtCut_{:0.0f}'.format(gJetEt_cut)][filenameEnd]['gTowerEtCut_{:0.0f}'.format(gTowerEt_cut)] = {}
      gTowersnsj_count = np.sum(all_gTower_Et > gTowerEt_cut, axis=0)
      # make a label for gTowerEt cut
      gTowerEtThrLabel = r'$N\left(E_T^\mathrm{{gFEX\ tower}} >\ {:0.2f}\ \mathrm{{GeV}}\right) \geq X$'.format(gTowerEt_cut)
      fig, ax = pl.subplots(figsize=helpers.figsize)

      for oJetnsj_cut, color in zip([0, 1, 2, 3, 4], helpers.colors):
        print "\t"*4, "oJetnsj_cut = {}".format(oJetnsj_cut)
        # add storage by oJetnsj_cut
        myStore = dataStorage['gJetEtCut_{:0.0f}'.format(gJetEt_cut)][filenameEnd]['gTowerEtCut_{:0.0f}'.format(gTowerEt_cut)]['oJetnsjCut_{:d}'.format(oJetnsj_cut)] = {}
        # find oJets with exactly N subjets with Pt >= 20.0 GeV
        cut = (oJetsubjets_count >= oJetnsj_cut)
        # make a label for oJet nsj cut
        oJetnsjCutLabel = r'$N(P_T^\mathrm{{oJet\ subjet}} >\ 20 \ \mathrm{{GeV}}) \geq {:d}$'.format(oJetnsj_cut)
        # make histograms, sum, and plot
        hist, bins = np.histogram(gTowersnsj_count[np.where(cut & offlineCut & triggerCut)], bins=np.arange(1, 6, 1), weights=data['weight'][np.where(tJet_exists)][np.where(cut & offlineCut & triggerCut)])
        cumsum = np.cumsum(hist[::-1].astype(float))[::-1]/np.sum(hist)
        ax.plot(bins[:-1], cumsum, linestyle='steps-mid', alpha=0.75, color=color, marker='x', ms=20, mew=10, linewidth=0, label=oJetnsjCutLabel)
        # store the data
        myStore['bins'] = bins[:-1]
        myStore['vals'] = cumsum

      ax.set_xlim((-1.0, 6.0))
      ax.set_ylim((0.0, 2.0))
      ax.set_xticks(np.arange(1, 5, 1))
      ax.set_yticks(np.linspace(0., 1., 6))
      helpers.add_legend(fig, ax, numpoints=1)
      helpers.add_labels(fig, ax, xlabel=gTowerEtThrLabel, ylabel=plotConfigs.yaxis_efficiencyLabel, title=gJetEtCutLabel)
      helpers.add_description(fig, ax, align='tl', strings=[helpers.dataSetStr, 'iso., $\Delta R(\mathrm{gJet},\mathrm{oJet})\leq 1$', helpers.seedCutStr, helpers.noiseCutStr, offlineLabel])
      helpers.add_grid(fig, ax)
      helpers.to_file(fig, ax, 'plots/subjets/{}_XgTower_gJetEt{:0.0f}_gTowerEt{:0.0f}{}.png'.format(filename_id, gJetEt_cut, gTowerEt_cut, filenameEnd))
      pl.close(fig)

      print "\t"*4, "-----"
      # make the denominator
      hist_den, bins_den = np.histogram(oJetsubjets_count[np.where(offlineCut)], bins=np.arange(1, 6, 1), weights=data['weight'][np.where(tJet_exists)][np.where(offlineCut)])
      # integral, inclusive
      hist_den = np.cumsum(hist_den[::-1])[::-1]

      fig, ax = pl.subplots(figsize=helpers.figsize)
      for gJetnsj_cut, color in zip([1, 2, 3, 4], helpers.colors):
        print "\t"*4, "gJetnsj_cut = {}".format(gJetnsj_cut)
        # add storage by oJetnsj_cut
        myStore = dataStorage['gJetEtCut_{:0.0f}'.format(gJetEt_cut)][filenameEnd]['gTowerEtCut_{:0.0f}'.format(gTowerEt_cut)]['gJetnsjCut_{:d}'.format(gJetnsj_cut)] = {}

        # find gJets with >= N gTowers with Et > gTowerEt_cut
        cut = (gTowersnsj_count >= gJetnsj_cut)

        # make a label for gJet nsj cut
        gJetnsjCutLabel = r'$N(E_T^\mathrm{{gJet\ gTower}} > {:0.2f} \ \mathrm{{GeV}}) \geq {:d}$'.format(gTowerEt_cut, gJetnsj_cut)

        hist_num, bins_num = np.histogram(oJetsubjets_count[np.where(cut & offlineCut & triggerCut)], bins=np.arange(1, 6, 1), weights=data['weight'][np.where(tJet_exists)][np.where(cut & offlineCut & triggerCut)])
        # integral, inclusive
        hist_num = np.cumsum(hist_num[::-1])[::-1]

        vals = hist_num.astype(float)/hist_den.astype(float)
        ax.plot(bins_den[:-1], vals, linestyle='steps-mid', alpha=0.75, color=color, marker='x', ms=20, mew=10, linewidth=0, label=gJetnsjCutLabel)

        # store the data
        myStore['bins'] = bins_den[:-1]
        myStore['vals'] = vals

      ax.set_xlim((0.0, 6.0))
      ax.set_ylim((0.0, 2.0))
      ax.set_xticks(np.arange(1, 6, 1))
      ax.set_yticks(np.linspace(0., 1., 6))
      helpers.add_legend(fig, ax, numpoints=1)
      helpers.add_labels(fig, ax, xlabel=r'$N(P_T^\mathrm{{oJet\ subjet}} >\ 20 \ \mathrm{{GeV}}) \geq X$', ylabel='Efficiency', title=gJetEtCutLabel)
      helpers.add_description(fig, ax, align='tl', strings=[helpers.dataSetStr, 'isolated, $\Delta R(\mathrm{gJet},\mathrm{oJet})\leq 1$', helpers.seedCutStr, helpers.noiseCutStr, offlineLabel])
      helpers.add_grid(fig, ax)
      helpers.to_file(fig, ax, 'plots/subjets/{}_XoJet_gJetEt{:0.0f}_gTowerEt{:0.0f}{}.png'.format(filename_id, gJetEt_cut, gTowerEt_cut, filenameEnd))
      pl.close(fig)

pickle.dump(dataStorage, file(helpers.write_file('plots/subjets/{}.pkl'.format(filename_id)), 'w+'))

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: {:0.2f} s \n\t Clock Time: {:0.2f} s".format((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
