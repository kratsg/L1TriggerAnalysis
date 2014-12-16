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

startTime_wall = time.time()
startTime_processor = time.clock()

# standard kinematic plots
fig, ax = pl.subplots(figsize=helpers.figsize)
n, bins, unweightedPatches = ax.hist(data['oJet.pt'], bins=np.arange(0, 500, 2), label='unweighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='r')
axt = ax.twinx()
n, bins, weightedPatches = axt.hist(data['oJet.pt'], weights=data['weight'], bins=np.arange(0, 500, 2), label='weighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='b')
# http://matplotlib.org/examples/api/two_scales.html
for tl in ax.get_yticklabels():
  tl.set_color('r')
for tl in axt.get_yticklabels():
  tl.set_color('b')

# make bigger
ax.xaxis.set_tick_params(width=2, which='both')
ax.yaxis.set_tick_params(width=2, which='both')
# make bigger
axt.xaxis.set_tick_params(width=2, which='both')
axt.yaxis.set_tick_params(width=2, which='both')

# lengthen
ax.yaxis.set_tick_params(length=10, which='major')
ax.yaxis.set_tick_params(length=7, which='minor')
axt.yaxis.set_tick_params(length=10, which='major')
axt.yaxis.set_tick_params(length=7, which='minor')


# http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
patches = [unweightedPatches[0], weightedPatches[0]]
# http://matplotlib.org/examples/pylab_examples/legend_auto.html
labels = [p.get_label() for p in patches]
ax.set_yscale('log', nonposy='clip')
axt.set_yscale('log', nonposy='clip')
legend = ax.legend(patches, labels, fancybox=True, framealpha=0.75, fontsize=helpers.labelsize)
legend.get_frame().set_facecolor(helpers.light_grey)
legend.get_frame().set_linewidth(0.0)
helpers.add_labels(fig, ax, xlabel=r'$p_T$ [GeV]', ylabel=r'unweighted count')
helpers.add_labels(fig, axt, ylabel=r'weighted count')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, "plots/offline_jet_kinematics/{}_oJet_Pt.png".format(filename_id))
pl.close(fig)

fig, ax = pl.subplots(figsize=helpers.figsize)
n, bins, unweightedPatches = ax.hist(data['oJet.eta'], bins=np.arange(-4.9, 4.9, 0.2), label='unweighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='r')
axt = ax.twinx()
n, bins, weightedPatches = axt.hist(data['oJet.eta'], weights=data['weight'], bins=np.arange(-4.9, 4.9, 0.2), label='weighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='b')
# http://matplotlib.org/examples/api/two_scales.html
for tl in ax.get_yticklabels():
  tl.set_color('r')
for tl in axt.get_yticklabels():
  tl.set_color('b')

# make bigger
ax.xaxis.set_tick_params(width=2, which='both')
ax.yaxis.set_tick_params(width=2, which='both')
# make bigger
axt.xaxis.set_tick_params(width=2, which='both')
axt.yaxis.set_tick_params(width=2, which='both')

# lengthen
ax.yaxis.set_tick_params(length=10, which='major')
ax.yaxis.set_tick_params(length=7, which='minor')
axt.yaxis.set_tick_params(length=10, which='major')
axt.yaxis.set_tick_params(length=7, which='minor')


# http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
patches = [unweightedPatches[0], weightedPatches[0]]
# http://matplotlib.org/examples/pylab_examples/legend_auto.html
labels = [p.get_label() for p in patches]
ax.set_yscale('log', nonposy='clip')
axt.set_yscale('log', nonposy='clip')
legend = ax.legend(patches, labels, fancybox=True, framealpha=0.75, fontsize=helpers.labelsize)
legend.get_frame().set_facecolor(helpers.light_grey)
legend.get_frame().set_linewidth(0.0)
helpers.add_labels(fig, ax, xlabel=r'$\eta^\mathrm{oJet}$', ylabel=r'unweighted count')
helpers.add_labels(fig, axt, ylabel=r'weighted count')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, "plots/offline_jet_kinematics/{}_oJet_Eta.png".format(filename_id))
pl.close(fig)

fig, ax = pl.subplots(figsize=helpers.figsize)
n, bins, unweightedPatches = ax.hist(data['oJet.phi'], bins=np.arange(-3.2, 3.2, 0.2), label='unweighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='r')
axt = ax.twinx()
n, bins, weightedPatches = axt.hist(data['oJet.phi'], weights=data['weight'], bins=np.arange(-3.2, 3.2, 0.2), label='weighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='b')
# http://matplotlib.org/examples/api/two_scales.html
for tl in ax.get_yticklabels():
  tl.set_color('r')
for tl in axt.get_yticklabels():
  tl.set_color('b')

# make bigger
ax.xaxis.set_tick_params(width=2, which='both')
ax.yaxis.set_tick_params(width=2, which='both')
# make bigger
axt.xaxis.set_tick_params(width=2, which='both')
axt.yaxis.set_tick_params(width=2, which='both')

# lengthen
ax.yaxis.set_tick_params(length=10, which='major')
ax.yaxis.set_tick_params(length=7, which='minor')
axt.yaxis.set_tick_params(length=10, which='major')
axt.yaxis.set_tick_params(length=7, which='minor')


# http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
patches = [unweightedPatches[0], weightedPatches[0]]
# http://matplotlib.org/examples/pylab_examples/legend_auto.html
labels = [p.get_label() for p in patches]
ax.set_yscale('log', nonposy='clip')
axt.set_yscale('log', nonposy='clip')
legend = ax.legend(patches, labels, fancybox=True, framealpha=0.75, fontsize=helpers.labelsize)
legend.get_frame().set_facecolor(helpers.light_grey)
legend.get_frame().set_linewidth(0.0)
helpers.add_labels(fig, ax, xlabel=r'$\phi^\mathrm{oJet}$', ylabel=r'unweighted count')
helpers.add_labels(fig, axt, ylabel=r'weighted count')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, "plots/offline_jet_kinematics/{}_oJet_Phi.png".format(filename_id))
pl.close(fig)

fig, ax = pl.subplots(figsize=helpers.figsize)
n, bins, unweightedPatches = ax.hist(data['oJet.m'], bins=np.arange(0, 500, 2), label='unweighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='r')
axt = ax.twinx()
n, bins, weightedPatches = axt.hist(data['oJet.m'], weights=data['weight'], bins=np.arange(0, 500, 2), label='weighted', stacked=True, fill=False, histtype='step', alpha=0.75, color='b')
# http://matplotlib.org/examples/api/two_scales.html
for tl in ax.get_yticklabels():
  tl.set_color('r')
for tl in axt.get_yticklabels():
  tl.set_color('b')

# make bigger
ax.xaxis.set_tick_params(width=2, which='both')
ax.yaxis.set_tick_params(width=2, which='both')
# make bigger
axt.xaxis.set_tick_params(width=2, which='both')
axt.yaxis.set_tick_params(width=2, which='both')

# lengthen
ax.yaxis.set_tick_params(length=10, which='major')
ax.yaxis.set_tick_params(length=7, which='minor')
axt.yaxis.set_tick_params(length=10, which='major')
axt.yaxis.set_tick_params(length=7, which='minor')


# http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
patches = [unweightedPatches[0], weightedPatches[0]]
# http://matplotlib.org/examples/pylab_examples/legend_auto.html
labels = [p.get_label() for p in patches]
ax.set_yscale('log', nonposy='clip')
axt.set_yscale('log', nonposy='clip')
legend = ax.legend(patches, labels, fancybox=True, framealpha=0.75, fontsize=helpers.labelsize)
legend.get_frame().set_facecolor(helpers.light_grey)
legend.get_frame().set_linewidth(0.0)
helpers.add_labels(fig, ax, xlabel=r'$m^\mathrm{oJet}$ [GeV]', ylabel=r'unweighted count')
helpers.add_labels(fig, axt, ylabel=r'weighted count')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, "plots/offline_jet_kinematics/{}_oJet_M.png".format(filename_id))
pl.close(fig)

# weights
fig, ax = pl.subplots(figsize=helpers.figsize)
ax.hist(data['weight'], bins=100, fill=False, alpha=0.75, linewidth=helpers.linewidth)
helpers.add_labels(fig, ax, xlabel=r'$\omega_{event} * \sigma * \varepsilon_{filter} / N_{events}$', ylabel='counts')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.add_grid(fig, ax)
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_weights.png' % (filename_id))
pl.close(fig)

# spatial correlations
x = data['oJet.eta']
y = data['oJet.phi']
bins_x = np.arange(-2.7, 2.7, 0.2)
bins_y = np.arange(-3.2, 3.2, 0.2)
label_x = r'$\eta^\mathrm{oJet}$'
label_y = r'$\phi^\mathrm{oJet}$'
fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, profile_x=True, profile_y=True, align='tl',
                         strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_angular_positions.png' % (filename_id))
pl.close(fig)

# y projection slices
fig, ax = pl.subplots(figsize=helpers.figsize)
selections = {'left':   np.where((-1.9 < data['oJet.eta']) & (data['oJet.eta'] < -1.5)),
              'middle': np.where((-0.2 < data['oJet.eta']) & (data['oJet.eta'] < 0.2)),
              'right':  np.where((1.5 < data['oJet.eta']) & (data['oJet.eta'] < 1.9))}
ax.hist(data['oJet.phi'][selections['left']], bins=100, weights=data['weight'][selections['left']], stacked=True, fill=False, histtype='step', linewidth=helpers.linewidth, alpha=0.75, color='b', label='$-1.9\ < \eta^\mathrm{oJet} <\ -1.5$')
ax.hist(data['oJet.phi'][selections['middle']], bins=100, weights=data['weight'][selections['middle']], stacked=True, fill=False, histtype='step', linewidth=helpers.linewidth, alpha=0.75, color='c', label='$-0.2\ < \eta^\mathrm{oJet} <\ 0.2$')
ax.hist(data['oJet.phi'][selections['right']], bins=100, weights=data['weight'][selections['right']], stacked=True, fill=False, histtype='step', linewidth=helpers.linewidth, alpha=0.75, color='r', label='$1.5\ < \eta^\mathrm{oJet} <\ 1.9$')
helpers.add_legend(fig, ax)
helpers.add_labels(fig, ax, xlabel=r'$\eta^\mathrm{oJet}$', ylabel='weighted counts', title='Y-Axis Projections of $\phi^\mathrm{oJet}$')
helpers.add_description(fig, ax, align='cr', strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
ax.set_yscale('log', nonposy='clip')
helpers.add_grid(fig, ax)
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_offline_jet_phi_projection_y.png' % (filename_id))
pl.close(fig)

# eta-pt correlations
x = data['oJet.eta']
y = data['oJet.pt']
bins_x = np.arange(-2.7, 2.7, 0.2)
bins_y = np.arange(0., 375., 10.)
label_x = r'$\eta^\mathrm{oJet}$'
label_y = r'$p_T^\mathrm{oJet}$'
fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, profile_x=True, profile_y=True, align='tl',
                         strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_eta_pt_corr.png' % (filename_id))
pl.close(fig)

# pt-phi correlations
x = data['oJet.pt']
y = data['oJet.phi']
bins_x = np.arange(0., 375., 10.)
bins_y = np.arange(-3.2, 3.2, 0.2)
label_x = r'$p_T^\mathrm{oJet}$'
label_y = r'$\phi^\mathrm{oJet}$'
fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, profile_x=True, profile_y=True, align='tl',
                         strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_pt_phi_corr.png' % (filename_id))
pl.close(fig)

# eta-m correlations
x = data['oJet.eta']
y = data['oJet.m']
bins_x = np.arange(-2.7, 2.7, 0.2)
bins_y = np.arange(0., 175., 10.)
label_x = r'$\eta^\mathrm{oJet}$'
label_y = r'$m^\mathrm{oJet}$'
fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, profile_x=True, profile_y=True, align='tl',
                         strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_eta_m_corr.png' % (filename_id))
pl.close(fig)

# m-phi correlations
x = data['oJet.m']
y = data['oJet.phi']
bins_x = np.arange(0., 175., 10.)
bins_y = np.arange(-3.2, 3.2, 0.2)
label_x = r'$m^\mathrm{oJet}$'
label_y = r'$\phi^\mathrm{oJet}$'
fig, ax = helpers.corr2d(x, y, bins_x, bins_y, label_x, label_y, profile_x=True, profile_y=True, align='tl',
                         strings=[helpers.dataSetStr, helpers.seedCutStr, helpers.noiseCutStr, helpers.towerThrStr])
helpers.to_file(fig, ax, 'plots/offline_jet_kinematics/%s_m_phi_corr.png' % (filename_id))
pl.close(fig)
endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))
