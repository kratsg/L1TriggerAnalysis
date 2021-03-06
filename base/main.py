import argparse

parser = argparse.ArgumentParser(description='Process TTbar events for Level 1 Trigger.')
parser.add_argument('--processNum',   type=int,   required=False, dest='process_num',     help='a process number for tagging files', default=0)
parser.add_argument('--file',         type=str,   required=True,  dest='filename',        help='a filename to read in')
parser.add_argument('--start',        type=int,   required=True,  dest='event_start',     help='start index of event')
parser.add_argument('--numEvents',    type=int,   required=True,  dest='num_events',      help='number of events to process')
parser.add_argument('--stepSize',     type=int,   required=False, dest='step_size',       help='chunking of events', default=10)
parser.add_argument('--tree',         type=str,   required=False, dest='treename',        help='treename', default='mytree')
parser.add_argument('--seedEt',       type=float, required=False, dest='seedEt_thresh',   help='seed Et Threshold (GeV)', default=6.)
parser.add_argument('--noiseFilter',  type=float, required=False, dest='noise_filter',    help='gTower noise filter (GeV)', default=1.0)
parser.add_argument('--towerThresh',  type=float, required=False, dest='tower_thresh',    help='gTower rho Et thresh (GeV)', default=3.)
parser.add_argument('--digitization', type=float, required=False, dest='digitization',    help='digitization for gTower Et (MeV)', default=0.)
# parse the arguments, throw errors if missing any
args = parser.parse_args()

print args

# load up everything else
from atlas_jets import *
import numpy as np
import root_numpy as rnp
from ROOT import TFile, TLorentzVector

from itertools import chain

import os
try:
  import cPickle as pickle
except:
  import pickle

import time


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def write_file(f):
  ensure_dir(f)
  return f


def distance(a, b):
  # a = (phi, eta); b = (phi, eta)
  delta = np.abs(a-b)
  delta = np.array([2*np.pi - delta[0] if delta[0] > np.pi else delta[0], delta[1]])  # deal with periodicity in phi
  return np.sqrt((delta**2.).sum(axis=-1))


# define helper functions - also a source of parallelization!
def compute_jetDistance(jet1, jet2):
  return distance(np.array(jet1.coord), np.array(jet2.coord))


def match_jets(oJets=[], tJets=[]):
  null_gTower = gTowers.Tower(et=0.0, etamin=0.0, etamax=0.0, phimin=0.0, phimax=0.0, num_cells=0)
  null_gJet = gTowers.Jet(TLorentzVector(), radius=1.0, towers=[null_gTower, null_gTower, null_gTower], seed=null_gTower, area=0.0)
  # dR as a change in distance
  dR = 1.0
  if len(tJets) == 0:
    return np.array([[oJet, null_gJet] for oJet in oJets])
  # we want to match the closest gTower jet for every offline jet
  matched_jets = []
  for oJet in oJets:
    distances = np.array(map(lambda tJet: compute_jetDistance(tJet, oJet), tJets))
    energies = np.array(map(lambda tJet: tJet.Et, tJets))
    # return jet with highest ET within dR
    if np.where(distances <= dR)[0].size == 0:
      closest_tJet = null_gJet
    else:
      max_energy_in_distance = np.amax(energies[np.where(distances <= 1.0)])
      index_jet = np.where(energies == max_energy_in_distance)[0][0]
      closest_tJet = tJets[index_jet]
    matched_jets.append([oJet, closest_tJet])
  return matched_jets


def isolation_condition(leading_jet, jets):
  for jet in jets:
    # skip leading jet
    if jet == leading_jet:
      continue
    if compute_jetDistance(jet, leading_jet) < 2.0:
      return False
  return True

# open up the file and grab the tree
f = TFile.Open(args.filename)
t = f.Get(args.treename)

# set seed cuts
seed_filter = gTowers.SeedFilter(et=args.seedEt_thresh, n=100)
# for offline rho and vxp_n
rho_column_name     = 'Eventshape_rhoKt4LC'
vertices_column_name = 'vxp_nTracks'
weight_column_name = 'weight'  # the weight stored in "weight" is just the "sample weight", or "cross-section weight"
event_weights_column_name = 'mcevt_weight'  # event weight

# offline and trigger jet column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
offline_column_names = ['jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_%s' % col for col in ['pt', 'm', 'eta', 'phi', 'TrimmedSubjetsPtFrac5SmallR30_nsj', 'Tau1', 'Tau2', 'Tau3', 'SPLIT12', 'SPLIT23', 'SPLIT34', 'TrimmedSubjetsPtFrac5SmallR30_index']] + ['jet_AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR30_pt']
gTower_column_names = ['gTower%s' % col for col in ['Et', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]
trig_L1_et_column_names = ['trig_L1_jet_et8x8', 'trig_L1_jet_eta']

# paired jets initialization
paired_jets = []

startTime_wall      = time.time()
startTime_processor = time.clock()

datatype = [('offline_rho', 'float32'),
            ('gFEX_rho_all', 'float32'),
            ('gFEX_rho_1', 'float32'),
            ('gFEX_rho_2', 'float32'),
            ('gFEX_rho_3', 'float32'),
            ('gFEX_rho_4', 'float32'),
            ('vxp_n', 'int32'),
            ('gTower_distribution', 'object'),
            ('weight', 'float32'),
            ('leading_trig_L1_jet_et8x8', 'float32'),
            ('sum_trig_L1_jet_et8x8', 'float32')]

dataLoadTimes = []
rowEvalTimes = []
rhoCalcTimes = []
oEventLoadTimes = []
tEventLoadTimes = []

allBranches = [rho_column_name, vertices_column_name, weight_column_name, event_weights_column_name] + offline_column_names + gTower_column_names + trig_L1_et_column_names

for event_num in xrange(args.event_start, args.event_start+args.num_events, args.step_size):
  t1 = time.time()
  data = rnp.tree2rec(t,
                      branches=allBranches,
                      start=(event_num),
                      stop=(event_num+args.step_size)
                      )
  dataLoadTimes.append(time.time() - t1)

  for row in data:
    t2 = time.time()
    # revive the record array
    event = np.array(row)
    t4 = time.time()
    # jetPt, jetM, jetEta, jetPhi, nsj, tau1, tau2, tau3, split12, split23, split34, subjetsIndex, subjetsPt = event
    oEvent = OfflineJets.Event(event=[event[col].tolist() for col in offline_column_names])

    # if there are no offline jets, we skip it
    if len(oEvent.jets) == 0:
      continue

    # satisfy isolation condition on the offline jets
    if not isolation_condition(oEvent.jets[0], oEvent.jets):
      continue

    oEventLoadTimes.append(time.time() - t4)

    t5 = time.time()
    # at this point, offline event filtering is done, so we go to trigger data
    gTowerData = []
    for col in gTower_column_names:
      # deal with digitization first, before anything else
      if col == 'gTowerEt':
        if args.digitization == 0:
          temp_gTowerEt = event[col].tolist()
        else:
          temp_gTowerEt = np.floor(event[col].tolist()/args.digitization)*args.digitization
        gTowerData.append(temp_gTowerEt)
        hist_gTowerMult, bins_gTowerMult = np.histogram(temp_gTowerEt/1000., bins=np.arange(0.0, 1024.0, 32.0))
        del temp_gTowerEt
      else:
        # no digitization
        gTowerData.append(event[col].tolist())

    # now do offline rho and num vertices
    offline_rho = event[rho_column_name].item()/1000.  # it's an array of one element, so just return the element
    num_vertices = np.sum(event[vertices_column_name].tolist() >= 2)
    event_weight = (event[weight_column_name].item() * event[event_weights_column_name].item()[0][0]) or 1  # it's an array of one element, so just return the element (if it's 0 or 0.0, return 1.0)

    # Michael Begel messaged me this algorithm
    # max([T['trig_L1_jet_et8x8'][k]/1000.0 for k in xrange(T['trig_L1_jet_n']) if abs(T['trig_L1_jet_eta'][k])<3.2]+[1.0])>100.0
    # sum([T['trig_L1_jet_et8x8'][k]/1000.0 for k in xrange(T['trig_L1_jet_n']) if T['trig_L1_jet_et8x8'][k]>20000. if abs(T['trig_L1_jet_eta'][k])<3.2])>200.00
    trig_L1_jet_et = np.append(event['trig_L1_jet_et8x8'].tolist()/1000., 0.0)
    trig_L1_jet_eta = np.append(event['trig_L1_jet_eta'].tolist(), 0.0)
    eta_mask = np.abs(trig_L1_jet_eta) < 3.2
    et_mask = trig_L1_jet_et > 20.

    del event  # extracted all necessary data from it

    tEvent = gTowers.TowerEvent(event=gTowerData, seed_filter=seed_filter)
    tEventLoadTimes.append(time.time() - t5)
    t3 = time.time()

    del gTowerData  # converted that information into tEvent

    gFEX_rho_holder = {1: [], 2: [], 3: [], 4: []}
    for tower in tEvent.towers_below(args.tower_thresh):
      gFEX_rho_holder[tower.region].append(tower.rho)
    gFEX_rho = {'all': np.mean(list(chain(*gFEX_rho_holder.values()))), 1: -1., 2: -1., 3: -1., 4: -1.}
    for region in [1, 2, 3, 4]:
      gFEX_rho[region] = np.mean(gFEX_rho_holder[region])
    rhoCalcTimes.append(time.time() - t3)

    towers = tEvent.towers_above(args.noise_filter)
    tEvent.get_event(towers=towers)

    # only do the leading jet for now
    matched_jets = match_jets(oJets=[oEvent.jets[0]], tJets=tEvent.event.jets)

    leading_oJet = matched_jets[0][0]
    matched_tJet = matched_jets[0][1]

    event_data = np.array([leading_oJet.as_rec.tolist()[0]
                           + matched_tJet.as_rec.tolist()[0]
                           + (offline_rho,
                              gFEX_rho['all'],
                              gFEX_rho[1],
                              gFEX_rho[2],
                              gFEX_rho[3],
                              gFEX_rho[4],
                              num_vertices,
                              hist_gTowerMult,
                              event_weight,
                              np.amax(trig_L1_jet_et[np.where(eta_mask)]),
                              np.sum(trig_L1_jet_et[np.where(et_mask & eta_mask)]))
                           ], dtype=(leading_oJet.as_rec.dtype.descr + matched_tJet.as_rec.dtype.descr + datatype))

    paired_jets.append(event_data)

    rowEvalTimes.append(time.time() - t2)

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job %d in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % (args.process_num, (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

'''at this point, we've processed all the data and we just need to dump it'''
filename_ending = 'unweighted_seed%d_noise%d_signal%d_digitization%d_process%d' % (args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization, args.process_num)
pickle.dump(np.hstack(paired_jets), file(write_file('data/seed%d/matched_jets_%s.pkl' % (args.seedEt_thresh, filename_ending)), 'w+'))
print len(paired_jets)
print "dataLoadTimes", np.mean(dataLoadTimes), np.std(dataLoadTimes), len(dataLoadTimes)
print "rowEvalTimes", np.mean(rowEvalTimes), np.std(rowEvalTimes), len(rowEvalTimes)
print "\toEventLoadTimes", np.mean(oEventLoadTimes), np.std(oEventLoadTimes), len(oEventLoadTimes)
print "\ttEventLoadTimes", np.mean(tEventLoadTimes), np.std(tEventLoadTimes), len(tEventLoadTimes)
print "\trhoCalcTimes", np.mean(rhoCalcTimes), np.std(rhoCalcTimes), len(rhoCalcTimes)
