from atlas_jets import *
import numpy as np
import root_numpy as rnp
from ROOT import TFile

from itertools import chain

import os
try:
  import cPickle as pickle
except:
  import pickle

import argparse
import time

parser = argparse.ArgumentParser(description='Process TTbar events for Level 1 Trigger.')
parser.add_argument('--processNum',   type=int,   required=False, dest='process_num',     help='a process number for tagging files', default=0)
parser.add_argument('--file',         type=str,   required=True,  dest='filename',        help='a filename to read in')
parser.add_argument('--start',        type=int,   required=True,  dest='event_start',     help='start index of event')
parser.add_argument('--numEvents',    type=int,   required=True,  dest='num_events',      help='number of events to process')
parser.add_argument('--stepSize',     type=int,   required=False, dest='step_size',       help='chunking of events', default=10)
parser.add_argument('--tree',         type=str,   required=False, dest='treename',        help='treename', default='mytree')
# parse the arguments, throw errors if missing any
args = parser.parse_args()

print args

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def write_file(f):
  ensure_dir(f)
  return f

# open up the file and grab the tree
f = TFile.Open(args.filename)
t = f.Get(args.treename)

#offline and trigger jet column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
offline_column_names = ['jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_%s' % col for col in ['pt', 'm', 'eta', 'phi', 'TrimmedSubjetsPtFrac5SmallR30_nsj', 'Tau1', 'Tau2', 'Tau3', 'SPLIT12', 'SPLIT23', 'SPLIT34', 'TrimmedSubjetsPtFrac5SmallR30_index']] + ['jet_AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR30_pt']

#paired jets initialization
out_data = []

startTime_wall      = time.time()
startTime_processor = time.clock()

datatype = ([('leading_offline_jet','object'),\
              ('leading_trig_L1_jet_et8x8','float32'),\
              ('sum_trig_L1_jet_et8x8', 'float32')])

for event_num in xrange(args.event_start, args.event_start+args.num_events, args.step_size):
  data = rnp.tree2rec(t, branches=offline_column_names + ['trig_L1_jet_et8x8','trig_L1_jet_eta'], start=(event_num), stop=(event_num+args.step_size))
  for row in data:
    #revive the record array
    event = np.array(row)
    # jetE, jetPt, jetM, jetEta, jetPhi, nsj, tau1, tau2, tau3, split12, split23, split34, subjetsIndex, subjetsPt = event
    oEvent = OfflineJets.Event(event=[event[col].tolist() for col in offline_column_names])

    # if there are no offline jets, we skip it
    if len(oEvent.jets) == 0:
      continue

    # max([T['trig_L1_jet_et8x8'][k]/1000.0 for k in xrange(T['trig_L1_jet_n']) if abs(T['trig_L1_jet_eta'][k])<3.2]+[1.0])>100.0
    # sum([T['trig_L1_jet_et8x8'][k]/1000.0 for k in xrange(T['trig_L1_jet_n']) if T['trig_L1_jet_et8x8'][k]>20000. if abs(T['trig_L1_jet_eta'][k])<3.2])>200.00

    trig_L1_jet_et = np.append(event['trig_L1_jet_et8x8'].tolist()/1000., 0.0)
    trig_L1_jet_eta = np.append(event['trig_L1_jet_eta'].tolist(), 0.0)
    eta_mask = np.abs(trig_L1_jet_eta) < 3.2
    et_mask = trig_L1_jet_et > 20.

    event_data = np.array([(oEvent.jets[0], np.amax(trig_L1_jet_et[np.where(eta_mask)]), np.sum(trig_L1_jet_et[np.where(et_mask&eta_mask)]) )], dtype=datatype)

    out_data.append(event_data)

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job %d in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % (args.process_num, (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

'''at this point, we've processed all the data and we just need to dump it'''
filename_ending = 'trig_L1_et_process%d' % (args.process_num)
pickle.dump(np.hstack(out_data), file( write_file('data/%s.pkl' % (filename_ending) ), 'w+') )
print len(out_data)
