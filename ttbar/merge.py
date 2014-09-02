import numpy as np
from glob import glob
import os
import argparse
import time

#randomize list of files
from random import shuffle

try:
  import cPickle as pickle
except:
  import pickle


parser = argparse.ArgumentParser(description='Process pickled data to make turn-on curves')
parser.add_argument('--seedEt',       type=float, required=True,  dest='seedEt_thresh',   help='seed Et Threshold (GeV)')
parser.add_argument('--noiseFilter',  type=float, required=True,  dest='noise_filter',    help='gTower noise filter (GeV)')
parser.add_argument('--towerThresh',  type=float, required=True,  dest='tower_thresh',    help='gTower rho Et thresh (GeV)')
parser.add_argument('--digitization', type=float, required=True,  dest='digitization',    help='digitization for gTower Et (MeV)')
parser.add_argument('--numFiles',     type=int,   required=False, dest='numFiles',        help='number of files to read in for merging', default=-1)

# parse the arguments, throw errors if missing any
args = parser.parse_args()

files = glob("data/seed%d/matched_jets_unweighted_seed%d_noise%d_signal%d_digitization%d_process*.pkl" % (args.seedEt_thresh, args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization) )
shuffle(files)


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def write_file(f):
  ensure_dir(f)
  return f

startTime_wall      = time.time()
startTime_processor = time.clock()

raw_data = np.array([])

numFiles = min(len(files), args.numFiles)

if numFiles > 0:
  filesToLoop = files[:numFiles]
else:
  filesToLoop = files

for filename in filesToLoop:
  data = pickle.load(file(filename))
  if raw_data.size > 0:
    raw_data = np.append(raw_data, data)
  else:
    raw_data = data
  print len(data), len(raw_data), filename

print "\nStats\n-------\n\t%d files\n\t%d events\n\t%d events per file" % (len(filesToLoop), len(raw_data), (len(raw_data)*1./len(filesToLoop)) )

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: %0.2f s \n\t Clock Time: %0.2f s" % ( (endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

pickle.dump(raw_data, file( write_file('data/seed%d/leading_jets_seed%d_noise%d_signal%d_digitization%d.pkl' % (args.seedEt_thresh, args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization)), 'w+') )
