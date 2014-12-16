import numpy as np
from glob import glob
import os
import argparse
import time
import re

# randomize list of files
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

files = []
for JZXW in range(4):
  files = files + glob("../dijets_{:d}W/data/seed{:0.0f}/matched_jets_unweighted_seed{:0.0f}_noise{:0.0f}_signal{:0.0f}_digitization{:0.0f}_process*.pkl".format(JZXW, args.seedEt_thresh, args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization))
shuffle(files)


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def write_file(f):
  ensure_dir(f)
  return f


def humansize(nbytes):
  suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
  if nbytes == 0:
    return '0 B'
  i = 0
  while nbytes >= 1024 and i < len(suffixes)-1:
    nbytes /= 1024.
    i += 1
  f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
  return '%s %s' % (f, suffixes[i])


startTime_wall      = time.time()
startTime_processor = time.clock()

numFiles = min(len(files), args.numFiles)

p = re.compile(r'(?<=dijets_\dW/).*(?=process)')

if numFiles > 0:
  filesToLoop = files[:numFiles]
else:
  filesToLoop = files

dataHolder = []
dataSize = 0
counter = 0
for filename in filesToLoop:
  dataHolder.append(pickle.load(file(filename)))
  rawDataSize = dataHolder[-1].size
  dataSize += rawDataSize
  counter += 1
  print "{:3d}/{:3d} | {:3d} | {:7d} | {}".format(counter, len(filesToLoop), rawDataSize, dataSize, p.sub('', filename))

print "\nStats\n-------\n\t%d files\n\t%d events\n\t%d events per file" % (len(filesToLoop), dataSize, (dataSize*1./len(filesToLoop)))

endTime_wall      = time.time()
endTime_processor = time.clock()
print "Finished job in:\n\t Wall time: {:0.2f}s \n\t Clock Time: {:0.2f}s".format((endTime_wall - startTime_wall), (endTime_processor - startTime_processor))

fileToWrite = write_file('data/seed{:0.0f}/leading_jets_seed{:0.0f}_noise{:0.0f}_signal{:0.0f}_digitization{:0.0f}.pkl'.format(args.seedEt_thresh, args.seedEt_thresh, args.noise_filter, args.tower_thresh, args.digitization))
pickle.dump(np.hstack(dataHolder), file(fileToWrite, 'w+'))
print "File Written: {}\nFile Size: {}".format(fileToWrite, humansize(os.path.getsize(fileToWrite)))
