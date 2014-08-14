from atlas_jets import *
import numpy as np
import root_numpy as rnp
from ROOT import TFile

for i in range(13,41):
  filename = 'root://faxbox.usatlas.org//user/kratsg/LArStudies/PileupSkim_PowhegPythia_P2011C_ttbar_%d.root' % i
  print "--- opening file %d" % i
  f = TFile.Open(filename)
  t = f.Get('mytree')
  offline_column_names = ['jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi', 'TrimmedSubjetsPtFrac5SmallR30_nsj', 'Tau1', 'Tau2', 'Tau3', 'SPLIT12', 'SPLIT23', 'SPLIT34', 'TrimmedSubjetsPtFrac5SmallR30_index']] + ['jet_AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR30_pt']
  for event_num in xrange(0,10000,100):
    data = rnp.tree2rec(t, branches = offline_column_names, start=(event_num), stop=(event_num+100))
    j = 0
    for row in data:
      event=np.array(row)
      oEvent = OfflineJets.Event(event=[event[col].tolist() for col in offline_column_names])
      for jet in oEvent.jets:
        if 274.7 < jet.Pt < 274.9 and jet.nsj == 3:
          print "i = %d, e = %d, j = %d" % (i,event_num, j)
      j += 1
  print "--- %d" % i


#example function of how to pull the correct event given reporting above
def print_jets(i, event_num, j):
  index = event_num + j
  filename = 'root://faxbox.usatlas.org//user/kratsg/LArStudies/PileupSkim_PowhegPythia_P2011C_ttbar_%d.root' % i
  f = TFile.Open(filename)
  t = f.Get('mytree')
  data = rnp.tree2rec(t, branches = offline_column_names, start=(index), stop=(index+1))
  for row in data:
    event=np.array(row)
    oEvent = OfflineJets.Event(event=[event[col].tolist() for col in offline_column_names])
    for jet in oEvent.jets:
      print jet


