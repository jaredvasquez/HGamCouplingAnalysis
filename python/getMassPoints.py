from ROOT import *
from HGamMoriondCats import CatLabels

nCats = len(CatLabels)
datasets=['data15_all','data16_DS1','data16_DS2','data16_DS4','data16_DS5','data16_DS6']
massPoints = { icat : [] for icat in xrange(1,nCats+1) }

for dataset in datasets:
  tf = TFile('output/Coupling_%s/hist-%s.root' % (dataset,dataset))
  nt = tf.Get('physics')
  for ievt in xrange( nt.GetEntries() ):
    nt.GetEntry( ievt )
    massPoints[nt.category].append( nt.m_yy )
    #print '%3d  %6.3f' % (nt.category, nt.m_yy)

# Output each categories' data points to text file.

print ''
print '   Data Summary:'
print '-'*40
for icat in xrange(nCats):
  print '  %15s  :  %8d events' % (CatLabels[icat], len(massPoints[icat+1]))
print ''
