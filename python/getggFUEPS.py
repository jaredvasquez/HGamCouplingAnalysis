from ROOT import *
from HGamMoriondCatsBDT import CatLabels

fnom = TFile('output/HGamCoupling_ggH_NNLOPS/hist-ggH_NNLOPS.root')
fsys = TFile('output/HGamCoupling_ggH_H7/hist-ggH_H7.root')

histName = 'h_catSTXS'
hnom = fnom.Get( histName )
hsys = fsys.Get( histName )

hdiff = hnom.Clone()
print '\n  UE/PS Uncertainty '
print '-'*45
for ibin in xrange(1,hdiff.GetNbinsX()+1):
  nom, sys = hnom.GetBinContent(ibin), hsys.GetBinContent(ibin)
  if not nom or not sys:
    hdiff.SetBinContent( ibin, 0. )
    continue
  print '%25s : %6.2f %%' % ( CatLabels[ibin-1], (sys/nom-1)*100. )
  hdiff.SetBinContent( ibin, (1-sys/nom)*100. )
print ''
