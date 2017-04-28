import json
from ROOT import *
from HGamMoriondCatsBDT import CatLabels

fnom = TFile('output/HGamCoupling_ttH/hist-ttH.root')
fsys = TFile('output/HGamCoupling_ttH_Hwpp/hist-ttH_Hwpp.root')

histName = 'h_catSTXS'
hnom = fnom.Get( histName )
hsys = fsys.Get( histName )


output = { c: {} for c in CatLabels }

hdiff = hnom.Clone()
print '\n  PS Generator Uncertainty '
print '-'*45
for ibin in xrange(1,hdiff.GetNbinsX()+1):
  nom, sys = hnom.GetBinContent(ibin), hsys.GetBinContent(ibin)
  if not nom or not sys:
    hdiff.SetBinContent( ibin, 0. )
    continue
  print '%25s : %6.2f %%' % ( CatLabels[ibin-1], (sys/nom-1)*100. )
  diff = (sys/nom-1)
  output[CatLabels[ibin-1]]['ATLAS_UEPS_ttH'] = [diff, -diff, 'logn']
  hdiff.SetBinContent( ibin, (1-sys/nom)*100. )
print ''

output = { 'ttH' : output }
json.dump(output, open('syst_ttHPS.json','wb'))
