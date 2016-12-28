from ROOT import *
from math import sqrt
from STXSFineIndex import FineIndexLabels

procName = {
   'ggH' : 'ggH (Powheg + Pythia8), with pTH rweighting',
   'VBF' : 'VBF (Powheg + Pythia8)',
   'WH'  : 'WH (Pythia8)',
   'ZH'  : 'ZH (Pythia8)',
   'ttH' : 'ttH (aMC@NLO + Pythia8)',
   'bbH' : 'bbH (aMC@NLO + Pythia8)',
   'tWH' : 'tWH (aMC@NLO + Herwig++)',
  'tHjb' : 'tHjb (MadGraph + Pythia8)',
}

procs=['ggH','VBF','WH','ZH','ttH','bbH','tWH','tHjb']
#procs=['ggH','VBF','WH','ZH','ttH','bbH']
histName = 'h_truthAcc_fineIndex_weightMC'

def getSumHist( h ):
  Ntot, Nerr = 0., 0.
  hsum = h.Clone('hsum')

  # Skip underflow and unkown, right?
  for ibin in xrange( 2, h.GetNbinsX()+1 ):
    Ntot += h.GetBinContent(ibin)
    Nerr += h.GetBinError(ibin)**2

  Nerr = sqrt(Nerr)
  for ibin in xrange( 2, h.GetNbinsX()+1 ):
    hsum.SetBinContent(ibin, Ntot)
    hsum.SetBinError(ibin, Nerr)

  return hsum


for proc in procs:
  tf = TFile( 'output/Coupling_%s/hist-%s.root' % (proc,proc) )
  hbin = tf.Get(histName)
  for ibin, binName in enumerate( FineIndexLabels ):
    hbin.GetXaxis().SetBinLabel( ibin+1, binName )
    #print ibin, binName

  hsum = getSumHist( hbin )
  hbin.Divide( hbin, hsum, 1.0, 1.0, 'b' ) #binomial errors

  print '\n//  %s' % procName[proc]
  print '-'*60
  for ibin, binName in enumerate( FineIndexLabels ):
    if (ibin < 1): continue
    acc = hbin.GetBinContent( ibin+1 )
    err = hbin.GetBinError( ibin+1 )
    if (acc != 0.0):
      print '%35s : %6.4f +/- %6.4f' % (binName, acc, err)
  print ''
