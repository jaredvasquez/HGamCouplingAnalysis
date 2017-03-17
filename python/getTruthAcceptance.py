from ROOT import *
from math import sqrt
from STXSFineIndex import FineIndexLabels

procName = {
   #'ggH' : 'ggH (Powheg + Pythia8), with pTH rweighting',
   #'VBF' : 'VBF (Powheg + Pythia8)',
   #'WH'  : 'WH (Pythia8)',
   #'ZH'  : 'ZH (Pythia8)',
   'ggH_NNLOPS' : 'ggH (NNLOPS)',
    'VBF_NNPDF' : 'VBF (Powheg + Pythia8, NNPDF)',
       'WH_NLO' : 'WH (Powheg + Pythia8)',
       'ZH_NLO' : 'ZH (Powheg + Pythia8)',
          'ttH' : 'ttH (aMC@NLO + Pythia8)',
          'bbH' : 'bbH (aMC@NLO + Pythia8)',
          'tHW' : 'tHW (aMC@NLO + Herwig++)',
         'tHjb' : 'tHjb (MadGraph + Pythia8)',
# tH samples
#    'tHW_plus2' : 'tHW, k=+2 (aMC@NLO + Herwig++)',
#   'tHjb_plus2' : 'tHjb, k=+2 (MadGraph + Pythia8)',
#   'tHW_minus1' : 'tHW, k=-1 (aMC@NLO + Herwig++)',
#  'tHjb_minus1' : 'tHjb, k=-1 (MadGraph + Pythia8)',
}


# Cross Sections and Branching Ratios
BR = 2.270E-03
xsec = {          # [pb]
  #       'ggH' : 4.852E+01,
  #       'VBF' : 3.779E+00,
  #       'WH'  : 1.369E+00,
  #       'ZH'  : 8.824E-01,
  'ggH_NNLOPS' : 4.852E+01,
   'VBF_NNPDF' : 3.779E+00,
      'WH_NLO' : 1.369E+00,
      'ZH_NLO' : 8.824E-01,
         'ttH' : 5.065E-01,
         'bbH' : 4.863E-01,
         'tHW' : 2.500E-02,
        'tHjb' : 7.425E-02,
   'tHW_plus2' : 9.700E-02,
  'tHjb_plus2' : 2.685E-01,
  'tHW_minus1' : 1.504E-01,
 'tHjb_minus1' : 7.320E-01,
}
binXS = []


#procs=['ggH','VBF','WH','ZH','ttH','bbH','tHW','tHjb','tHW_plus2','tHjb_plus2','tHW_minus1','tHjb_minus1']
#procs=['ggH','VBF','WH','ZH','ttH','bbH','tHW','tHjb']
procs=['ggH_NNLOPS','VBF_NNPDF','WH_NLO','ZH_NLO','ttH','bbH','tHW','tHjb']
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
  tf = TFile( 'output/HGamPDF_%s/hist-%s.root' % (proc,proc) )
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
      binXS.append( ( binName, acc*BR*xsec[proc]*1000) )
  print ''
print '\n\n'

# Get Cross-Sections
import collections
counts = collections.Counter(binXS)
dups = [i for i in counts if counts[i]>1] # Check for duplicates
if dups: print 'Uh oh! Duplicates of bins:', dups
print '   STXS Cross-Sections [fb]'
print '-'*45
for binName, XS in binXS:
  print '%30s : %8.3E,' % ('\'%s\'' % binName, XS)
