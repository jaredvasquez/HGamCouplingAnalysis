from ROOT import *
from math import sqrt

binNames = ["UNKNOWN","GG2H_FWDH","GG2H_VBFTOPO_JET3VETO","GG2H_VBFTOPO_JET3","GG2H_0J","GG2H_1J_PTH_0_60","GG2H_1J_PTH_60_120","GG2H_1J_PTH_120_200","GG2H_1J_PTH_GT200","GG2H_GE2J_PTH_0_60","GG2H_GE2J_PTH_60_120","GG2H_GE2J_PTH_120_200","GG2H_GE2J_PTH_GT200","VBF_QQ2HQQ_FWDH","QQ2HQQ_VBFTOPO_JET3VETO","QQ2HQQ_VBFTOPO_JET3","QQ2HQQ_VH2JET","QQ2HQQ_REST","QQ2HQQ_PTJET1_GT200","QQ2HLNU_FWDH","QQ2HLNU_PTV_0_150","QQ2HLNU_PTV_150_250_0J","QQ2HLNU_PTV_150_250_GE1J","QQ2HLNU_PTV_GT250","QQ2HLL_FWDH","QQ2HLL_PTV_0_150","QQ2HLL_PTV_150_250_0J","QQ2HLL_PTV_150_250_GE1J","QQ2HLL_PTV_GT250","GG2HLL_FWDH","GG2HLL_PTV_0_150","GG2HLL_PTV_GT150_0J","GG2HLL_PTV_GT150_GE1J","TTH_FWDH","TTH","BBH_FWDH","BBH","THJB_FWDH","THJB","TWH_FWDH","TWH"]

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
procs=['tWH','tHjb']
histName = 'h_truthAcc_weightMC'

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
  for ibin, binName in enumerate( binNames ):
    hbin.GetXaxis().SetBinLabel( ibin+1, binName )
    #print ibin, binName

  hsum = getSumHist( hbin )
  hbin.Divide( hbin, hsum, 1.0, 1.0, 'b' ) #binomial errors

  print '\n  %s' % procName[proc]
  print '-'*60
  for ibin, binName in enumerate( binNames ):
    if (ibin < 1): continue
    acc = hbin.GetBinContent( ibin+1 )
    err = hbin.GetBinError( ibin+1 )
    if (acc != 0.0):
      print '%35s : %6.4f +/- %6.4f' % (binName, acc, err)
  print ''
