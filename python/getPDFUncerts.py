import ROOT
import HGamMoriondCatsBDT as HG

def getDiff(nom, sys):
  if (nom==0.): return 0.
  return (sys/nom-1.0)


def getPDFUncerts(tf):
  hnom  = tf.Get('h_catSTXS')
  othersys = '_alphaS_up _alphaS_dn'.split()
  hsys  = [ tf.Get('h_catSTXS_PDF%d' % i) for i in xrange(30) ]
  hsys += [ tf.Get('h_catSTXS'+suf) for suf in othersys ]
  systMap = { cat : {} for cat in HG.CatLabels }

  for ipdf in xrange(30):
    for icat, catName in enumerate(HG.CatLabels):
      nom = hnom.GetBinContent(icat+1)
      sys = hsys[ipdf].GetBinContent(icat+1)
      systMap[catName]['ATLAS_PDF4LHC_NLO_30_EV%d' % (ipdf+1)] = getDiff(nom, sys)

  for catName in reversed(HG.CatLabels):
    print '\n  '+catName
    for ipdf in xrange(30):
      sysName = 'ATLAS_PDF4LHC_NLO_30_EV%d' % (ipdf+1)
      print '    --> %s : %.2f%%' % (sysName, abs(systMap[catName][sysName]*100.))
      #print '    -->', sysName, ':', abs(systMap[catName][sysName]

if __name__ == "__main__":
  # execute only if run as a script
  procs = ['VBF_PDF']
  procs = ['WH_NLO']
  #procs = ['ggH_NNLOPS']

  for proc in procs:
    banner = '===      %s      ===' % proc
    print '='*len(banner)
    print banner
    print '='*len(banner)
    tf = ROOT.TFile('output/HGamCoupling_{proc}/hist-{proc}.root'.format(proc=proc))
    getPDFUncerts(tf)
