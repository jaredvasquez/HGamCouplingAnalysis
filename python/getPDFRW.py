import json
import ROOT
import HGamMoriondCatsBDT as HG

NPnames  = ['nominal']
NPnames += [ 'ATLAS_PDF4LHC_NLO_30_EV%d' % (i+1) for i in xrange(30) ]
NPnames += [ 'ATLAS_PDF4LHC_NLO_30_alphaS'+suf for suf in '__1up __1dn'.split() ]

nSys = 33
procs = ['ttH','bbH','tHW','tHjb']

# -----------------------------------------------------------------------
def fixPrecision( N ):
  return float( '%.6f' % N )

# -----------------------------------------------------------------------
def getDiff(nom, sys):
  if (nom==0): return 0.0
  return abs(sys/nom - 1.0)


# -----------------------------------------------------------------------
def pruneSysts( allSys ):
  for proc in allSys:
    for cat in allSys[proc]:
      sysNames = allSys[proc][cat].keys()
      for sysName in sysNames:
        pruneSys = False
        hi, lo, form = allSys[proc][cat][sysName]
        allSys[proc][cat][sysName] = (abs(hi), -1.0*abs(lo), form) ### Assume positive correlations
        if 'EG_' in sysName:
          pruneSys = True
        if (hi == 0. or lo == 0.):
          pruneSys = True
        elif (abs(lo/hi) > 5 or abs(hi/lo) > 5):
          pruneSys = True
        elif lo*hi > 0:
           # assume positive correlations when up/down have same trends
           #  --> could also symmetrize with max?
           #  --> many stat dominated systs, need a way to check this... (N RAW)
          allSys[proc][cat][sysName] = (abs(hi), -1.0*abs(lo), form)
        #elif (sysTOT - sysNM1 < 0.0001):
        #  pruneSys = True
        if pruneSys:
          del  allSys[proc][cat][sysName]
  return allSys


# -----------------------------------------------------------------------
def getSys():
  sysMap = { p: { cat : {} for cat in HG.CatLabels } for p in procs }
  for proc in procs:
    tf = ROOT.TFile('output/HGamPDF_{0}/hist-{0}.root'.format(proc))
    hnom = tf.Get('h_category_PDFVar0')
    for ipdf in xrange(1, nSys):
      hsys = tf.Get('h_category_PDFVar{}'.format(ipdf))
      for icat, cat in enumerate(HG.CatLabels):
        nom = hnom.GetBinContent(icat+1)
        sys = hsys.GetBinContent(icat+1)
        err = fixPrecision( getDiff(nom, sys) )
        sysMap[proc][cat][NPnames[ipdf]] = [ err, err, 'logn' ]

    # collet asymm uncerts
    for NPup in NPnames:
      if '__1up' not in NPup: continue
      NP = NPup.replace('__1up','')
      NPdn = NP+'__1dn'
      for cat in HG.CatLabels:
        errHI = sysMap[proc][cat][NPup][0]
        errLO = sysMap[proc][cat][NPdn][0]
        sysMap[proc][cat][NP] = [ errHI, errLO, 'asym' ]
        del sysMap[proc][cat][NPup]
        del sysMap[proc][cat][NPdn]

  return sysMap


# -----------------------------------------------------------------------
if __name__ == "__main__":
  sysMap = getSys()
  #sysMap = pruneSysts(sysAll)
  json.dump(sysMap, open('pdfSysRW.json','wb'))
  for proc in procs:
    cat = 'ttHlep'
    print proc
    for NP in sorted(sysMap[proc][cat]):
      print '  -->', NP, ':', sysMap[proc][cat][NP]

