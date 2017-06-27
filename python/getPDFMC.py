import os
import sys
import json
import ROOT
import tabulate
from math import sqrt

tabulate.LATEX_ESCAPE_RULES={}
sys.dont_write_bytecode = True
import HGamMoriondCatsBDT as HG

NPnames  = [ 'ATLAS_PDF4LHC_NLO_30_EV%d' % (i+1) for i in xrange(30) ]
NPnames += [ 'ATLAS_PDF4LHC_NLO_30_alphaS'+suf for suf in '_up _dn'.split() ]

nSys = 30
procs = ['ggH_NNLOPS','VBF_NNPDF','WH_NLO','ZH_NLO'] #ttH','bbH','tHW','tHjb']

# -----------------------------------------------------------------------
def fixPrecision( N ):
  return float( '%.4f' % N )


# -----------------------------------------------------------------------
def getDiff(nom, sys):
  if (nom==0): return 0.0
  return (sys/nom - 1.0)


# -----------------------------------------------------------------------
def formatCell(d, cat, sys):
  if cat not in d:
    return '---'
  if sys not in d[cat]:
    return '---'
  else:
    hi, lo, form = d[cat][sys]
    if form == 'logn':
      if float('%.1f'%(hi*100)) == 0.:
        return '---'
      return '%+.1f' % ((hi*100))
    return '${}^{%+.1f}_{%+.1f}$' % (hi*100, lo*100)


# -----------------------------------------------------------------------
def pruneSysts( allSys ):
  for proc in allSys:
    for cat in allSys[proc]:
      sysNames = allSys[proc][cat].keys()
      for sysName in sysNames:
        pruneSys = False
        hi, lo, form = allSys[proc][cat][sysName]
        hi = float('%.3f' % hi)
        lo = float('%.3f' % lo)
        allSys[proc][cat][sysName] = ( hi, lo, form )
        #allSys[proc][cat][sysName] = (abs(hi), -1.0*abs(lo), form) ### Assume positive correlations
        if (hi == 0. or lo == 0.):
          pruneSys = True
        elif (abs(lo/hi) > 5 or abs(hi/lo) > 5):
          pruneSys = True
        # ****** IMPORTANT *********
        #elif (abs(hi) < 0.0025 or abs(lo) < 0.0025):
        #elif (abs(hi) < 0.005 or abs(lo) < 0.005):
        #  pruneSys = True
        elif lo*hi > 0:
          # assume positive correlations when up/down have same trends
          #  --> could also symmetrize with max?
          #  --> many stat dominated systs, need a way to check this... (N RAW)
          allSys[proc][cat][sysName] = ( hi, -1.0*lo, form)
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
    hnom = tf.Get('h_catSTXS')
    for ipdf in xrange(nSys):
      hsys = tf.Get('h_catSTXS_PDF{}'.format(ipdf))
      for icat, cat in enumerate(HG.CatLabels):
        nom = hnom.GetBinContent(icat+1) #/ hnom.GetBinContent(0)
        sys = hsys.GetBinContent(icat+1) #/ hsys.GetBinContent(0)
        err = fixPrecision( getDiff(nom, sys) )
        sysMap[proc][cat][NPnames[ipdf]] = [ err, err, 'logn' ]

    for suffix in [ '_alphaS_up', '_alphaS_dn']:
      hsys = tf.Get('h_catSTXS'+suffix)
      for icat, cat in enumerate(HG.CatLabels):
        nom = hnom.GetBinContent(icat+1) #/ hnom.GetBinContent(0)
        sys = hsys.GetBinContent(icat+1) #/ hsys.GetBinContent(0)
        err = fixPrecision( getDiff(nom, sys) )
        sysMap[proc][cat]['ATLAS_PDF4LHC_NLO_30'+suffix] = [ err, err, 'logn' ]

    # collet asymm uncerts
    for NPup in NPnames:
      if '_up' not in NPup: continue
      NP = NPup.replace('_up','')
      NPdn = NP+'_dn'
      for cat in HG.CatLabels:
        errHI = sysMap[proc][cat][NPup][0]
        errLO = sysMap[proc][cat][NPdn][0]
        sysMap[proc][cat][NP] = [ errHI, errLO, 'asym' ]
        del sysMap[proc][cat][NPup]
        del sysMap[proc][cat][NPdn]

  return sysMap


# -----------------------------------------------------------------------
if __name__ == "__main__":
  sysAll = getSys()
  #sysAll.update( json.load(open('pdfSysRW.json')) )
  sysAll.update( json.load(open('pdfSysRW_yield.json')) )

  #sysMap = sysAll
  sysMap = pruneSysts(sysAll)
  json.dump(sysMap, open('pdfSys.json','wb'))

  uniqNPs = set()
  for p in sysMap:
    for c in sysMap[p]:
      for NP in sysMap[p][c]:
        uniqNPs.add(NP)

  print ' Unique NPs ( %d )' % len(uniqNPs)
  print '-'*35
  for NP in sorted(uniqNPs):
    print NP

  print ''
  for p in sysMap: print p

  #  Prepare table
  # ----------------------------------------------------------------------
  os.system('mkdir -p tablesPDF')
  NPnames  = [ 'ATLAS_PDF4LHC_NLO_30_EV%d' % (i+1) for i in xrange(30) ]
  NPnames += [ 'ATLAS_PDF4LHC_NLO_30_alphaS' ]

  for proc in sysMap:
    table = []
    headers = ['Systematic'] + [ i+1 for i in xrange(0,17) ]
    for sys in NPnames:
      line = [sys.replace('_','\_')] + [ formatCell(sysMap[proc], cat, sys) for cat in HG.CatLabels[:17] ]
      table.append(line)
    log = open('tablesPDF/table_%s1.tex'%proc,'w+')
    print >> log, tabulate.tabulate( table, headers=headers, tablefmt='latex' )
    #if 'ggH' in proc:
    #  print tabulate.tabulate( table, headers=headers )

    table = []
    headers = ['Systematic'] + [ i+1 for i in xrange(17,len(HG.CatLabels)+1) ]
    for sys in NPnames:
      line = [sys.replace('_','\_')] + [ formatCell(sysMap[proc], cat, sys) for cat in HG.CatLabels[17:] ]
      table.append(line)
    log = open('tablesPDF/table_%s2.tex'%proc,'w+')
    print >> log, tabulate.tabulate( table, headers=headers, tablefmt='latex' )

  import sys
  sys.exit()
  print '\n\n'
  for cat in HG.CatLabels:
    print ''
    print cat
    print '-'*35
    for proc in sysMap:
      #cat = 'GGH_0J_CEN'
      print proc
      tot = 0
      tot2 = 0
      for NP in sorted(sysMap[proc][cat]):
        syst = max(sysMap[proc][cat][NP][0:1])
        tot += syst**2
        if abs(syst) >= 0.0025:
          tot2 += syst**2
          #print '  -->', NP, ':', sysMap[proc][cat][NP]
      print '  --> TOTAL  : ', sqrt(tot)
      print '  --> PRUNED : ', sqrt(tot2)


  #for icat, cat in enumerate(HG.CatLabels):
  #  print icat+1, cat
