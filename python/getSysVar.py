import os, sys, math
from sysTools import *
from ROOT import *
import HGamMoriondCats as HG

# Parameters
_NSYS = 120
procs = ['ggH','VBF','WH','ZH','ttH','bbH','tHjb','tWH']
strPATH = 'output/Systematics/HGamSys_%s/hist-%s.root'

#_NSYS = 5
#procs = ['ggH','WH','ZH']
#strPATH = 'output/PHSFtest_%s/hist-%s.root'


# List all Categories
listCats = False
if (listCats):
  print "Categories (%d Total):" % (len(HG.CatLabels))
  for label in HG.CatLabels:
    print "  ", label
  print ""


# Get Systematics Variations for each process
missingFiles = False
allSys = {}
for proc in procs:
  tfs = []
  for i in xrange(0,_NSYS+1):
    sampName = "%s_%d" % (proc,i)
    filePATH = strPATH % (sampName, proc)
    #filePATH = strPATH % (sampName, sampName)
    if not os.path.isfile(filePATH):
      print "ERROR: Could not find", sampName
      missingFiles = True
      continue
    tfs.append( TFile(filePATH) )
  if missingFiles: sys.exit()
  allSys[proc] = getVariation(tfs)


# Organize All Systematics By: ( Category | Proc | Systematic )
allSysByCat = { cat: { proc: allSys[proc][cat] for proc in procs } for cat in HG.CatLabels }


# Print results
skipSys = ['Trig', 'EG_']

# Print results
for cat in HG.CatLabels:
  print '\n\nCategory: %s' % cat
  print '-'*40 + '\n'
  for proc in procs: #allSysByCat[cat]:
    print " ", proc
    mcstat = abs( allSysByCat[cat][proc]['MCstat'][0] )*100.
    minsys = math.sqrt( (mcstat+1.0)**2 - mcstat**2 )
    syslist = sorted(list(allSysByCat[cat][proc]))
    #for sys in sorted(list(allSysByCat[cat][proc])):
    for sys in syslist + [ syslist.pop(syslist.index('MCstat')) ]:
      if 'EG_' in sys: continue
      if 'Trig' in sys: continue
      #if not 'PH_' in sys: continue

      sigup, sigdn = allSysByCat[cat][proc][sys]
      sigup *= 100.; sigdn *= 100.

      tagasym = ''
      sigabs = (abs(sigup), abs(sigdn))

      # if one or more systs are zero, assume low stats and prune
      if not min(sigabs): continue

      # if +/- variations are more than 100% diff
      if ( max(sigabs)/min(sigabs) > 2 ): tagasym = '***'

      if ( sys != 'MCstat' and max(sigabs) < minsys ): continue
      #if not float( '%.2f' % sum( map( abs, list(allSysByCat[cat][proc][sys]) ) ) ): continue
      print "    -->", sys, "+/- ( %+.2f %%, %+.2f %% )" % (sigup,sigdn), tagasym

print ""
