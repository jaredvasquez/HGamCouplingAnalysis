import os, sys, math
from sysTools import *
from ROOT import *
import HGamMoriondCats as HG

doPruning = 0

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
    #minsys = math.sqrt( (mcstat+1.0)**2 - mcstat**2 ) # Add at least 1.0% to the total uncertainty w/ stat uncert
    minsys = math.sqrt( (1.05*mcstat)**2 - mcstat**2 ) # Impact the total uncertainty by at least 5%

    syslist = sorted(list(allSysByCat[cat][proc]))
    for sys in syslist + [ syslist.pop(syslist.index('MCstat')) ]:
      if 'EG_' in sys: continue
      if 'Trig' in sys: continue
      #if not 'PH_' in sys: continue

      sigup, sigdn, form = allSysByCat[cat][proc][sys]
      sigup *= 100.; sigdn *= 100.

      tagasym = ''
      sigabs = (abs(sigup), abs(sigdn))

      # tag if +/- variations are more than 100% diff
      if ( min(sigabs) and max(sigabs)/min(sigabs) > 2 ): tagasym = '***'

      if (doPruning):
        # if one or more systs are zero, assume low stats and prune
        if not min(sigabs): continue

        # prune if +/- variations are more than 500% diff
        if ( max(sigabs)/min(sigabs) > 5 ): continue

        # if up/down are same sign variation, assume stat dominated and prune
        #if ((sigup*sigdn) > 0): continue

        # prune if does not add 5% relative change to total uncert
        if ( sys != 'MCstat' and max(sigabs) < minsys ): continue

      # Rename some things
      if (sys == 'MCstat'): sys = 'MCstat_%s_%s' % (cat,proc)
      else: sys = 'ATLAS_'+sys

      #if not float( '%.2f' % sum( map( abs, list(allSysByCat[cat][proc][sys]) ) ) ): continue
      print "    -->", sys, "+/- ( %+.2f %%, %+.2f %% )" % (sigup,sigdn), tagasym

print ""
