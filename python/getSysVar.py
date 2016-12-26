import os, sys
from sysTools import *
from ROOT import *
import HGamMoriondCats as HG

# Parameters
_NSYS = 120
procs = ["WH"]
strPATH = 'output/SystCoupling/SystCoupling_%s/hist-%s.root'


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
    sampName = "%s-%d" % (proc,i)
    filePATH = strPATH % (sampName,sampName)
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
for cat in HG.CatLabels:
  break
  print "\n",cat
  for proc in allSysByCat[cat]:
    print " ", proc
    for sys in allSysByCat[cat][proc]:
      print "      -->", sys
print ""

