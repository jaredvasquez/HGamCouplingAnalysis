import os
from sysTools import *
from sys import exit
from ROOT import *
import HGamMoriondCats as HG

# Parameters
_NSYS = 120
procs = ["HHDMmr260mx60br90_DMgamgam"]
strPATH = "output/SystTetst_%s/hist-%s.root"


# List all Categories
listCats = False
if (listCats):
  print "Categories (%d Total):" % (len(HG.CatLabels))
  for label in HG.CatLabels:
    print "  ", label
  print ""


# Get Systematics Variations for each process
allSys = {}
for proc in procs:
  tfs = []
  for i in xrange(0,_NSYS+1):
    sampName = "%s-%d" % (proc,i)
    filePATH = strPATH % (sampName,sampName)
    if not os.path.isfile(filePATH):
      print "ERROR: Could not find", sampName
      sys.exit()
    tfs.append( TFile(filePATH) )
  allSys[proc] = getVariation(tfs)


# Organize All Systematics By: ( Category | Proc | Systematic )
allSysByCat = { cat: { proc: allSys[proc][cat] for proc in procs } for cat in HG.CatLabels }


# Print results
for cat in HG.CatLabels:
  print "\n",cat
  for proc in allSysByCat[cat]:
    print " ", proc
    for sys in allSysByCat[cat][proc]:
      print "      -->", sys
print ""

