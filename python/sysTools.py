import os
from sys import exit
from ROOT import *
import HGamMoriondCats as HG

suffixes = { '__1up': '__1down', '_ScaleUp': '_ScaleDown' }

def getStat( hnom ):
  diff = [ abs(hnom.GetBinError(i+1)/hnom.GetBinContent(i+1)) if hnom.GetBinContent(i+1) else 0. for i in xrange(len(HG.CatLabels)) ]
  return diff

def getDiff( hnomOrig, hsysOrig ):
  hsys = hsysOrig.Clone()
  hsys.Add( hnomOrig, -1 )
  hsys.Divide( hsysOrig )
  diff = [ hsys.GetBinContent(i+1) for i in xrange(len(HG.CatLabels)) ]
  return diff


def groupNPs( sysByCat ):
  newSysByCat = { cat: {} for cat in sysByCat }
  for cat in sysByCat:
    systs = sysByCat[cat]
    for sysName in systs:
      for upsuf in suffixes:
        dnsuf = suffixes[upsuf]
        if upsuf in sysName:
          npName = sysName.replace(upsuf,'')
          if not (npName+dnsuf) in systs:  # has up suffix but no down variation, e.g. JER
            newSysByCat[cat][npName] = ( systs[npName+upsuf], -1.0*systs[npName+upsuf], 'logn' )
          else:
            newSysByCat[cat][npName] = ( systs[npName+upsuf], systs[npName+dnsuf], 'asym' )
          break

        elif dnsuf in sysName: # ignore _down suffixes, after quick check
          npName = sysName.replace(dnsuf,'')
          if not (npName+upsuf) in systs:
            print "WARNING: %s exists but %s does not" % (npName+dnsuf, npName+upsuf)
          break

      else: # Symmetric systematics
        newSysByCat[cat][sysName] = ( systs[sysName], -1.0*systs[sysName], 'logn' )

  return newSysByCat


def getVariation( tfs, histName='h_catSTXS' ):
  systVars = {}
  hnom = tfs[0].Get(histName)
  systVars['MCstat'] = getStat( hnom )
  for ifile in xrange(1,len(tfs)):
    tf = tfs[ifile]
    for k in tf.GetListOfKeys():
      keyName = k.GetName()
      if histName+'_' in keyName:
        sysName = keyName.replace(histName+'_','')
        hsys = tfs[ifile].Get(keyName)
        systVars[sysName] = getDiff( hnom, hsys )
        break
    else:
      print "*** Never found systematic histogram"

  # Organize systs by category
  sysByCat = { cat: { sys: systVars[sys][i] for sys in systVars } for i, cat in enumerate(HG.CatLabels) }
  sysByCat = groupNPs( sysByCat )
  return sysByCat
