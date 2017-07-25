import os
import sys
import json
import tabulate
import numpy as np
from ROOT import *

tabulate.LATEX_ESCAPE_RULES={}
sys.dont_write_bytecode = True
import HGamMoriondCatsBDT as HG

#binNames = ['gg2H_fwdH','gg2H_VBFtopo_jet3veto','gg2H_VBFtopo_jet3','gg2H_0J','gg2H_1J_ptH_0_60','gg2H_1J_ptH_60_120','gg2H_1J_ptH_120_200','gg2H_1J_ptH_gt200','gg2H_ge2J_ptH_0_60','gg2H_ge2J_ptH_60_120','gg2H_ge2J_ptH_120_200','gg2H_ge2J_ptH_gt200']
from STXSFineIndex import FineIndexLabels as binNames
#binNames = binNames[1:]

for ibin, binName in enumerate(binNames):
  print ibin, binName
print ''
print ''



mergeBinSet_WEAK = [
                [2, 3],   # WEAK MERGING
                [14, 15, 20, 21, 26, 27],
                [16, 17, 22, 23, 28, 29],
                [18, 24, 30],
                [32, 33, 34, 35, 37, 38, 39, 40, 42, 43, 44],
              ]
mergeBinSet_STRONG = [
                [8, 12],  # STRONG MERGING
                [2, 3, 9, 10, 11],
                list(range(14,19) + range(20,25) + range(26,31)),
                list(range(32,36) + range(37,41) + range(42,45)),
              ]
mergeBinSet_ALL = [ range(1,len(binNames)) ] # uncomment for production XS uncerts


mergeBinSet = mergeBinSet_WEAK
outputfile = 'AccQCDScale_WEAK.json'

#mergeBinSet = mergeBinSet_STRONG
#outputfile = 'AccQCDScale_STRONG.json'

mergeBinSet = mergeBinSet_ALL
outputfile = 'AccQCDScale_prodXS.json'




print 'Merging bins:'
for mergeBins in mergeBinSet:
  print '-'*30
  for ibin in mergeBins:
    print '   -', binNames[ibin]
print '-'*30


fmt = '%5.1f'
def formatCell(sys):
  if float(fmt % sys) == 0:
    return '---'
  return fmt % sys

def formatSys(sys):
  return fmt % sys

def printTable( table, header=None ):
  if (header):
    head = ' | '.join(header)
    print head
    print '-'*len(head)
  for row in table: print ' | '.join(row)

def mergeBinsH2( h2 ):
  for mergeBins in mergeBinSet:
    for ibin in xrange(0,h2.GetNbinsX()+1):
      sumBins = sum( [ h2.GetBinContent(ibin, im+2) for im in mergeBins ] )
      #print ibin
      for jbin in mergeBins:
        #print '  -', jbin, h2.GetBinContent( ibin, jbin+2 ), sumBins
        h2.SetBinContent( ibin, jbin+2, sumBins )
  return h2

def getAccHist( h2 ):
  h2 = mergeBinsH2( h2 )
  for ibin in xrange(1,h2.GetNbinsX()+1):
    for jbin in xrange(1,h2.GetNbinsY()+1):
      err = h2.GetBinContent(ibin,jbin)
      if (err == 0.): continue
      h2.SetBinContent( ibin, jbin, err / h2.GetBinContent(0,jbin) )
  return h2

procs = 'VBF_NNPDF WH_NLO ZH_NLO ggZH'.split()
NPproc = {'VBF_NNPDF': 'qqH', 'WH_NLO': 'VH', 'ZH_NLO': 'VH', 'ggZH': 'VH'}
#procs = 'VBF_NNPDF'.split()
#procs = 'WH_NLO'.split()

systMap = { c: {} for c in HG.CatLabels }

for proc in procs:
  tf = TFile('output/HGamPDF_{0}/hist-{0}.root'.format(proc))

  nCats, nBins = 31, len(binNames)

  tables = {}
  tables2 = {}
  #histName = 'h2_catSTXS'
  histName = 'h2_fineIndex'

  #systMap = { cat: { tbin: {} for tbin in binNames } for cat in HG.CatLabels }

  hnom = tf.Get( histName )
  hnom = mergeBinsH2( hnom )
  hnom = getAccHist( hnom ) # uncomment for acceptance uncerts

  NQCD = 8
  if 'VBF' in proc: NQCD = 6
  maps = []

  for isys in xrange(NQCD):
    hSystName = '%s_QCD%d' % (histName, isys)
    hsys = tf.Get( hSystName )
    hsys = mergeBinsH2( hsys )
    hsys = getAccHist( hsys ) # uncomment for acceptance uncerts

    hsys.Add( hnom, -1 )
    hsys.Divide( hnom )

    sysmap = [[hsys.GetBinContent(ibin+1, jbin+1) for ibin in xrange(hsys.GetNbinsX()+1)] for jbin in xrange(hsys.GetNbinsY()+1)]
    sysmap = np.array( sysmap )
    maps.append(sysmap)

  sys3D = np.dstack(tuple(maps))
  sys2D = np.max(np.abs(sys3D), axis=2)

  print sys3D.shape, '-->', sys2D.shape
  print len(binNames)

  for j, catName in enumerate(HG.CatLabels):
    for i, binName in enumerate(binNames):
      if 'Unknown' in binName: continue
      if '_fwdH' in binName: continue
      uncert = round(sys2D[i][j],4)
      if (uncert < 0.005): continue
      if binName not in systMap[catName]:
        systMap[catName][binName] = {}
      systMap[catName][binName]['ATLAS_QCDscale_'+NPproc[proc]] = [uncert, uncert, 'logn']
      print '%25s - %35s  :  %5s  %4.1f' % (catName, binName, NPproc[proc], uncert*100.)
    print ''

Nuncerts = 0
for c in systMap:
  for b in systMap[c]:
    Nuncerts += 1
    #print c, b, systMap[c][b]
print '# of uncerts :', Nuncerts
json.dump(systMap, open(outputfile,'wb'))

  #sys2D = np.amax(sys3D, axis=2, key='abs')

    #table = []
    #table2 = []
    #for ibin, catName in enumerate(HG.CatLabels):
    #  row = [ '%22s' % catName.replace('_','\_') ]
    #  row2 = [ '%22s' % catName.replace('_','\_') ]
    #  for jbin, binName in enumerate(binNames):
    #    sysUnc = hsys.GetBinContent( ibin+1, jbin+3 )*100.
    #    row.append( formatCell(sysUnc) )
    #    row2.append( formatSys(sysUnc) )
    #  table.append(row)
    #  table2.append(row2)
    #tables[systName] = table

    #for ibin, catName in enumerate(HG.CatLabels):
    #  for jbin, binName in enumerate(binNames):
    #    err = float('%.3f' % hsys.GetBinContent( ibin+1, jbin+3 ))
    #    if (err==0.): continue
    #    systMap[catName][binName]['ATLAS_QCDscale_ggF_'+systName] = [ err, err, 'logn' ]

    #print ''
    #print ''
    #print 'Systematic Table:', systName
    #print ''
    #header = [ '%22s' % 'Truth Bin' ] + [ '%4d '%(icat+1) for icat in xrange(nBins) ]
    #printTable( table2, header=header )

    #log = open('tablesQCD/table_%s.tex' % systName,'w+')
    #print >> log, tabulate.tabulate( table, headers=header, tablefmt='latex' )
    #print tabulate.tabulate( table, headers=header, tablefmt='latex' )


  #print ''
  #table = []
  #for ibin, binName in enumerate(binNames):
  #  print ' %2d : %s ' % (ibin+1, binName)
  #  table.append( [ibin+1,binName.replace('_','\_')] )
  #log = open('tablesQCD/table_Legend.tex', 'w+')
  #print >> log, tabulate.tabulate( table, headers=['Index','Truth Bin'], tablefmt='latex' )

  ## clean up ggF uncerts
  #for catName in HG.CatLabels:
  #  for binName in binNames:
  #    if len(systMap[catName][binName]) == 0:
  #      del systMap[catName][binName]
  #  if len(systMap[catName]) == 0:
  #    del systMap[catName]
  #json.dump(systMap, open('ggFQCD.json','wb'))
  #
  #print ''
