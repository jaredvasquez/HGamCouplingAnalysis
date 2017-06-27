import os
import sys
import json
import tabulate
from ROOT import *

tabulate.LATEX_ESCAPE_RULES={}
sys.dont_write_bytecode = True
import HGamMoriondCatsBDT as HG

binNames = ['gg2H_fwdH','gg2H_VBFtopo_jet3veto','gg2H_VBFtopo_jet3','gg2H_0J','gg2H_1J_ptH_0_60','gg2H_1J_ptH_60_120','gg2H_1J_ptH_120_200','gg2H_1J_ptH_gt200','gg2H_ge2J_ptH_0_60','gg2H_ge2J_ptH_60_120','gg2H_ge2J_ptH_120_200','gg2H_ge2J_ptH_gt200']
binNames = binNames[1:]

mergeBins = [ 1, 2 ]
#mergeBins = range(1,len(binNames)) # uncomment for production XS uncerts
print 'Merging bins...'
for ibin in mergeBins:
  print '  -', binNames[ibin-1]


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
  for ibin in xrange(0,h2.GetNbinsX()+1):
    sumBins = sum( [ h2.GetBinContent(ibin, im+2) for im in mergeBins ] )
    #print ibin
    for jbin in mergeBins:
      #print '  -', jbin, h2.GetBinContent( ibin, jbin+2 ), sumBins
      h2.SetBinContent( ibin, jbin+2, sumBins )
  return h2

def getAccHist(h2):
  h2 = mergeBinsH2(h2)
  for ibin in xrange(1, h2.GetNbinsX()+1):
    for jbin in xrange(1, h2.GetNbinsY()+1):
      err = h2.GetBinContent(ibin, jbin)
      if (err == 0.): continue
      h2.SetBinContent(ibin, jbin, err/h2.GetBinContent(0,jbin))
  return h2

tf = TFile('output/HGamPDF_ggH_NNLOPS/hist-ggH_NNLOPS.root')
os.system('mkdir -p tablesQCD')

nCats, nBins = 31, len(binNames)

tables = {}
tables2 = {}
histName = 'h2_catSTXS'
systNames = [ 'mu', 'res', 'qm_t', 'pTH60', 'pTH120', 'mig01', 'mig12', 'vbf2j', 'vbf3j' ]

systMap = { cat: { tbin: {} for tbin in binNames } for cat in HG.CatLabels }

hnom = tf.Get( histName )
#hnom = getAccHist( hnom ) # uncomment for acceptance uncerts
for systName in systNames:
  hSystName = '%s_QCD_2017_%s' % (histName, systName)

  hsys = tf.Get( hSystName )
  #hsys = getAccHist( hsys ) # uncomment for acceptance uncerts
  hsys.Add( hnom, -1 )
  hsys.Divide( hnom )

  table = []
  table2 = []
  for ibin, catName in enumerate(HG.CatLabels):
    row = [ '%22s' % catName.replace('_','\_') ]
    row2 = [ '%22s' % catName.replace('_','\_') ]
    for jbin, binName in enumerate(binNames):
      sysUnc = hsys.GetBinContent( ibin+1, jbin+3 )*100.
      row.append( formatCell(sysUnc) )
      row2.append( formatSys(sysUnc) )
    table.append(row)
    table2.append(row2)
  tables[systName] = table

  for ibin, catName in enumerate(HG.CatLabels):
    for jbin, binName in enumerate(binNames):
      err = float('%.3f' % hsys.GetBinContent( ibin+1, jbin+3 ))
      if (err==0.): continue
      systMap[catName][binName]['ATLAS_QCDscale_ggF_'+systName] = [ err, err, 'logn' ]

  print ''
  print ''
  print 'Systematic Table:', systName
  print ''
  header = [ '%22s' % 'Truth Bin' ] + [ '%4d '%(icat+1) for icat in xrange(nBins) ]
  printTable( table2, header=header )

  log = open('tablesQCD/table_%s.tex' % systName,'w+')
  print >> log, tabulate.tabulate( table, headers=header, tablefmt='latex' )


print ''
table = []
for ibin, binName in enumerate(binNames):
  print ' %2d : %s ' % (ibin+1, binName)
  table.append( [ibin+1,binName.replace('_','\_')] )
log = open('tablesQCD/table_Legend.tex', 'w+')
print >> log, tabulate.tabulate( table, headers=['Index','Truth Bin'], tablefmt='latex' )


# clean up ggF uncerts
for catName in HG.CatLabels:
  for binName in binNames:
    if len(systMap[catName][binName]) == 0:
      del systMap[catName][binName]
  if len(systMap[catName]) == 0:
    del systMap[catName]
json.dump(systMap, open('ggFQCD.json','wb'))

print ''
