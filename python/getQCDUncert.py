from ROOT import *

binNames = ["GG2H_FWDH","GG2H_VBFTOPO_JET3VETO","GG2H_VBFTOPO_JET3","GG2H_0J","GG2H_1J_PTH_0_60","GG2H_1J_PTH_60_120","GG2H_1J_PTH_120_200","GG2H_1J_PTH_GT200","GG2H_GE2J_PTH_0_60","GG2H_GE2J_PTH_60_120","GG2H_GE2J_PTH_120_200","GG2H_GE2J_PTH_GT200"]

def printTable( table, header=None ):
  if (header):
    head = " | ".join(header)
    print head
    print "-"*len(head)
  for row in table: print " | ".join(row)

tf = TFile("output/Coupling_ggH/hist-ggH.root")

nCats, nBins = 10, len(binNames)

tables = {}
histName = "h2_catSTXS"
systNames = [ 'QCDyield', 'QCDres', 'QCDcut01', 'QCDcut12' ]

hnom = tf.Get( histName )

for systName in systNames:
  hSystName = "%s_%s" % (histName, systName)
  hsys = tf.Get( hSystName )
  hsys.Add( hnom, -1 )
  hsys.Divide( hnom )

  table = []
  for jbin in xrange(2, nBins+2):
    row = [ "%22s" % binNames[jbin-2] ]
    for ibin in xrange(1, nCats+1):
      sysUnc = hsys.GetBinContent( ibin, jbin )*100.
      row.append( "%+5.1f" % sysUnc )
    table.append( row )
  tables[systName] = table

  print ""
  print ""
  print "Systematic Table:", systName
  print ""
  header = [ "%22s" % "Truth Bin" ] + [ "%4d "%(icat+1) for icat in xrange(nCats) ]
  printTable( table, header=header )

print ""
