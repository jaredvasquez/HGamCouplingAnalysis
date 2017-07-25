import math
import ROOT
import tabulate
import numpy as np
import HGamMoriondCatsBDT2 as HG

#--------------------------------------------------------------------
def fixName(p):
  p = p.replace('_NNLOPS','')
  p = p.replace('_NNPDF', '')
  p = p.replace('_NLO',   '')
  return '$%s$' % p

def sigfigs(num, sig_figs):
  if num == 0: return 0
  num = round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
  if (num > 10): return int(num)
  return num

def myround(num, digs):
  num = round(num, digs)
  #if num >= 10: num = int(num)
  if num == 0: return '---'
  return num

def printTable(tab, headers=None):
  ncols = len(tab[0])
  print '\\begin{tabular}{r|c|%s}' % ('r'*(ncols-2))
  print ' & & \multicolumn{9}{c}{Composition [\%]} \\\\'
  tmp = '%32s' + ' & %6s'*(ncols-1) + ' \\\\'
  if headers:
    print tmp % tuple(headers)
    print '\hline'
  for row in tab:
    print tmp % tuple(row)
  print '\\end{tabular}'

#--------------------------------------------------------------------
histName = 'h_catSTXS'
procs = ['ggH_NNLOPS','VBF_NNPDF','WH_NLO','ZH_NLO','ggZH','ttH','bbH','tHjb','tHW']
tfs = [ ROOT.TFile('output/HGamCoupling_%s/hist-%s.root' % (p,p)) for p in procs ]
hs  = [ tf.Get(histName) for tf in tfs ]

table = []
for ibin, catName in enumerate(HG.CatLabels):
  Ns = np.array([ h.GetBinContent(ibin+1) for h in hs ]) * 36.067
  Ns = np.array([ N if N > 0 else 0 for N in Ns ])
  Nsum = np.sum(Ns)
  fracs = 100 * Ns / Nsum
  #row = [ catName.replace('_','\_'), round(Nsum, 3) ] + [ round(f, 1) for f in fracs ]
  nsig = 2
  #if (Nsum > 1): nsig = 3
  row = [ catName.replace('_','\_'), sigfigs(Nsum, nsig) ] + [ myround(f, 1) for f in fracs ]
  table.append( row )

headers = ['Category', '$N_H$'] + [fixName(p) for p in procs]
#print tabulate.tabulate( table, headers=headers )
printTable(table, headers=headers)
