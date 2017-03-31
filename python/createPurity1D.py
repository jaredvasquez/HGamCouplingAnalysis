from ROOT import *
import prettyplots
import HGamMoriondCatsBDT as HG

gStyle.SetOptStat(0)
#prettyplots.setPalette("rainbow")

#CatLabels = ["NONE","Rest, ==0j (Cen)","Rest, ==0j (Fwd)","Rest, ==1j, pT [0,60]","Rest, ==1j, pT [60,120]","Rest, ==1j, pT [120,200]","Rest, ==1j, pT [200,inf]","Rest, >=2j, pT [0,60]","Rest, >=2j, pT [60,120]","Rest, >=2j, pT [120,200]","Rest, >=2j, pT [200,inf]","VBF LOW, loose", "VBF LOW, tight", "VBF HIGH, loose", "VBF HIGH, tight", "VH had loose", "VHhad tight","qqH BSM","VHMET, pT [0,150]","VHMET, pT [150,250]","VHMET, pT [250,inf]","VHlep, pT [0,150]","VHlep, pT [150,inf]","VHdilep, pT [0,150]","VHdilep, pT [150,inf]","ttHhad_6j2b","ttHhad_6j1b","ttHhad_5j2b","ttHhad_5j1b","ttHhad_4j2b","ttHhad_4j1b","ttHlep","tHlep_1fwd","tHlep_0fwd"]

def DrawNDC(self, x, y, text): self.DrawLatexNDC( x, y, '#bf{%s}' % text )
TLatex.DrawNDC = DrawNDC

#--------------------------------------------------------------------
def decorateHist( hist, ztitle=None ):
  for ibin, catName in enumerate(CatLabels):
    hist.GetXaxis().SetBinLabel( ibin, catName )
  for jbin, binName in enumerate(Stage1_Labels):
    hist.GetYaxis().SetBinLabel( jbin, binName )
  hist.GetXaxis().LabelsOption("v")
  hist.GetXaxis().SetTitle("Reco Category")
  hist.GetYaxis().SetTitle("STXS Truth Bin")
  if (ztitle): hist.GetZaxis().SetTitle(ztitle)
  hist.GetXaxis().SetTitleOffset(4.5)
  hist.GetYaxis().SetTitleOffset(4.5)
  hist.SetMinimum(0.0)
  hist.GetXaxis().SetRangeUser(  0.5, len(CatLabels)+0.5     )
  hist.GetYaxis().SetRangeUser( -0.5, len(Stage1_Labels)-0.5 )
  return hist

#--------------------------------------------------------------------
def sumHist( histName, tfs ):
  hsum = tfs[0].Get( histName ).Clone()
  for i in xrange(1,len(tfs)):
    hsum.Add( tfs[i].Get( histName ) )
  return hsum

#--------------------------------------------------------------------
def purityHistX( hist ):
  hPurity = hist.Clone()
  for ibin in xrange( 1, hist.GetNbinsX() ):
    NX = 0.0
    for jbin in xrange( 1, hist.GetNbinsY() ):
      NX += hPurity.GetBinContent( ibin, jbin )
    if (NX == 0.): continue
    for jbin in xrange( 1, hist.GetNbinsY() ):
      hPurity.SetBinContent( ibin, jbin, hPurity.GetBinContent( ibin, jbin )/float(NX) )
  #return decorateHist( hPurity, "Truth Purity / Reco Category" )
  return hPurity

#--------------------------------------------------------------------
color = {}
color[ 'ggH'] = 416
color[ 'VBF'] = 418
color[  'WH'] = 432
color[  'ZH'] = 600
color[ 'ttH'] = 880
color[ 'bbH'] = 632
color[ 'tHW'] = 800
color['tHjb'] = 401
color['ggH_NNLOPS'] = 416
color[ 'VBF_NNPDF'] = 418
color[    'WH_NLO'] = 432
color[    'ZH_NLO'] = 600


#--------------------------------------------------------------------
histName = "h_catSTXS"
#procs = ["ggH","VBF","WH","ZH","ttH","tHjb","tWH","bbH"]
procs = ["ggH_NNLOPS","VBF_NNPDF","WH_NLO","ZH_NLO","ttH","bbH","tHjb","tHW"]
tfs = [ TFile("output/HGamCoupling_%s/hist-%s.root" % (p,p)) for p in procs ]
hs  = [ tf.Get(histName) for tf in tfs ]
hsum = sumHist( histName, tfs )

for ibin in xrange(1,hsum.GetNbinsX()+1):
  sum = hsum.GetBinContent(ibin)
  if (sum==0.): continue
  for h in hs: h.SetBinContent( ibin, h.GetBinContent(ibin)/float(sum) )

hstack = THStack("hstack","")
for i, h in enumerate(hs):
  h.SetFillColor(color[procs[i]])
  hstack.Add(h)

can = TCanvas('can','can',800,1200); can.cd()
can.SetTopMargin(0.12)
can.SetRightMargin(0.02)
can.SetLeftMargin(0.20)
can.SetBottomMargin(0.10)

hsum.GetYaxis().SetTitle("Fraction of Signal Process / Category")
hsum.GetYaxis().CenterTitle()
for ibin, catName in enumerate(HG.CatLabels,1):
  hsum.GetXaxis().SetBinLabel(ibin, catName)
hsum.SetMaximum(1.0)
hsum.GetXaxis().SetRangeUser(1,len(HG.CatLabels))
#hsum.GetYaxis().SetLabelOffset(0.0005)
hsum.GetYaxis().SetTitleOffset(1.1)
hsum.Draw("HIST HBAR")
hstack.Draw("HBAR SAME")

line = TLine(); line.SetLineWidth(2)
for i in xrange(len(HG.CatLabels)):
  line.DrawLine(0,i+.5,1,i+.5);

tl = TLatex()
tl.SetNDC()
#tl.DrawLatex(0.22, 0.90, "#bf{#it{#bf{ATLAS}} Internal}")
tl.DrawNDC(0.21, 0.90, "#it{#bf{ATLAS}} Internal")
tl.SetTextSize(0.032)
tl.DrawNDC(0.57, 0.90, "#it{H #rightarrow #gamma#gamma,  m_{H} = 125.09} GeV")

paves = []
def drawLabel(x, y, proc):
  tp = TPave()
  tp.SetFillColor(color[proc])
  tp.SetX1NDC(x)
  tp.SetY1NDC(y)
  tp.SetX2NDC(x+0.04)
  tp.SetY2NDC(y+0.03)
  tp.Draw('NDC')
  paves.append(tp)
  tl.DrawNDC(x+0.05, y+0.007, proc)

drawLabel(0.04,0.95,'ggH')
drawLabel(0.16,0.95,'VBF')
drawLabel(0.28,0.95,'WH')
drawLabel(0.40,0.95,'ZH')
drawLabel(0.51,0.95,'ttH')
drawLabel(0.62,0.95,'bbH')
drawLabel(0.75,0.95,'tHjb')
drawLabel(0.87,0.95,'tHW')

gPad.RedrawAxis()
gPad.SetTicks(1,1)

to = TFile("stack.root","RECREATE")
hstack.Write("purity1D")

can.SaveAs("plots/purity1d.pdf")
can.SaveAs("plots/purity1d.png")
