from ROOT import *
import prettyplots
import HGamMoriondCats as HG

gStyle.SetOptStat(0)
prettyplots.setPalette("rainbow")

CatLabels = ["NONE","Rest, ==0j (Cen)","Rest, ==0j (Fwd)","Rest, ==1j, pT [0,60]","Rest, ==1j, pT [60,120]","Rest, ==1j, pT [120,200]","Rest, ==1j, pT [200,inf]","Rest, >=2j, pT [0,60]","Rest, >=2j, pT [60,120]","Rest, >=2j, pT [120,200]","Rest, >=2j, pT [200,inf]","VBF LOW, loose", "VBF LOW, tight", "VBF HIGH, loose", "VBF HIGH, tight", "VH had loose", "VHhad tight","qqH BSM","VHMET, pT [0,150]","VHMET, pT [150,250]","VHMET, pT [250,inf]","VHlep, pT [0,150]","VHlep, pT [150,inf]","VHdilep, pT [0,150]","VHdilep, pT [150,inf]","ttHhad_6j2b","ttHhad_6j1b","ttHhad_5j2b","ttHhad_5j1b","ttHhad_4j2b","ttHhad_4j1b","ttHlep","tHlep_1fwd","tHlep_0fwd"]

def decorateHist( hist, ztitle=None ):
  for ibin, catName in enumerate(CatLabels):      hist.GetXaxis().SetBinLabel( ibin, catName )
  for jbin, binName in enumerate(Stage1_Labels):  hist.GetYaxis().SetBinLabel( jbin, binName )
  hist.GetXaxis().LabelsOption("v")
  hist.GetXaxis().SetTitle("Reco Category")
  hist.GetYaxis().SetTitle("STXS Truth Bin")
  if (ztitle): hist.GetZaxis().SetTitle(ztitle)
  hist.GetXaxis().SetTitleOffset(4.0)
  hist.GetYaxis().SetTitleOffset(4.5)
  hist.SetMinimum(0.0)
  hist.GetXaxis().SetRangeUser(  0.5, len(CatLabels)+0.5     )
  hist.GetYaxis().SetRangeUser( -0.5, len(Stage1_Labels)-0.5 )
  return hist

def sumHist( histName, tfs ):
  hsum = tfs[0].Get( histName ).Clone()
  for i in xrange(1,len(tfs)):
    hsum.Add( tfs[i].Get( histName ) )
  return hsum

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

color = {}
color[ 'ggH'] = 416
color[ 'VBF'] = 418
color[  'WH'] = 432
color[  'ZH'] = 600
color[ 'ttH'] = 880
color[ 'bbH'] = 632
color[ 'tWH'] = 800
color['tHjb'] = 401

histName = "h_catSTXS"
procs = ["ggH","VBF","WH","ZH","ttH","tHjb","tWH","bbH"]
tfs = [ TFile("output/Coupling_%s/hist-%s.root" % (p,p)) for p in procs ]
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

can = TCanvas('can','can',800,1600); can.cd()
can.SetTopMargin(0.08)
can.SetRightMargin(0.02)
can.SetLeftMargin(0.20)
can.SetBottomMargin(0.08)

hsum.GetYaxis().SetTitle("Fraction of Signal Process / Category")
hsum.GetYaxis().CenterTitle()
for ibin, catName in enumerate(HG.CatLabels,1):
  hsum.GetXaxis().SetBinLabel(ibin, catName)
hsum.SetMaximum(1.0)
hsum.GetXaxis().SetRangeUser(1,len(HG.CatLabels))
hsum.Draw("HIST HBAR")
hstack.Draw("HBAR SAME")

line = TLine(); line.SetLineWidth(2)
for i in xrange(len(HG.CatLabels)):
  line.DrawLine(0,i+.5,1,i+.5);

tl = TLatex()
tl.SetNDC()
tl.DrawLatex(0.22, 0.95, "#bf{#it{#bf{ATLAS}} Internal}")
tl.SetTextSize(0.032)
tl.DrawLatex(0.55, 0.95, "#bf{#it{H #rightarrow #gamma#gamma,  m_{H} = 125.09} GeV}")

gPad.RedrawAxis()
gPad.SetTicks(1,1)

can.SaveAs("plots/purity1d.pdf")
can.SaveAs("plots/purity1d.png")
