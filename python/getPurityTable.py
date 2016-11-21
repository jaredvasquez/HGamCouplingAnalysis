from ROOT import *
import prettyplots

gStyle.SetOptStat(0)
prettyplots.setPalette("rainbow")

Stage1_Labels = ["ERROR","UNKNOWN","GG2H_FWDH","GG2H_VBFTOPO_JET3VETO","GG2H_VBFTOPO_JET3","GG2H_0J","GG2H_1J_PTH_0_60","GG2H_1J_PTH_60_120","GG2H_1J_PTH_120_200","GG2H_1J_PTH_GT200","GG2H_GE2J_PTH_0_60","GG2H_GE2J_PTH_60_120","GG2H_GE2J_PTH_120_200","GG2H_GE2J_PTH_GT200","VBF_QQ2HQQ_FWDH","QQ2HQQ_VBFTOPO_JET3VETO","QQ2HQQ_VBFTOPO_JET3","QQ2HQQ_VH2JET","QQ2HQQ_REST","QQ2HQQ_PTJET1_GT200","QQ2HLNU_FWDH","QQ2HLNU_PTV_0_150","QQ2HLNU_PTV_150_250_0J","QQ2HLNU_PTV_150_250_GE1J","QQ2HLNU_PTV_GT250","QQ2HLL_FWDH","QQ2HLL_PTV_0_150","QQ2HLL_PTV_150_250_0J","QQ2HLL_PTV_150_250_GE1J","QQ2HLL_PTV_GT250"]
#,"GG2HLL_FWDH","GG2HLL_PTV_0_150","GG2HLL_PTV_GT150_0J","GG2HLL_PTV_GT150_GE1J","TTH_FWDH","TTH","BBH_FWDH","BBH","TH_FWDH","TH"]

HCombStage1_Labels = ["ERROR","UNKNOWN","GG2H_FWDH","GG2H_VBFTOPO_JET3VETO","GG2H_VBFTOPO_JET3","GG2H_0J","GG2H_1J_PTH_0_60","GG2H_1J_PTH_60_120","GG2H_1J_PTH_120_200","GG2H_1J_PTH_GT200","GG2H_GE2J_PTH_0_60","GG2H_GE2J_PTH_60_120","GG2H_GE2J_PTH_120_200","GG2H_GE2J_PTH_GT200","VBF_QQ2HQQ_FWDH","VBF_QQ2HQQ_VBFTOPO_JET3VETO","VBF_QQ2HQQ_VBFTOPO_JET3","VBF_QQ2HQQ_VH2JET","VBF_QQ2HQQ_REST","VBF_QQ2HQQ_PTJET1_GT200","WH_QQ2HQQ_FWDH","WH_QQ2HQQ_VBFTOPO_JET3VETO","WH_QQ2HQQ_VBFTOPO_JET3","WH_QQ2HQQ_VH2JET","WH_QQ2HQQ_REST","WH_QQ2HQQ_PTJET1_GT200","ZH_QQ2HQQ_FWDH","ZH_QQ2HQQ_VBFTOPO_JET3VETO","ZH_QQ2HQQ_VBFTOPO_JET3","ZH_QQ2HQQ_VH2JET","ZH_QQ2HQQ_REST","ZH_QQ2HQQ_PTJET1_GT200","QQ2HLNU_FWDH","QQ2HLNU_PTV_0_150","QQ2HLNU_PTV_150_250_0J","QQ2HLNU_PTV_150_250_GE1J","QQ2HLNU_PTV_GT250","QQ2HLL_FWDH","QQ2HLL_PTV_0_150","QQ2HLL_PTV_150_250_0J","QQ2HLL_PTV_150_250_GE1J","QQ2HLL_PTV_GT250","GG2HLL_FWDH","GG2HLL_PTV_0_150","GG2HLL_PTV_GT150_0J","GG2HLL_PTV_GT150_GE1J","TTH_FWDH","TTH","BBH_FWDH","BBH","THQB_FWDH","THQB","TWH_FWDH","TWH"]

CatLabels = ["NONE","Rest, ==0j (Cen)","Rest, ==0j (Fwd)","Rest, ==1j, pT [0,60]","Rest, ==1j, pT [60,120]","Rest, ==1j, pT [120,200]","Rest, ==1j, pT [200,inf]","Rest, >=2j, pT [0,60]","Rest, >=2j, pT [60,120]","Rest, >=2j, pT [120,200]","Rest, >=2j, pT [200,inf]","VBF loose, LOW", "VBF loose, HIGH", "VBF tight, LOW","VBF tight, HIGH", "VH had loose","VH had tight","VBF BSM","VHMET, pT [0,150]","VHMET, pT [150,250]","VHMET, pT [250,inf]","VHlep, pT [0,150]","VHlep, pT [150,inf]","VHdilep, pT [0,150]","VHdilep, pT [150,inf]", "ttHhad_4j1b","ttHhad_4j2b","ttHhad_5j1b","ttHhad_5j2b","ttHhad_6j1b","ttHhad_6j2b","tHlep_0fwd","tHlep_1fwd","ttHlep"]

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
  hsum = tfs[0].Get( histName )
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
  return decorateHist( hPurity, "Truth Purity / Reco Category" )


histName = "h2_catSTXS"
procs = ["ggH","VBF","WH","ZH"]
tfs = [ TFile("output/Coupling_%s/hist-%s.root" % (p,p)) for p in procs ]
hsum = sumHist( histName, tfs )
hpur = purityHistX( hsum )

can = TCanvas(); can.cd()
can.SetTopMargin(0.02)
can.SetRightMargin(0.13)
can.SetLeftMargin(0.28)
can.SetBottomMargin(0.25)

hpur.Draw("colz")

can.SaveAs("plots/purity.pdf")
can.SaveAs("plots/purity.png")
