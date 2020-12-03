import ROOT
import math
file_ = ROOT.TFile("CMS_Hemu_13TeV_multipdf.root ")
ws = file_.Get("multipdf")
bkg = ws.var("env_pdf_vbf_bern1_p0")

bkg.printValue(ROOT.cout)
