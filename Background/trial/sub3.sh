#!/bin/bash
cd /afs/hep.wisc.edu/user/kaho/CMSSW_10_2_16_UL/src/flashggFinalFit/Background
eval `scramv1 runtime -sh`
$CMSSW_BASE/src/flashggFinalFit/Background/bin/makeBkgPlots -f ggcat0,ggcat1,ggcat2,vbf -b CMS_Hemu_13TeV_multipdf.root -o trial/BkgPlots_cat3.root -d trial -c 3 -l "vbf" --sqrts 13  --intLumi 137.000000  --year 2016  --higgsResolution 1.0 -s ../../UWHiggs2016/em/results/Data2016JEC/AnalyzeEMSys/Signal.root --isMultiPdf --useBinnedData
