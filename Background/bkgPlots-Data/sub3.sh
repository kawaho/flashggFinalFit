#!/bin/bash
cd /afs/hep.wisc.edu/user/kaho/CMSSW_10_2_16_UL/src/flashggFinalFit/Background
eval `scramv1 runtime -sh`
$CMSSW_BASE/src/flashggFinalFit/Background/bin/makeBkgPlots -f vbf,gg0,gg1,gg2 -b example.root -o bkgPlots-Data/BkgPlots_cat3.root -d bkgPlots-Data -c 3 -l "gg2" --sqrts 13  --intLumi 12.900000  --year 2016  --doBands --massStep 1.000 --nllTolerance 0.050 -L 160 -H 110 --higgsResolution 1.0 --isMultiPdf --useBinnedData
