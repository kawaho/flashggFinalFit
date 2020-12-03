#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <typeinfo>

#include "TROOT.h"
#include <TStyle.h>
#include "TFile.h"
#include "TMath.h"
#include "TH2.h"
#include "RooCBShape.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TKey.h"
#include "TMacro.h"
#include "TClass.h"
#include "TIterator.h"
#include "TRandom3.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"
#include "TCanvas.h"

#include "../interface/SimultaneousFit.h"
#include "../interface/InitialFit.h"
#include "../interface/LinearInterp.h"
#include "../interface/FinalModelConstruction.h"
#include "../interface/Packager.h"
#include "../interface/WSTFileWrapper.h"
#include "../interface/ReplacementMap.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = boost::program_options;

typedef map<int,map<string,RooRealVar*> > parlist_t;

string flashggCatsStr_ = "gg0,gg1,gg2,gg3,vbf";//,0JetEE,1JetEB-MB,1JetEB-ME,1JetEE,2JetEB-MB,2JetEB-ME,2JetEE,2JetVBF";
vector<string> flashggCats_;

int  nBins_;
vector<int> massList_;
int mh = 125;
int mhLow_=115;
int mhHigh_=135;
vector<int> skipMasses_;
RooRealVar *mass_;

vector<string> map_cat_;

int ncpu_=1;

RooAbsPdf* build_CBpGaus_Pdf(string cat_, parlist_t fitParams, RooRealVar *mass, RooRealVar *MH, bool verbosity_){
 
  string catname;
  catname = Form("%s",cat_.c_str());

  map<string,RooRealVar*> fitParams_125 = fitParams[125];

  RooAbsReal *a1_cb = fitParams_125["a1_mh125_cb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter a1_cb"  << " ? " << a1_cb << std::endl;
  a1_cb->SetName(Form("a1_%s",cat_.c_str()));

  RooAbsReal *n1_cb = fitParams_125["n1_mh125_cb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter n1_cb"  << " ? " << n1_cb << std::endl;
  n1_cb->SetName(Form("n1_%s",cat_.c_str()));

  RooAbsReal *dm_cb = fitParams_125["dm_mh125_cb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter dm_cb"  << " ? " << dm_cb << std::endl;
  dm_cb->SetName(Form("dm_%s",cat_.c_str()));

  RooAbsReal *mean_cb = new RooFormulaVar(Form("mean_cb_%s",cat_.c_str()),Form("mean_cb_%s",cat_.c_str()),"@0+@1",RooArgList(*MH,*dm_cb));

  RooAbsReal *sig_fit_cb = fitParams_125["sigma_mh125_cb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter sigma_cb"  << " ? " << sig_fit_cb << std::endl;
  sig_fit_cb->SetName(Form("#sigma_%s",cat_.c_str()));

  RooCBShape *cb = new RooCBShape(Form("cb_%s",cat_.c_str()),Form("cb_%s",cat_.c_str()), *mass,*mean_cb,*sig_fit_cb, *a1_cb, *n1_cb);

  RooAbsReal *sigma_gaus = fitParams_125["sigma_mh125_gaus"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter sigma_gaus"  << " ? " << sigma_gaus << std::endl;
  sigma_gaus->SetName(Form("sigma_gaus_%s",cat_.c_str()));

  RooAbsReal *dm_gaus = fitParams_125["dm_mh125_gaus"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter dm_gaus"  << " ? " << dm_gaus << std::endl;
  dm_gaus->SetName(Form("dm_%s",cat_.c_str()));

  RooAbsReal *mean_gaus = new RooFormulaVar(Form("mean_gaus_%s",cat_.c_str()),Form("mean_gaus_%s",cat_.c_str()),"@0+@1",RooArgList(*MH,*dm_gaus));

  RooGaussian *gaus = new RooGaussian(Form("gaus_%s",cat_.c_str()),Form("gaus_%s",cat_.c_str()),*mass,*mean_gaus,*sigma_gaus);
  RooAbsReal *frac = fitParams_125["frac_mh125"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter frac_constrained"  << " ? " << frac<< std::endl;
  frac->SetName(Form("frac_%s",cat_.c_str()));

  RooArgList *pdfs_holder = new RooArgList();
  pdfs_holder->add(*cb);
  pdfs_holder->add(*gaus);
  RooArgList *coeffs_holder= new RooArgList();
  coeffs_holder->add(*frac);
  // add the DCB and the Gaussian
  RooAddPdf *pdf = new RooAddPdf(Form("%s",cat_.c_str()),Form("%s",cat_.c_str()),*pdfs_holder,*coeffs_holder,true);

  if (verbosity_){
    std::cout << " [INFO] build a Crystal Ball + 1 Gaussian called " << pdf->GetName() << " with the following parameters" << std::endl;
    pdf->Print();
    a1_cb->Print();
    n1_cb->Print();
    mean_cb->Print();
    sig_fit_cb->Print();
    sigma_gaus->Print();
    mean_gaus->Print();
    frac->Print();
  }
  return pdf;
}

RooAbsPdf* build_DCBpGaus_Pdf(string cat_, parlist_t fitParams, RooRealVar *mass, RooRealVar *MH, bool verbosity_){
 
  string catname;
  catname = Form("%s",cat_.c_str());

  map<string,RooRealVar*> fitParams_125 = fitParams[125];

  RooAbsReal *a1_dcb = fitParams_125["a1_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter a1_dcb"  << " ? " << a1_dcb << std::endl;
  a1_dcb->SetName(Form("a1_%s",cat_.c_str()));

  RooAbsReal *n1_dcb = fitParams_125["n1_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter n1_dcb"  << " ? " << n1_dcb << std::endl;
  n1_dcb->SetName(Form("n1_%s",cat_.c_str()));

  RooAbsReal *a2_dcb = fitParams_125["a2_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter a2_dcb"  << " ? " << a2_dcb << std::endl;
  a2_dcb->SetName(Form("a2_%s",cat_.c_str()));

  RooAbsReal *n2_dcb = fitParams_125["n2_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter n2_dcb"  << " ? " << n2_dcb << std::endl;
  n2_dcb->SetName(Form("n2_%s",cat_.c_str()));

  RooAbsReal *dm_dcb = fitParams_125["dm_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter dm_dcb"  << " ? " << dm_dcb << std::endl;
  dm_dcb->SetName(Form("dm_%s",cat_.c_str()));

  RooAbsReal *mean_dcb = new RooFormulaVar(Form("mean_dcb_%s",cat_.c_str()),Form("mean_dcb_%s",cat_.c_str()),"@0+@1",RooArgList(*MH,*dm_dcb));

  RooAbsReal *sig_fit_dcb = fitParams_125["sigma_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter sigma_dcb"  << " ? " << sig_fit_dcb << std::endl;
  sig_fit_dcb->SetName(Form("#sigma_%s",cat_.c_str()));

  RooAbsPdf *dcb = new RooDoubleCBFast(Form("dcb_%s",cat_.c_str()),Form("dcb_%s",cat_.c_str()), *mass,*mean_dcb,*sig_fit_dcb, *a1_dcb, *n1_dcb, *a2_dcb, *n2_dcb);

  RooAbsReal *sigma_gaus = fitParams_125["sigma_mh125_gaus"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter sigma_gaus"  << " ? " << sigma_gaus << std::endl;
  sigma_gaus->SetName(Form("sigma_gaus_%s",cat_.c_str()));

  RooGaussian *gaus = new RooGaussian(Form("gaus_%s",cat_.c_str()),Form("gaus_%s",cat_.c_str()),*mass,*mean_dcb,*sigma_gaus);
  RooAbsReal *frac = fitParams_125["frac_mh125_gaus"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter frac_constrained"  << " ? " << frac<< std::endl;
  frac->SetName(Form("frac_%s",cat_.c_str()));

  RooArgList *pdfs_holder = new RooArgList();
  pdfs_holder->add(*dcb);
  pdfs_holder->add(*gaus);
  RooArgList *coeffs_holder= new RooArgList();
  coeffs_holder->add(*frac);
  // add the DCB and the Gaussian
  RooAddPdf *pdf = new RooAddPdf(Form("%s",cat_.c_str()),Form("%s",cat_.c_str()),*pdfs_holder,*coeffs_holder,true);

  if (verbosity_){
    std::cout << " [INFO] build a Double Crystal Ball + 1 Gaussian (same mean) called " << pdf->GetName() << " with the following parameters" << std::endl;
    pdf->Print();
    a1_dcb->Print();
    a2_dcb->Print();
    n1_dcb->Print();
    n2_dcb->Print();
    mean_dcb->Print();
    sig_fit_dcb->Print();
    sigma_gaus->Print();
    frac->Print();
  }
  return pdf;
}

RooAbsPdf* build_DCB_Pdf(string cat_, parlist_t fitParams, RooRealVar *mass, RooRealVar *MH, bool verbosity_){
 
  string catname;
  catname = Form("%s",cat_.c_str());

  map<string,RooRealVar*> fitParams_125 = fitParams[125];

  RooAbsReal *a1_dcb = fitParams_125["a1_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter a1_dcb"  << " ? " << a1_dcb << std::endl;
  a1_dcb->SetName(Form("a1_%s",cat_.c_str()));

  RooAbsReal *n1_dcb = fitParams_125["n1_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter n1_dcb"  << " ? " << n1_dcb << std::endl;
  n1_dcb->SetName(Form("n1_%s",cat_.c_str()));

  RooAbsReal *a2_dcb = fitParams_125["a2_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter a2_dcb"  << " ? " << a2_dcb << std::endl;
  a2_dcb->SetName(Form("a2_%s",cat_.c_str()));

  RooAbsReal *n2_dcb = fitParams_125["n2_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter n2_dcb"  << " ? " << n2_dcb << std::endl;
  n2_dcb->SetName(Form("n2_%s",cat_.c_str()));

  RooAbsReal *dm_dcb = fitParams_125["dm_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter dm_dcb"  << " ? " << dm_dcb << std::endl;
  dm_dcb->SetName(Form("dm_%s",cat_.c_str()));

  RooAbsReal *mean_dcb = new RooFormulaVar(Form("mean_dcb_%s",cat_.c_str()),Form("mean_dcb_%s",cat_.c_str()),"@0+@1",RooArgList(*MH,*dm_dcb));

  RooAbsReal *sig_fit_dcb = fitParams_125["sigma_mh125_dcb"];
  if (verbosity_) std::cout << "[INFO] retrieved parameter sigma_dcb"  << " ? " << sig_fit_dcb << std::endl;
  sig_fit_dcb->SetName(Form("#sigma_%s",cat_.c_str()));

  RooAbsPdf *pdf = new RooDoubleCBFast(Form("dcb_%s",cat_.c_str()),Form("dcb_%s",cat_.c_str()), *mass,*mean_dcb,*sig_fit_dcb, *a1_dcb, *n1_dcb, *a2_dcb, *n2_dcb);


  if (verbosity_){
    std::cout << " [INFO] build a Double Crystal Ball (same mean) called " << pdf->GetName() << " with the following parameters" << std::endl;
    pdf->Print();
    a1_dcb->Print();
    a2_dcb->Print();
    n1_dcb->Print();
    n2_dcb->Print();
    mean_dcb->Print();
    sig_fit_dcb->Print();
  }
  
  return pdf;
}

RooAbsPdf* build_SumofGaus_Pdf(string cat_, int nGaussians, parlist_t fitParams, RooRealVar *mass, RooRealVar *MH, bool verbosity_, bool recursive){
  
  string ext = cat_;
  
  map<string,RooRealVar*> fitParams_125 = fitParams[125]; 
 
  RooArgList *gaussians = new RooArgList();
  RooArgList *coeffs = new RooArgList();
 
  for (int g=0; g<nGaussians; g++){
    RooAbsReal *dm = fitParams_125[Form("dm_mh125_g%d",g)];
    if (verbosity_) std::cout << "[INFO] retrieved parameter " << Form("dm_g%d",g) << " ? " << dm << std::endl;
    dm->SetName(Form("dm_g%d_%s",g,ext.c_str()));
 
    RooAbsReal *mean = new RooFormulaVar(Form("mean_mh%d_g%d",mh,g),Form("mean_mh%d_g%d",mh,g),"@0+@1",RooArgList(*MH,*dm));
    mean->SetName(Form("mean_g%d_%s",g,ext.c_str()));

    RooAbsReal *sigma = fitParams_125[Form("sigma_mh125_g%d",g)];
    if (verbosity_) std::cout << "[INFO] retrieved parameter " << Form("sigma_g%d",g) << " ? " << sigma << std::endl;
    sigma->SetName(Form("sigma_g%d_%s",g,ext.c_str()));

    RooGaussian *gaus = new RooGaussian(Form("gaus_g%d_%s",g,ext.c_str()),Form("gaus_g%d_%s",g,ext.c_str()),*mass,*mean,*sigma);
    gaussians->add(*gaus);

    if (g < nGaussians-1) {
      RooAbsReal *frac = fitParams_125[Form("frac_mh125_g%d",g)];
      frac->SetName(Form("frac_g%d_%s", g, ext.c_str()));
      coeffs->add(*frac);
    }
  }
  RooAbsPdf *pdf = new RooAddPdf(Form("%s",ext.c_str()),Form("%s",ext.c_str()),*gaussians,*coeffs,recursive);
  return pdf;
}

int main(int argc, char *argv[]){

  gROOT->SetBatch();
  RooWorkspace *inWS;
  massList_.push_back(125);
  split(flashggCats_,flashggCatsStr_,boost::is_any_of(","));

  TFile *inFile = TFile::Open("/afs/hep.wisc.edu/user/kaho/CMSSW_10_2_16_UL/src/UWHiggs2016/em/signalws.root");
  inWS = (RooWorkspace*)inFile->Get("CMS_emu_workspace");

  TFile *outFile = new TFile("/afs/hep.wisc.edu/user/kaho/CMSSW_10_2_16_UL/src/UWHiggs2016/em/trial.root","RECREATE");
  RooWorkspace *outWS;
  outWS = new RooWorkspace("wsig_13TeV");
  
  for (unsigned int itag =0 ; itag < flashggCats_.size() ; itag++){
    map_cat_.push_back(flashggCats_[itag]);
  }

  for (unsigned int iLine = 0 ; iLine < map_cat_.size() ; iLine++){
    string cat = map_cat_[iLine];
    cout << cat << endl;
    RooDataSet *data0 = (RooDataSet*)inWS->data(Form("Signal_13TeV_%s",cat.c_str()));
    
    cout << Form("Signal_13TeV_%s",cat.c_str()) << endl;
    data0->Print("V");
    map<int,RooDataSet*> FITdatasets;
    FITdatasets.insert(pair<int,RooDataSet*>(125,data0));

    mass_ = (RooRealVar*)inWS->var("CMS_emu_Mass");
    RooRealVar *MH = new RooRealVar("MH","m_{H}",mhLow_,mhHigh_);
      
    InitialFit initFit(mass_,MH,mhLow_,mhHigh_,skipMasses_,true,nBins_,massList_);
    initFit.setVerbosity(true);
//    initFit.buildDCB(cat);
//    initFit.buildDCBplusGaussian(cat);
//    initFit.buildSumOfGaussians(cat, 2, true);
    initFit.buildCBplusGaussian(cat);
    initFit.setDatasets(FITdatasets);
    initFit.runFits(ncpu_);
    parlist_t fitParams = initFit.getFitParams();
    RooAbsPdf* finalpdf = build_CBpGaus_Pdf(cat, fitParams, mass_, MH, true);
//    RooAbsPdf* finalpdf = build_DCB_Pdf(cat, fitParams, mass_, MH, true);
//    RooAbsPdf* finalpdf = build_DCBpGaus_Pdf(cat, fitParams, mass_, MH, true);
//    RooAbsPdf* finalpdf = build_SumofGaus_Pdf(cat, 2, fitParams, mass_, MH, true);
    outWS->import(*data0);
    outWS->import(*finalpdf);
    
    initFit.printFitParams();
  }

  outFile->cd();
  outWS->Write();
  outFile->Close();
}
