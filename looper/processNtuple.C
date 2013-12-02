
#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "dilepStudyLooper.h"

#include <iostream>
#endif

void pickSkimIfExists( TChain *ch, const std::string& base, const std::string& skimPrefix = "" )
{
  TChain *dummy = new TChain("Events");

  int nFiles = 0;
  if (dummy->Add(base.c_str())) {
    nFiles = ch->Add(base.c_str());
    std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
              << nFiles << " files" << std::endl;
  } else
    std::cout << "Couldn't find " << base << std::endl;

  // be paranoid
  if (nFiles == 0) {
    std::cout << "ERROR: expected to read files " 
              << base.c_str() << "  but found none" << std::endl;
    //assert(0);
    exit(0);
  }

  return;
}

void processNtuple( TString outfileid = "me_test", TString infile = "/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/MuEG_Run2012D-PromptReco-v1_AOD/merged/merged_ntuple_6.root", int sign = 0, int em = 0 )
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-00-02";
  const char* jsonfile   = "jsons/Cert_198050-207279_8TeV_19p47ifb_Collisions12_JSON_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  // Load Everything
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libdilepStudyCORE.so");
  gSystem->Load("libdilepStudyLooper.so");

  dilepStudyLooper* looper = new dilepStudyLooper();

  //make baby ntuple
  //  looper->set_createTree(1);
  //set version
  looper->set_version(version);
  //set json
  looper->set_json( jsonfile );
  //batch mode
  //  looper->SetBatchMode(true);

  TChain *chain = new TChain("Events");
  pickSkimIfExists(chain, infile.Data());

  //-------------------------------------
  //set name to get comprehensible output
  //-------------------------------------
  char* sample;
  //MC
  if (infile.Contains("TTJets_MassiveBinDECAY_TuneZ2star_8TeV"))     sample = Form("ttall_%s",  	 outfileid.Data());
  else if (infile.Contains("TTbar_40PU")) sample = Form("ttall_40PU_%s",        outfileid.Data());
  else if (infile.Contains("TTbar_60PU")) sample = Form("ttall_60PU_%s",        outfileid.Data());
  else if (infile.Contains("GluGluToHToGG_M-125")) sample = Form("hgg_m125_%s",  	 outfileid.Data());
  else if (infile.Contains("GluGluToHToWWTo2LAndTau2Nu_M-125")) sample = Form("hww_m125_%s",  	 outfileid.Data());
  else if (infile.Contains("DYJetsToLL"))                            sample = Form("DYtot_%s",           outfileid.Data());
  //Data
  //single mu-had
  else if (infile.Contains("MuHad_Run2012A-recover-06Aug2012-v1"))          sample =  Form("MuHad2012A_recover06Aug2012v1V532_%s",     outfileid.Data());
  else if (infile.Contains("MuHad_Run2012A-13Jul2012-v1"))                  sample =  Form("MuHad2012A_13Jul2012v1V532_%s",            outfileid.Data());
  //single muon
  else if (infile.Contains("SingleMu_Run2012A-recover-06Aug2012-v1"))       sample =  Form("SingleMu2012A_recover06Aug2012v1V532_%s",     outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012A-13Jul2012-v1"))       	sample =  Form("SingleMu2012A_13Jul2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012B-13Jul2012-v1"))       	sample =  Form("SingleMu2012B_13Jul2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012C-24Aug2012-v1"))      		sample =  Form("SingleMu2012C_24Aug2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012C-PromptReco-v2"))      	sample =  Form("SingleMu2012C_PromptRecov2V532_%s",           outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012D-PromptReco-v1"))      	sample =  Form("SingleMu2012D_PromptRecov1V532_%s",           outfileid.Data());
  //single electron
  else if (infile.Contains("SingleElectron_Run2012A-recover-06Aug2012-v1")) sample =  Form("SingleElectron2012A_recover06Aug2012V532_%s", outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012A-13Jul2012-v1"))       	sample =  Form("SingleElectron2012A_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012B-13Jul2012-v1"))       	sample =  Form("SingleElectron2012B_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012C-24Aug2012-v1"))      	sample =  Form("SingleElectron2012C_24Aug2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012C-PromptReco-v2"))      	sample =  Form("SingleElectron2012C_PromptRecov2V532_%s",     outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012D-PromptReco-v1"))      	sample =  Form("SingleElectron2012D_PromptRecov1V532_%s",     outfileid.Data());
  //dimuon 
  else if (infile.Contains("DoubleMu_Run2012A-recover-06Aug2012-v1")) 	sample =  Form("DoubleMu2012A_recover06Aug2012V532_%s", outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012A-13Jul2012-v1"))     		sample =  Form("DoubleMu2012A_13Jul2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012B-13Jul2012-v4"))     		sample =  Form("DoubleMu2012B_13Jul2012v4V532_%s",            outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012C-24Aug2012-v1"))    		sample =  Form("DoubleMu2012C_24Aug2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012C-PromptReco-v2"))    		sample =  Form("DoubleMu2012C_PromptRecov2V532_%s",           outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012D-PromptReco-v1"))    		sample =  Form("DoubleMu2012D_PromptRecov1V532_%s",           outfileid.Data());
  //electron+muon
  else if (infile.Contains("MuEG_Run2012A-recover-06Aug2012-v1"))      	sample =  Form("MuEG2012A_recover06Aug2012V532_%s",           outfileid.Data());
  else if (infile.Contains("MuEG_Run2012A-13Jul2012-v1"))      		sample =  Form("MuEG2012A_13Jul2012v1V532_%s",      	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012B-13Jul2012-v1"))      		sample =  Form("MuEG2012B_13Jul2012v1V532_%s",      	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012C-24Aug2012-v1"))     		sample =  Form("MuEG2012C_24Aug2012v1V532_%s",     	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012C-PromptReco-v2"))     		sample =  Form("MuEG2012C_PromptRecov2V532_%s",     	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012D-PromptReco-v1"))     		sample =  Form("MuEG2012D_PromptRecov1V532_%s",     	      outfileid.Data());
  //dielectron
  else if (infile.Contains("DoubleElectron_Run2012A-recover-06Aug2012-v1")) sample =  Form("DoubleElectron2012A_recover06Aug2012V532_%s", outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012A-13Jul2012-v1"))      	sample =  Form("DoubleElectron2012A_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012B-13Jul2012-v1"))      	sample =  Form("DoubleElectron2012B_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012C-24Aug2012-v1"))     	sample =  Form("DoubleElectron2012C_24Aug2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012C-PromptReco-v2"))     	sample =  Form("DoubleElectron2012C_PromptRecov2V532_%s",     outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012D-PromptReco-v1"))     	sample =  Form("DoubleElectron2012D_PromptRecov1V532_%s",     outfileid.Data());
  //otherwise
  else sample = Form("boiade_%s", outfileid.Data());

  cout<<"sample is "<<sample<<endl;

  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  looper->ScanChain(chain, TString(sample), sign, em);

  //  gSystem->Exit(0);
 
}
