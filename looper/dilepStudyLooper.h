#ifndef dilepStudyLooper_h
#define dilepStudyLooper_h

#include <vector>
#include <list>
#include <string>
#include <map>
#include <set>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

struct isovals {
  float chiso00;
  float chiso04;
  float chiso07;
  float chiso10;
};


class dilepStudyLooper
{
public: 
  dilepStudyLooper();
  ~dilepStudyLooper() {}

  int  ScanChain(TChain *chain, const TString& prefix = "" );
  void BookHistos (const TString& prefix);
  void InitBaby();
  float electronPFiso(const unsigned int index, const bool cor = false);
  float muonPFiso(const unsigned int imu, const bool cor = false);
  float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );
  void labelAxis(TH2F* h, int axis, int lep);

  //
  // WW Electron Id
  //
  bool   ww_elBase(unsigned int i);
  bool   ww_eld0PV(unsigned int i);
  bool   ww_eldZPV(unsigned int i);
  bool ElectronFOV4(unsigned int i);
  bool ElectronFOIdV4(unsigned int i);
  double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);
  int primaryVertex();

  // Set globals
  void set_createTree   (bool  b)    { g_createTree   = b; }
  void set_version      (const char* v)    { g_version      = v; }
  void set_json         (const char* v)    { g_json         = v; }        

  // Baby ntuple methods
  void makeOutput (const TString& prefix);
  void closeOutput ();

  // trigger-object matching
  int findTriggerIndex(const TString& trigName);
  TString triggerName(const TString& triggerPattern);
  bool objectPassTrigger(const LorentzVector &obj, const TString& trigname, int type, float drmax = 0.1 );

  // utils
  float getdphi( float phi1 , float phi2 );
  isovals muonChIsoValuePF2012 (const unsigned int imu, const float R = 0.3, const int ivtx = 0);

private:

  // Globals
  bool  g_createTree;
  const char* g_version;
  const char* g_json;      
  bool initialized;
  bool isdata_;
  TFile* outFile;

  float ptthresh_high;
  float ptthresh_low;

  // histograms
  TH1F* h_nvtx;

  TH1F* h_el_pt;
  TH1F* h_el_lead_pt;
  TH1F* h_el_subl_pt;
  TH1F* h_el_lead_eta;
  TH1F* h_el_subl_eta;
  TH1F* h_ee_mll;
  TH1F* h_ee_denom_mll;
  TH1F* h_ee_mediso_mll;

  TH1F* h_mu_pt;
  TH1F* h_mu_lead_pt;
  TH1F* h_mu_subl_pt;
  TH1F* h_mu_lead_eta;
  TH1F* h_mu_subl_eta;
  TH1F* h_mm_mll;
  TH1F* h_mu_dr;
  TH2F* h_mm_mll_vs_dr;
  TH1F* h_mm_chargeprod;

  TH1F* h_mu_tiso_lead_pt;
  TH1F* h_mu_tiso_subl_pt;
  TH1F* h_mu_tiso_lead_eta;
  TH1F* h_mu_tiso_subl_eta;
  TH1F* h_mm_tiso_mll;

  TH1F* h_mu_dup_dpt;
  TH1F* h_mu_dup_deta;
  TH1F* h_mu_dup_dphi;
  TH1F* h_mu_dup_dr;
  TH1F* h_mu_dup_lead_type;
  TH1F* h_mu_dup_subl_type;

  TH1F* h_mu_pfchiso;
  TH1F* h_mu_pfchiso_minus_chiso;
  TH1F* h_mu_chiso04;
  TH1F* h_mu_chiso07;
  TH1F* h_mu_chiso10;

  TH1F* h_em_el_pt;
  TH1F* h_em_mu_pt;
  TH1F* h_em_el_eta;
  TH1F* h_em_mu_eta;
  TH1F* h_em_mll;
  TH1F* h_em_dr;

  TH1F* h_me_mu_pt;
  TH1F* h_me_el_pt;
  TH1F* h_me_mu_eta;
  TH1F* h_me_el_eta;
  TH1F* h_me_mll;
  TH1F* h_me_dr;

  TH1F* h_em_el_tiso_pt;
  TH1F* h_em_mu_tiso_pt;
  TH1F* h_em_el_tiso_eta;
  TH1F* h_em_mu_tiso_eta;
  TH1F* h_em_tiso_mll;
  TH1F* h_em_tiso_dr;

  TH1F* h_me_mu_tiso_pt;
  TH1F* h_me_el_tiso_pt;
  TH1F* h_me_mu_tiso_eta;
  TH1F* h_me_el_tiso_eta;
  TH1F* h_me_tiso_mll;
  TH1F* h_me_tiso_dr;

  TH1F* h_em_el_hlt_pt;
  TH1F* h_em_mu_hlt_pt;
  TH1F* h_em_el_hlt_eta;
  TH1F* h_em_mu_hlt_eta;

  TH1F* h_me_mu_hlt_pt;
  TH1F* h_me_el_hlt_pt;
  TH1F* h_me_mu_hlt_eta;
  TH1F* h_me_el_hlt_eta;

  TH1F* h_em_el_hlt_noreco_pt;
  TH1F* h_em_mu_hlt_noreco_pt;
  TH1F* h_em_el_hlt_noreco_eta;
  TH1F* h_em_mu_hlt_noreco_eta;

  TH1F* h_me_mu_hlt_noreco_pt;
  TH1F* h_me_el_hlt_noreco_pt;
  TH1F* h_me_mu_hlt_noreco_eta;
  TH1F* h_me_el_hlt_noreco_eta;

  TH1F* h_mm_overlap;

  TH2F* h_ee_events;
  TH2F* h_em_events;
  TH2F* h_me_events;
  TH2F* h_mm_events;
  TH2F* h_mmtk_events;
  TH2F* h_mmor_events;

};

#endif
