#include "dilepStudyLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
//#include "../Tools/pfjetMVAtools.h"

#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/jetSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetSmearingTools.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"

bool verbose              = false;
bool doTenPercent         = false;
bool doLowerPtThresh      = false;
bool doEM                 = true;
bool doME                 = false;
bool requireTrigMatch     = true;
bool doSS                 = false;
bool doOS                 = true;

using namespace std;
using namespace tas;

//--------------------------------------------------------------------

int getMotherIndex(int motherid){
  for(unsigned int i = 0; i < genps_id().size() ; i++){
    if( motherid == genps_id().at(i) ) return i;
  }

  return -1;
}

//--------------------------------------------------------------------

float dilepStudyLooper::dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 ) { 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}

//--------------------------------------------------------------------

dilepStudyLooper::dilepStudyLooper()
{

  std::cout << " construct " << std::endl;
  g_createTree   = false;
  initialized = false;
}

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

// void dilepStudyLooper::InitBaby(){
// }

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void dilepStudyLooper::closeOutput()
{
  outFile->cd();
  //  outTree->Write();
  outFile->Write();
  outFile->Close();
  delete outFile;
}

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

int dilepStudyLooper::ScanChain(TChain* chain, const TString& prefix)

{

  //  cout << "ciao " << isData << endl;

  bool isData = false;
  bool isEE = false;
  bool isEM = false;
  bool isMM = false;
  bool isDYMC = false;
  if( prefix.Contains("data") || prefix.Contains("2012") 
      || prefix.Contains("dimu") || prefix.Contains("diel")
      || prefix.Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
    if (prefix.Contains("DoubleElectron")) {
      isEE = true;
      std::cout << "DoubleElectron data" << std::endl;
    }
    else if (prefix.Contains("DoubleMu")) {
      isMM = true;
      std::cout << "DoubleMu data" << std::endl;
    }
    else if (prefix.Contains("MuEG")) {
      isEM = true;
      std::cout << "MuEG data" << std::endl;
    }
  }

  cout << "IS DATA: " << isData << endl;

  if ( prefix.Contains("DYtot") ) {
    std::cout << "Running on DY MC" << std::endl;
    isDYMC = true;
    //    isEE = true;
    isMM = true;
  }

  if( doTenPercent ) cout << "Processing 10% of MC" << endl;

  ptthresh_high = 17.;
  ptthresh_low = 8.;
  if (doLowerPtThresh) {
    ptthresh_high = 15.;
    ptthresh_low = 5.;
  }

  const TString trigname_mm = "HLT_Mu17_Mu8_v";
  const TString trigname_mmtk = "HLT_Mu17_TkMu8_v";
  const TString trigname_me = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
  const TString trigname_em = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
  const TString trigname_ee = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";

  //------------------------------------------------------------------------------------------------------
  // set json\, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){

    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    //    if( prefix.Contains("ttall_massivebin") ) 
    //    set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS10_19fb_Zselection.root",true);

    initialized = true;
  }

  // -----------------------------------------------------------

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  makeOutput(prefix);
  //  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  //  BookHistos(prefix);
  BookHistos("h");

  cout << " done with initialization "  << endl;
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nEventsPreReco = 0;
  unsigned int nEventsPass = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  // test chain
  if (!chain)
    {
      throw std::invalid_argument("at::ScanChain: chain is NULL!");
    }
  if (chain->GetListOfFiles()->GetEntries()<1)
    {
      throw std::invalid_argument("at::ScanChain: chain has no files!");
    }
  if (not chain->GetFile())
    {
      throw std::invalid_argument("at::ScanChain: chain has no files or file path is invalid!");
    }
  int nSkip_els_conv_dist = 0;

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    cout << currentFile->GetTitle() << endl;

    if (!f || f->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain is invalid or corrupt: %s", currentFile->GetTitle()));
    }
    
    // get the trees in each file
    // TTree *tree = (TTree*)f->Get("Events");
    TTree *tree = dynamic_cast<TTree*>(f->Get("Events"));
    if (!tree || tree->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain has an invalid TTree or is corrupt: %s", currentFile->GetTitle()));
    }

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();
    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      /////////      cout << nEventsTotal << endl;

      if( doTenPercent ){
	if( !(nEventsTotal%10==0) ) continue;
      }

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //Matevz
      tree->LoadTree(z);

      cms2.GetEntry(z);

      if( evt_ww_rho_vor() != evt_ww_rho_vor() ){
	cout << "Skipping event with rho = nan!!!" << endl;
	continue;
      }

      //      InitBaby();

      isdata_ = isData ? 1 : 0;

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      TString datasetname(evt_dataset().at(0));

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      //      if( !cleaning_goodVertexApril2011() )                          continue;
      if( isData && !goodrun(evt_run(), evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------

      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
          skipEvent = true;
        }
      }
             
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }

      //---------------------------------------------
      // count vertices
      //---------------------------------------------
      int nvtx = 0;
      for (size_t v = 0; v < vtxs_position().size(); ++v){
	if(isGoodVertex(v)) ++nvtx;
      }

      //---------------------------------------------
      // gen selection for DY MC
      //---------------------------------------------

      if (isDYMC) {
	int nel = 0;
	int nmu = 0;
	// require 2 gen muons, to study Z->mm
	for (unsigned int igen = 0; igen < genps_p4().size(); ++igen) {
	  if (genps_status().at(igen) != 3) continue;
	  int id = genps_id().at(igen);
	  if ((abs(id) != 11) && (abs(id) != 13)) continue;
	  if (abs(genps_id_mother().at(igen)) != 23) continue;
	  if (abs(id) == 11) ++nel;
	  else if (abs(id) == 13) ++nmu;
	}
	if (nmu < 2) continue;
      }

      //---------------------------------------------
      // trigger selection
      //---------------------------------------------
      bool pass_trig_mm        = passUnprescaledHLTTriggerPattern(trigname_mm.Data())                                     ? 1 : 0;
      bool pass_trig_mmtk      = passUnprescaledHLTTriggerPattern(trigname_mmtk.Data())                                   ? 1 : 0;
      bool pass_trig_me        = passUnprescaledHLTTriggerPattern(trigname_me.Data()) ? 1 : 0;
      bool pass_trig_em        = passUnprescaledHLTTriggerPattern(trigname_em.Data()) ? 1 : 0;
      bool pass_trig_ee        = passUnprescaledHLTTriggerPattern(trigname_ee.Data()) ? 1 : 0;

      // require the event to pass one of these triggers (only for data)
      if (isData) {
	if (isEE && !pass_trig_ee) continue;
	if (isEM) {
	  bool pass = false;
	  if (doEM && pass_trig_em) pass = true;
	  if (doME && pass_trig_me) pass = true;
	  if (!pass) continue;
	}
	if (isMM && !pass_trig_mm && !pass_trig_mmtk) continue;
	if (!pass_trig_mm && !pass_trig_mmtk && !pass_trig_me && !pass_trig_em && !pass_trig_ee) continue;
      }

      //---------------------------------------------
      // reco electron selection
      //---------------------------------------------

      int minlep = 1;
      if (isEE || isMM) minlep = 2;

      std::vector<int> el_lead_flags;
      int nel20_cand = 0;
      int nel20_wwdenom = 0;
      int nel20_loose = 0;
      int nel20_med = 0;
      int nel20_iso = 0;
      int nel20_loose_iso = 0;
      int nel20_med_iso = 0;

      std::vector<int> el_subl_flags;
      int nel10_cand = 0;
      int nel10_wwdenom = 0;
      int nel10_loose = 0;
      int nel10_med = 0;
      int nel10_iso = 0;
      int nel10_loose_iso = 0;
      int nel10_med_iso = 0;

      float pt_lead_el = -1.;
      int idx_lead_el = -1;
      float pt_subl_el = -1.;
      int idx_subl_el = -1;

      float pt_lead_denom_el = -1.;
      int idx_lead_denom_el = -1;
      float pt_subl_denom_el = -1.;
      int idx_subl_denom_el = -1;

      float pt_lead_mediso_el = -1.;
      int idx_lead_mediso_el = -1;
      float pt_subl_mediso_el = -1.;
      int idx_subl_mediso_el = -1;

      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) {
	// require match to trigger object
	if (requireTrigMatch) {
	  bool matched = false;
	  if (isEE) {
	    matched = objectPassTrigger(els_p4().at(iel),trigname_ee,82);
	  } else if (isEM) {
	    if (doEM) matched |= objectPassTrigger(els_p4().at(iel),trigname_em,82);
	    if (doME) matched |= objectPassTrigger(els_p4().at(iel),trigname_me,82);
	  }
	  if (!matched) continue;
	}

	// cut on pt and eta
	//	if (fabs(els_p4().at(iel).eta()) > 2.4) continue;
	float pt = els_p4().at(iel).pt();
	h_el_pt->Fill(pt);
	if (pt > pt_lead_el) {
	  pt_subl_el = pt_lead_el;
	  idx_subl_el = idx_lead_el;
	  pt_lead_el = pt;
	  idx_lead_el = iel;
	} else if (pt > pt_subl_el) {
	  pt_subl_el = pt;
	  idx_subl_el = iel;
	}
	if (pt < ptthresh_low) continue;
	// WW denom def
        bool pass_wwdenom = ElectronFOV4(iel);
	// loose ID, no iso
        bool pass_loose = passElectronSelection_ZMet2012_v3_NoIso( iel,true,true,false);
	// med ID, no iso
        bool pass_med = passElectronSelection_Stop2012_v3_NoIso( iel,true,true,false);
	// iso cut: pfiso/pt < 0.15
	float iso_cor = electronPFiso(iel,true);
	bool pass_iso = bool(iso_cor < 0.15);

	if (pass_wwdenom) {
	  if (pt > pt_lead_denom_el) {
	    pt_subl_denom_el = pt_lead_denom_el;
	    idx_subl_denom_el = idx_lead_denom_el;
	    pt_lead_denom_el = pt;
	    idx_lead_denom_el = iel;
	  } else if (pt > pt_subl_denom_el) {
	    pt_subl_denom_el = pt;
	    idx_subl_denom_el = iel;
	  }
	}

	if (pass_med && pass_iso) {
	  if (pt > pt_lead_mediso_el) {
	    pt_subl_mediso_el = pt_lead_mediso_el;
	    idx_subl_mediso_el = idx_lead_mediso_el;
	    pt_lead_mediso_el = pt;
	    idx_lead_mediso_el = iel;
	  } else if (pt > pt_subl_mediso_el) {
	    pt_subl_mediso_el = pt;
	    idx_subl_mediso_el = iel;
	  }
	}

	// count after various selections
	// pt20
	if (pt > ptthresh_high) {
	  ++nel20_cand;
	  if (pass_wwdenom) ++nel20_wwdenom;
	  if (pass_loose) ++nel20_loose;
	  if (pass_med) ++nel20_med;
	  if (pass_iso) ++nel20_iso;
	  if (pass_loose && pass_iso) ++nel20_loose_iso;
	  if (pass_med && pass_iso) ++nel20_med_iso;
	}

	// pt 10-20
	else {
	  ++nel10_cand;
	  if (pass_wwdenom) ++nel10_wwdenom;
	  if (pass_loose) ++nel10_loose;
	  if (pass_med) ++nel10_med;
	  if (pass_iso) ++nel10_iso;
	  if (pass_loose && pass_iso) ++nel10_loose_iso;
	  if (pass_med && pass_iso) ++nel10_med_iso;
	}
      }

      // set flags based on observed candidates
      el_lead_flags.push_back(0);
      if (nel20_cand > 0) el_lead_flags.push_back(1);
      if (nel20_wwdenom > 0) el_lead_flags.push_back(2);
      if (nel20_loose > 0) el_lead_flags.push_back(3);
      if (nel20_med > 0) el_lead_flags.push_back(4);
      if (nel20_iso > 0) el_lead_flags.push_back(5);
      if (nel20_loose_iso > 0) el_lead_flags.push_back(6);
      if (nel20_med_iso > 0) el_lead_flags.push_back(7);

      el_subl_flags.push_back(0);
      if ((nel20_cand >= minlep) || (nel10_cand > 0))  el_subl_flags.push_back(1);
      if ((nel20_wwdenom >= minlep) || (nel10_wwdenom > 0)) el_subl_flags.push_back(2);
      if ((nel20_loose >= minlep) || (nel10_loose > 0)) el_subl_flags.push_back(3);
      if ((nel20_med >= minlep) || (nel10_med > 0)) el_subl_flags.push_back(4);
      if ((nel20_iso >= minlep) || (nel10_iso > 0)) el_subl_flags.push_back(5);
      if ((nel20_loose_iso >= minlep) || (nel10_loose_iso > 0)) el_subl_flags.push_back(6);
      if ((nel20_med_iso >= minlep) || (nel10_med_iso > 0)) el_subl_flags.push_back(7);

      // check for SS if that's what we're running
      if ((doSS || doOS) && isEE) {
	if (pt_lead_el <= 0. || pt_subl_el <= 0.) continue;
	bool isOS = bool(els_charge().at(idx_lead_el) != els_charge().at(idx_subl_el));
	  if (doSS && isOS) continue;
	  else if (doOS && !isOS) continue;
      }

      // pt hists
      if (pt_lead_el > 0.) {
	h_el_lead_pt->Fill(pt_lead_el);
	h_el_lead_eta->Fill(els_p4().at(idx_lead_el).eta());
	if (pass_trig_em) {
	  h_em_el_pt->Fill(pt_lead_el);
	  h_em_el_eta->Fill(els_p4().at(idx_lead_el).eta());
	}
	if (pass_trig_me) {
	  h_me_el_pt->Fill(pt_lead_el);
	  h_me_el_eta->Fill(els_p4().at(idx_lead_el).eta());
	}
      }
      if (pt_subl_el > 0.) {
	h_el_subl_pt->Fill(pt_subl_el);
	h_el_subl_eta->Fill(els_p4().at(idx_subl_el).eta());
	LorentzVector dilep = els_p4().at(idx_lead_el) + els_p4().at(idx_subl_el);
	h_ee_mll->Fill(dilep.M());
      }

      if (pt_lead_denom_el > 0. && pt_subl_denom_el > 0.) {
	LorentzVector dilep = els_p4().at(idx_lead_denom_el) + els_p4().at(idx_subl_denom_el);
	h_ee_denom_mll->Fill(dilep.M());
      }

      if (pt_lead_mediso_el > 0. && pt_subl_mediso_el > 0.) {
	LorentzVector dilep = els_p4().at(idx_lead_mediso_el) + els_p4().at(idx_subl_mediso_el);
	h_ee_mediso_mll->Fill(dilep.M());
      }

      //---------------------------------------------
      // reco muon selection
      //---------------------------------------------

      std::vector<int> mu_lead_flags;
      int nmu20_cand = 0;
      int nmu20_tight = 0;
      int nmu20_iso = 0;
      int nmu20_iso04 = 0;
      int nmu20_iso07 = 0;
      int nmu20_iso10 = 0;
      int nmu20_tight_iso = 0;

      std::vector<int> mu_subl_flags;
      int nmu10_cand = 0;
      int nmu10_tight = 0;
      int nmu10_iso = 0;
      int nmu10_iso04 = 0;
      int nmu10_iso07 = 0;
      int nmu10_iso10 = 0;
      int nmu10_tight_iso = 0;

      float pt_lead_mu = -1.;
      int idx_lead_mu = -1;
      float pt_subl_mu = -1.;
      int idx_subl_mu = -1;

      float pt_lead_tiso_mu = -1.;
      int idx_lead_tiso_mu = -1;
      float pt_subl_tiso_mu = -1.;
      int idx_subl_tiso_mu = -1;

      // loop once on muons to remove duplicates
      std::vector<int> good_mu_idx;
      std::vector<int> dup_mu_idx;
      for (unsigned int imu = 0; imu < mus_p4().size(); ++imu) {
	// first check duplicate list
	bool duplicate = false;
	for (unsigned int jmu = 0; jmu < dup_mu_idx.size(); ++jmu) {
	  if (dup_mu_idx.at(jmu) == imu) {
	    duplicate = true;
	    break;
	  }
	}
	if (duplicate) continue;

	// require match to trigger object (if requested)
	if (requireTrigMatch) {
	  bool matched = false;
	  if (isMM) {
	    matched = objectPassTrigger(mus_p4().at(imu),trigname_mm,83) || objectPassTrigger(mus_p4().at(imu),trigname_mmtk,83);
	  } else if (isEM) {
	    if (doEM) matched |= objectPassTrigger(mus_p4().at(imu),trigname_em,83);
	    if (doME) matched |= objectPassTrigger(mus_p4().at(imu),trigname_me,83);
	  }
	  if (!matched) continue;
	}

	// check ID, no iso -- used in duplicate rejection
	bool pass_id_i = muonIdNotIsolated(imu, ZMet2012_v1);

	// // check for duplicate muons using dR 0.1
	duplicate = false;
	for (unsigned int jmu = imu+1; jmu < mus_p4().size(); ++jmu) {
	  if (dRbetweenVectors(mus_p4().at(imu),mus_p4().at(jmu)) < 0.1) {
	    // check whether each passes tight ID
	    bool pass_id_j = muonIdNotIsolated(jmu, ZMet2012_v1);
	    // imu passes, jmu doesn't: keep imu
	    if (pass_id_i && !pass_id_j) {
	      dup_mu_idx.push_back(jmu);
	      continue;
	    } 
	    // jmu passes, imu doesn't: keep jmu
	    else if (!pass_id_i && pass_id_j) {
	      dup_mu_idx.push_back(imu);
	      duplicate = true;
	      break;
	    }
	    // if both or neither pass tight, take the highest pt muon
	    else {
	      if (mus_p4().at(imu).pt() > mus_p4().at(jmu).pt()) {
		dup_mu_idx.push_back(jmu);
		continue;
	      } else {
		dup_mu_idx.push_back(imu);
		duplicate = true;
		break;
	      }
	    }
	  } // if dR < 0.1 
	} // loop over jmu
	if (duplicate) continue;

	good_mu_idx.push_back(imu);
      } // first loop over muons

      // then loop again to save info and make hists
      for (unsigned int jmu = 0; jmu < good_mu_idx.size(); ++jmu) {
	int imu = good_mu_idx.at(jmu);

	// cut on pt and eta
	//	if (fabs(mus_p4().at(imu).eta()) > 2.4) continue;
	float pt = mus_p4().at(imu).pt();

	h_mu_pt->Fill(pt);
	if (pt > pt_lead_mu) {
	  pt_subl_mu = pt_lead_mu;
	  idx_subl_mu = idx_lead_mu;
	  pt_lead_mu = pt;
	  idx_lead_mu = imu;
	} else if (pt > pt_subl_mu) {
	  pt_subl_mu = pt;
	  idx_subl_mu = imu;
	}

	// require ID, no iso
	bool pass_id = muonIdNotIsolated(imu, ZMet2012_v1);
	// iso cut, pfiso/pt < 0.15
	float iso_cor = muonPFiso(imu,true);
	bool pass_iso = bool(iso_cor < 0.15);

	// different iso variables
	float pfchiso = cms2.mus_isoR03_pf_ChargedHadronPt().at(imu);
	isovals chisovals = muonChIsoValuePF2012(imu);

	h_mu_pfchiso->Fill(pfchiso/pt);
	h_mu_pfchiso_minus_chiso->Fill(pfchiso - chisovals.chiso00);
	h_mu_chiso04->Fill(chisovals.chiso04/pt);
	h_mu_chiso07->Fill(chisovals.chiso07/pt);
	h_mu_chiso10->Fill(chisovals.chiso10/pt);

	bool pass_iso04 = (chisovals.chiso07/pt < 0.4);
	bool pass_iso07 = (chisovals.chiso07/pt < 0.7);
	bool pass_iso10 = (chisovals.chiso07/pt < 1.0);

	if (pt < ptthresh_low) continue;

	//	if (pass_iso04) {
	if (pass_id && pass_iso) {
	  if (pt > pt_lead_tiso_mu) {
	    pt_subl_tiso_mu = pt_lead_tiso_mu;
	    idx_subl_tiso_mu = idx_lead_tiso_mu;
	    pt_lead_tiso_mu = pt;
	    idx_lead_tiso_mu = imu;
	  } else if (pt > pt_subl_tiso_mu) {
	    pt_subl_tiso_mu = pt;
	    idx_subl_tiso_mu = imu;
	  }
	}

	// check selections
	// pt20
	if (pt > ptthresh_high) {
	  ++nmu20_cand;
	  if (pass_id)  ++nmu20_tight;
	  if (pass_iso10)  ++nmu20_iso10;
	  if (pass_iso07)  ++nmu20_iso07;
	  if (pass_iso04)  ++nmu20_iso04;
	  if (pass_iso)  ++nmu20_iso;
	  if (pass_id && pass_iso)  ++nmu20_tight_iso;
	}

	// pt 10-20
	else {
	  ++nmu10_cand;
	  if (pass_id) ++nmu10_tight;
	  if (pass_iso10)  ++nmu10_iso10;
	  if (pass_iso07)  ++nmu10_iso07;
	  if (pass_iso04)  ++nmu10_iso04;
	  if (pass_iso) ++nmu10_iso;
	  if (pass_id && pass_iso) ++nmu10_tight_iso;
	}
      }

      // set flags based on observed candidates
      mu_lead_flags.push_back(0);
      if (nmu20_cand > 0) mu_lead_flags.push_back(1);
      if (nmu20_tight > 0) mu_lead_flags.push_back(2);
      if (nmu20_iso10 > 0) mu_lead_flags.push_back(3);
      if (nmu20_iso07 > 0) mu_lead_flags.push_back(4);
      if (nmu20_iso04 > 0) mu_lead_flags.push_back(5);
      if (nmu20_iso > 0) mu_lead_flags.push_back(6);
      if (nmu20_tight_iso > 0) mu_lead_flags.push_back(7);

      mu_subl_flags.push_back(0);
      if ((nmu20_cand >= minlep) || (nmu10_cand > 0)) mu_subl_flags.push_back(1);
      if ((nmu20_tight >= minlep) || (nmu10_tight > 0)) mu_subl_flags.push_back(2);
      if ((nmu20_iso10 >= minlep) || (nmu10_iso10 > 0)) mu_subl_flags.push_back(3);
      if ((nmu20_iso07 >= minlep) || (nmu10_iso07 > 0)) mu_subl_flags.push_back(4);
      if ((nmu20_iso04 >= minlep) || (nmu10_iso04 > 0)) mu_subl_flags.push_back(5);
      if ((nmu20_iso >= minlep) || (nmu10_iso > 0)) mu_subl_flags.push_back(6);
      if ((nmu20_tight_iso >= minlep) || (nmu10_tight_iso > 0)) mu_subl_flags.push_back(7);

      // check for SS if that's what we're running
      if (doSS || doOS) {
	if (isMM) {
	  if (pt_lead_mu <= 0. || pt_subl_mu <= 0.) continue;
	  bool isOS = bool(mus_charge().at(idx_lead_mu) != mus_charge().at(idx_subl_mu));
	  if (doSS && isOS) continue;
	  else if (doOS && !isOS) continue;
	}
	else if (isEM) {
	  if (pt_lead_el <= 0. || pt_lead_mu <= 0.) continue;
	  bool isOS = bool(els_charge().at(idx_lead_el) != mus_charge().at(idx_lead_mu));
	  if (doSS && isOS) continue;
	  else if (doOS && !isOS) continue;
	}
      }

      // pt hists
      if (pt_lead_mu > 0.) {
	h_mu_lead_pt->Fill(pt_lead_mu);
	h_mu_lead_eta->Fill(mus_p4().at(idx_lead_mu).eta());
	if (pass_trig_em) {
	  h_em_mu_pt->Fill(pt_lead_mu);
	  h_em_mu_eta->Fill(mus_p4().at(idx_lead_mu).eta());
	  if (pt_lead_el > 0.) {
	    LorentzVector dilep =  els_p4().at(idx_lead_el) + mus_p4().at(idx_lead_mu);
	    h_em_mll->Fill(dilep.M());
	    h_em_dr->Fill(dRbetweenVectors(els_p4().at(idx_lead_el),mus_p4().at(idx_lead_mu)));
	  }
	  if ((nmu10_tight_iso > 0 || nmu20_tight_iso > 0) && (nel20_cand > 0)) {
	    h_em_el_tiso_pt->Fill(pt_lead_el);
	    h_em_el_tiso_eta->Fill(els_p4().at(idx_lead_el).eta());
	    h_em_mu_tiso_pt->Fill(pt_lead_mu);
	    h_em_mu_tiso_eta->Fill(mus_p4().at(idx_lead_mu).eta());
	    LorentzVector dilep =  els_p4().at(idx_lead_el) + mus_p4().at(idx_lead_mu);
	    h_em_tiso_mll->Fill(dilep.M());
	    h_em_tiso_dr->Fill(dRbetweenVectors(els_p4().at(idx_lead_el),mus_p4().at(idx_lead_mu)));
	  }
	}
	if (pass_trig_me) {
	  h_me_mu_pt->Fill(pt_lead_mu);
	  h_me_mu_eta->Fill(mus_p4().at(idx_lead_mu).eta());
	  if (pt_lead_el > 0.) {
	    LorentzVector dilep =  els_p4().at(idx_lead_el) + mus_p4().at(idx_lead_mu);
	    h_me_mll->Fill(dilep.M());
	    h_me_dr->Fill(dRbetweenVectors(mus_p4().at(idx_lead_mu),els_p4().at(idx_lead_el)));
	  }
	  if ((nmu20_tight_iso > 0) && (nel10_cand > 0 || nel20_cand > 0)) {
	    h_me_el_tiso_pt->Fill(pt_lead_el);
	    h_me_el_tiso_eta->Fill(els_p4().at(idx_lead_el).eta());
	    h_me_mu_tiso_pt->Fill(pt_lead_mu);
	    h_me_mu_tiso_eta->Fill(mus_p4().at(idx_lead_mu).eta());
	    LorentzVector dilep =  els_p4().at(idx_lead_el) + mus_p4().at(idx_lead_mu);
	    h_me_tiso_mll->Fill(dilep.M());
	    h_me_tiso_dr->Fill(dRbetweenVectors(mus_p4().at(idx_lead_mu),els_p4().at(idx_lead_el)));
	  }
	}
      }
      if (pt_subl_mu > 0.) {
	h_mu_subl_pt->Fill(pt_subl_mu);
	h_mu_subl_eta->Fill(mus_p4().at(idx_subl_mu).eta());
	LorentzVector dilep = mus_p4().at(idx_lead_mu) + mus_p4().at(idx_subl_mu);
	h_mm_mll->Fill(dilep.M());
	float dr = dRbetweenVectors(mus_p4().at(idx_lead_mu),mus_p4().at(idx_subl_mu));
	h_mu_dr->Fill(dr);
	h_mm_mll_vs_dr->Fill(dr,dilep.M());
	h_mm_chargeprod->Fill(mus_charge().at(idx_lead_mu) * mus_charge().at(idx_subl_mu));

	//	if (dilep.M() < 1.0) {
	if (dr < 0.1) {
	  h_mu_dup_dpt->Fill(mus_p4().at(idx_lead_mu).pt() - mus_p4().at(idx_subl_mu).pt());
	  h_mu_dup_deta->Fill(mus_p4().at(idx_lead_mu).eta() - mus_p4().at(idx_subl_mu).eta());
	  h_mu_dup_dphi->Fill(getdphi(mus_p4().at(idx_lead_mu).phi(),mus_p4().at(idx_subl_mu).phi()));
	  h_mu_dup_dr->Fill(dr);
	  // muon type, only look at bits for global, tracker, standalone
	  h_mu_dup_lead_type->Fill(mus_type().at(idx_lead_mu) & 0xe);
	  h_mu_dup_subl_type->Fill(mus_type().at(idx_subl_mu) & 0xe);
	  // std::cout << "-- WARNING: dimuon mass of " << dilep.M() << std::endl
	  // 	    << "--- muon 1: pt: " << mus_p4().at(idx_lead_mu).pt() << ", eta: " << mus_p4().at(idx_lead_mu).eta()
	  // 	    << ", phi: " << mus_p4().at(idx_lead_mu).phi();
	  // if ((mus_type().at(idx_lead_mu) & (1<<1)) == 1<<1) std::cout << ", global";
	  // if ((mus_type().at(idx_lead_mu) & (1<<2)) == 1<<2) std::cout << ", tracker";
	  // if ((mus_type().at(idx_lead_mu) & (1<<3)) == 1<<3) std::cout << ", standalone";
	  // std::cout << std::endl
	  // 	    << "--- muon 2: pt: " << mus_p4().at(idx_subl_mu).pt() << ", eta: " << mus_p4().at(idx_subl_mu).eta()
	  // 	    << ", phi: " << mus_p4().at(idx_subl_mu).phi();
	  // if ((mus_type().at(idx_subl_mu) & (1<<1)) == 1<<1) std::cout << ", global";
	  // if ((mus_type().at(idx_subl_mu) & (1<<2)) == 1<<2) std::cout << ", tracker";
	  // if ((mus_type().at(idx_subl_mu) & (1<<3)) == 1<<3) std::cout << ", standalone";
	  // std::cout << std::endl
	  // 	    << "--- dr: " << dr << std::endl;
	}

      }

      if (pt_lead_tiso_mu > 0. && pt_subl_tiso_mu > 0.) {
	h_mu_tiso_lead_pt->Fill(pt_lead_tiso_mu);
	h_mu_tiso_subl_pt->Fill(pt_subl_tiso_mu);
	h_mu_tiso_lead_eta->Fill(mus_p4().at(idx_lead_tiso_mu).eta());
	h_mu_tiso_subl_eta->Fill(mus_p4().at(idx_subl_tiso_mu).eta());
	LorentzVector dilep = mus_p4().at(idx_lead_tiso_mu) + mus_p4().at(idx_subl_tiso_mu);
	h_mm_tiso_mll->Fill(dilep.M());
      }

      //---------------------------------------------
      // putting together event selections
      //---------------------------------------------

      h_nvtx->Fill(nvtx);

      if (isEE && pass_trig_ee) {
	for (unsigned int ilead=0; ilead < el_lead_flags.size(); ++ilead) {
	  for (unsigned int isubl=0; isubl < el_subl_flags.size(); ++isubl) {
	    h_ee_events->Fill(el_lead_flags.at(ilead),el_subl_flags.at(isubl));
	  }
	}
      }

      if (isEM && pass_trig_em) {
	for (unsigned int ilead=0; ilead < el_lead_flags.size(); ++ilead) {
	  for (unsigned int isubl=0; isubl < mu_subl_flags.size(); ++isubl) {
	    h_em_events->Fill(el_lead_flags.at(ilead),mu_subl_flags.at(isubl));
	  }
	}

	// plot hlt obj kinematics
	int trigidx = findTriggerIndex(triggerName(trigname_em));
	std::vector<LorentzVector> trigp4 = hlt_trigObjs_p4()[trigidx];
	std::vector<int> trigid = hlt_trigObjs_id()[trigidx];

	for (unsigned int ihlt = 0; ihlt < trigp4.size(); ++ihlt){
	  if (trigid.at(ihlt) == 82) { // HLT code for electron
	    h_em_el_hlt_pt->Fill(trigp4.at(ihlt).pt());
	    h_em_el_hlt_eta->Fill(trigp4.at(ihlt).eta());
	    if (pt_lead_el < 0.) {
	      h_em_el_hlt_noreco_pt->Fill(trigp4.at(ihlt).pt());
	      h_em_el_hlt_noreco_eta->Fill(trigp4.at(ihlt).eta());
	    }
	  } else if (trigid.at(ihlt) == 83) { // HLT code for muon
	    h_em_mu_hlt_pt->Fill(trigp4.at(ihlt).pt());
	    h_em_mu_hlt_eta->Fill(trigp4.at(ihlt).eta());
	    if (pt_lead_mu < 0.) {
	      h_em_mu_hlt_noreco_pt->Fill(trigp4.at(ihlt).pt());
	      h_em_mu_hlt_noreco_eta->Fill(trigp4.at(ihlt).eta());
	    }
	  }
	} // loop on hlt objects

      } // if em trigger

      if (isEM && pass_trig_me) {
	for (unsigned int ilead=0; ilead < mu_lead_flags.size(); ++ilead) {
	  for (unsigned int isubl=0; isubl < el_subl_flags.size(); ++isubl) {
	    h_me_events->Fill(mu_lead_flags.at(ilead),el_subl_flags.at(isubl));
	  }
	}

	// plot hlt obj kinematics
	int trigidx = findTriggerIndex(triggerName(trigname_me));
	std::vector<LorentzVector> trigp4 = hlt_trigObjs_p4()[trigidx];
	std::vector<int> trigid = hlt_trigObjs_id()[trigidx];

	for (unsigned int ihlt = 0; ihlt < trigp4.size(); ++ihlt){
	  if (abs(trigid.at(ihlt)) == 82) { // HLT code for electron
	    h_me_el_hlt_pt->Fill(trigp4.at(ihlt).pt());
	    h_me_el_hlt_eta->Fill(trigp4.at(ihlt).eta());
	    if (pt_lead_el < 0.) {
	      h_me_el_hlt_noreco_pt->Fill(trigp4.at(ihlt).pt());
	      h_me_el_hlt_noreco_eta->Fill(trigp4.at(ihlt).eta());
	    }
	  } else if (abs(trigid.at(ihlt)) == 83) { // HLT code for muon
	    h_me_mu_hlt_pt->Fill(trigp4.at(ihlt).pt());
	    h_me_mu_hlt_eta->Fill(trigp4.at(ihlt).eta());
	    if (pt_lead_mu < 0.) {
	      h_me_mu_hlt_noreco_pt->Fill(trigp4.at(ihlt).pt());
	      h_me_mu_hlt_noreco_eta->Fill(trigp4.at(ihlt).eta());
	    }
	  }
	} // loop on hlt objects

      } // if me trigger

      if (isMM && pass_trig_mm) {
	for (unsigned int ilead=0; ilead < mu_lead_flags.size(); ++ilead) {
	  for (unsigned int isubl=0; isubl < mu_subl_flags.size(); ++isubl) {
	    h_mm_events->Fill(mu_lead_flags.at(ilead),mu_subl_flags.at(isubl));
	  }
	}
      }

      if (isMM && pass_trig_mmtk) {
	for (unsigned int ilead=0; ilead < mu_lead_flags.size(); ++ilead) {
	  for (unsigned int isubl=0; isubl < mu_subl_flags.size(); ++isubl) {
	    h_mmtk_events->Fill(mu_lead_flags.at(ilead),mu_subl_flags.at(isubl));
	  }
	}
      }

      if (isMM && (pass_trig_mm || pass_trig_mmtk)) {
	for (unsigned int ilead=0; ilead < mu_lead_flags.size(); ++ilead) {
	  for (unsigned int isubl=0; isubl < mu_subl_flags.size(); ++isubl) {
	    h_mmor_events->Fill(mu_lead_flags.at(ilead),mu_subl_flags.at(isubl));
	  }
	}

	// plot trigger overlap
	if (pass_trig_mm && !pass_trig_mmtk) h_mm_overlap->Fill(0);
	else if (!pass_trig_mm && pass_trig_mmtk) h_mm_overlap->Fill(1);
	else if (pass_trig_mm && pass_trig_mmtk) h_mm_overlap->Fill(2);
      }

      ++nEventsPass;

      //      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << endl;
  cout << "Processed events: " << nEventsTotal << endl;
  cout << "Events before reco cuts: " << nEventsPreReco << endl;
  cout << "Passed events: " << nEventsPass << endl;
  cout << endl;

  closeOutput();
  //  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal) 
    std::cout << "ERROR: number of events from files (" << nEventsChain 
	      << ") is not equal to total number of processed events (" << nEventsTotal << ")" << std::endl;
  
  return 0;

}


//--------------------------------------------------------------------
 
void dilepStudyLooper::BookHistos(const TString& prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  // rootdir->cd();
  if (outFile) outFile->cd();

  const int max_nvtx = 40;
  const int nbins_el = 8;
  const int nbins_mu = 8;

  h_nvtx = new TH1F(Form("%s_nvtx",prefix.Data()),";N(vtx)",max_nvtx,0,max_nvtx);

  h_el_pt = new TH1F(Form("%s_el_pt",prefix.Data()),";electron p_{T} [GeV]",100,0.,100.);
  h_el_lead_pt = new TH1F(Form("%s_el_lead_pt",prefix.Data()),";leading electron p_{T} [GeV]",100,0.,100.);
  h_el_subl_pt = new TH1F(Form("%s_el_subl_pt",prefix.Data()),";subleading electron p_{T} [GeV]",100,0.,100.);
  h_el_lead_eta = new TH1F(Form("%s_el_lead_eta",prefix.Data()),";leading electron #eta",100,-3.,3.);
  h_el_subl_eta = new TH1F(Form("%s_el_subl_eta",prefix.Data()),";subleading electron #eta",100,-3.,3.);
  h_ee_mll = new TH1F(Form("%s_ee_mll",prefix.Data()),";m_{ee} [GeV]",150,0.,150.);
  h_ee_denom_mll = new TH1F(Form("%s_ee_denom_mll",prefix.Data()),";m_{ee} [GeV]",150,0.,150.);
  h_ee_mediso_mll = new TH1F(Form("%s_ee_mediso_mll",prefix.Data()),";m_{ee} [GeV]",150,0.,150.);

  h_mu_pt = new TH1F(Form("%s_mu_pt",prefix.Data()),";muon p_{T} [GeV]",100,0.,100.);
  h_mu_lead_pt = new TH1F(Form("%s_mu_lead_pt",prefix.Data()),";leading muon p_{T} [GeV]",100,0.,100.);
  h_mu_subl_pt = new TH1F(Form("%s_mu_subl_pt",prefix.Data()),";subleading muon p_{T} [GeV]",100,0.,100.);
  h_mu_lead_eta = new TH1F(Form("%s_mu_lead_eta",prefix.Data()),";leading muon #eta",100,-3.,3.);
  h_mu_subl_eta = new TH1F(Form("%s_mu_subl_eta",prefix.Data()),";subleading muon #eta",100,-3.,3.);
  h_mm_mll = new TH1F(Form("%s_mm_mll",prefix.Data()),";m_{#mu#mu} [GeV]",150,0.,150.);
  h_mu_dr = new TH1F(Form("%s_mu_dr",prefix.Data()),";muon #Delta R",650,0.,6.5);
  h_mm_mll_vs_dr = new TH2F(Form("%s_mm_mll_vs_dr",prefix.Data()),";muon #Delta R; m_{#mu#mu} [GeV]",650,0.,6.5,150,0.,150.);
  h_mm_chargeprod = new TH1F(Form("%s_mm_chargeprod",prefix.Data()),";Q(mu1)*Q(mu2)",3,-1.,2.);

  h_mu_tiso_lead_pt = new TH1F(Form("%s_mu_tiso_lead_pt",prefix.Data()),";leading muon p_{T} [GeV]",100,0.,100.);
  h_mu_tiso_subl_pt = new TH1F(Form("%s_mu_tiso_subl_pt",prefix.Data()),";subleading muon p_{T} [GeV]",100,0.,100.);
  h_mu_tiso_lead_eta = new TH1F(Form("%s_mu_tiso_lead_eta",prefix.Data()),";leading muon #eta",100,-3.,3.);
  h_mu_tiso_subl_eta = new TH1F(Form("%s_mu_tiso_subl_eta",prefix.Data()),";subleading muon #eta",100,-3.,3.);
  h_mm_tiso_mll = new TH1F(Form("%s_mm_tiso_mll",prefix.Data()),";m_{#mu#mu} [GeV]",150,0.,150.);

  h_mu_dup_dpt = new TH1F(Form("%s_mu_dup_dpt",prefix.Data()),";muon #Delta p_{T} [GeV]",100,-50.,50.);
  h_mu_dup_deta = new TH1F(Form("%s_mu_dup_deta",prefix.Data()),";muon #Delta #eta",50,-0.2,0.2);
  h_mu_dup_dphi = new TH1F(Form("%s_mu_dup_dphi",prefix.Data()),";muon #Delta #phi",50,0.0,0.5);
  h_mu_dup_dr = new TH1F(Form("%s_mu_dup_dr",prefix.Data()),";muon #Delta R",50,0.,0.5);
  h_mu_dup_lead_type = new TH1F(Form("%s_mu_dup_lead_type",prefix.Data()),";leading muon type",15,-0.5,14.5);
  h_mu_dup_subl_type = new TH1F(Form("%s_mu_dup_subl_type",prefix.Data()),";subleading muon type",15,-0.5,14.5);

  h_mu_pfchiso = new TH1F(Form("%s_mu_pfchiso",prefix.Data()),";muon PF charged iso / pt",1000,0,2.0);
  h_mu_pfchiso_minus_chiso = new TH1F(Form("%s_mu_pfchiso_minus_chiso",prefix.Data()),";muon PF charged iso - recomputed charged iso [GeV]",1000,-5.0,5.0);
  h_mu_chiso04 = new TH1F(Form("%s_mu_chiso04",prefix.Data()),";muon PF charged iso / pt",1000,0,2.0);
  h_mu_chiso07 = new TH1F(Form("%s_mu_chiso07",prefix.Data()),";muon PF charged iso / pt",1000,0,2.0);
  h_mu_chiso10 = new TH1F(Form("%s_mu_chiso10",prefix.Data()),";muon PF charged iso / pt",1000,0,2.0);

  h_em_el_pt = new TH1F(Form("%s_em_el_pt",prefix.Data()),";electron p_{T} [GeV]",100,0.,100.);
  h_em_mu_pt = new TH1F(Form("%s_em_mu_pt",prefix.Data()),";muon p_{T} [GeV]",100,0.,100.);
  h_em_el_eta = new TH1F(Form("%s_em_el_eta",prefix.Data()),";electron #eta",100,-3.,3.);
  h_em_mu_eta = new TH1F(Form("%s_em_mu_eta",prefix.Data()),";muon #eta",100,-3.,3.);
  h_em_mll = new TH1F(Form("%s_em_mll",prefix.Data()),";m_{e#mu} [GeV]",150,0.,150.);
  h_em_dr = new TH1F(Form("%s_em_dr",prefix.Data()),";#DeltaR(e,#mu)",600,0.,6.);

  h_me_mu_pt = new TH1F(Form("%s_me_mu_pt",prefix.Data()),";muon p_{T} [GeV]",100,0.,100.);
  h_me_el_pt = new TH1F(Form("%s_me_el_pt",prefix.Data()),";electron p_{T} [GeV]",100,0.,100.);
  h_me_mu_eta = new TH1F(Form("%s_me_mu_eta",prefix.Data()),";muon #eta",100,-3.,3.);
  h_me_el_eta = new TH1F(Form("%s_me_el_eta",prefix.Data()),";electron #eta",100,-3.,3.);
  h_me_mll = new TH1F(Form("%s_me_mll",prefix.Data()),";m_{#mue} [GeV]",150,0.,150.);
  h_me_dr = new TH1F(Form("%s_me_dr",prefix.Data()),";#DeltaR(#mu,e)",600,0.,6.);

  h_em_el_tiso_pt = new TH1F(Form("%s_em_el_tiso_pt",prefix.Data()),";electron p_{T} [GeV]",100,0.,100.);
  h_em_mu_tiso_pt = new TH1F(Form("%s_em_mu_tiso_pt",prefix.Data()),";muon p_{T} [GeV]",100,0.,100.);
  h_em_el_tiso_eta = new TH1F(Form("%s_em_el_tiso_eta",prefix.Data()),";electron #eta",100,-3.,3.);
  h_em_mu_tiso_eta = new TH1F(Form("%s_em_mu_tiso_eta",prefix.Data()),";muon #eta",100,-3.,3.);
  h_em_tiso_mll = new TH1F(Form("%s_em_tiso_mll",prefix.Data()),";m_{e#mu} [GeV]",150,0.,150.);
  h_em_tiso_dr = new TH1F(Form("%s_em_tiso_dr",prefix.Data()),";#DeltaR(e,#mu)",600,0.,6.);

  h_me_mu_tiso_pt = new TH1F(Form("%s_me_mu_tiso_pt",prefix.Data()),";muon p_{T} [GeV]",100,0.,100.);
  h_me_el_tiso_pt = new TH1F(Form("%s_me_el_tiso_pt",prefix.Data()),";electron p_{T} [GeV]",100,0.,100.);
  h_me_mu_tiso_eta = new TH1F(Form("%s_me_mu_tiso_eta",prefix.Data()),";muon #eta",100,-3.,3.);
  h_me_el_tiso_eta = new TH1F(Form("%s_me_el_tiso_eta",prefix.Data()),";electron #eta",100,-3.,3.);
  h_me_tiso_mll = new TH1F(Form("%s_me_tiso_mll",prefix.Data()),";m_{#mue} [GeV]",150,0.,150.);
  h_me_tiso_dr = new TH1F(Form("%s_me_tiso_dr",prefix.Data()),";#DeltaR(#mu,e)",600,0.,6.);

  h_em_el_hlt_pt = new TH1F(Form("%s_em_el_hlt_pt",prefix.Data()),";HLT electron p_{T} [GeV]",100,0.,100.);
  h_em_mu_hlt_pt = new TH1F(Form("%s_em_mu_hlt_pt",prefix.Data()),";HLT muon p_{T} [GeV]",100,0.,100.);
  h_em_el_hlt_eta = new TH1F(Form("%s_em_el_hlt_eta",prefix.Data()),";HLT electron #eta",100,-3.,3.);
  h_em_mu_hlt_eta = new TH1F(Form("%s_em_mu_hlt_eta",prefix.Data()),";HLT muon #eta",100,-3.,3.);

  h_me_mu_hlt_pt = new TH1F(Form("%s_me_mu_hlt_pt",prefix.Data()),";HLT muon p_{T} [GeV]",100,0.,100.);
  h_me_el_hlt_pt = new TH1F(Form("%s_me_el_hlt_pt",prefix.Data()),";HLT electron p_{T} [GeV]",100,0.,100.);
  h_me_mu_hlt_eta = new TH1F(Form("%s_me_mu_hlt_eta",prefix.Data()),";HLT muon #eta",100,-3.,3.);
  h_me_el_hlt_eta = new TH1F(Form("%s_me_el_hlt_eta",prefix.Data()),";HLT electron #eta",100,-3.,3.);

  h_em_el_hlt_noreco_pt = new TH1F(Form("%s_em_el_hlt_noreco_pt",prefix.Data()),";HLT electron p_{T} [GeV]",100,0.,100.);
  h_em_mu_hlt_noreco_pt = new TH1F(Form("%s_em_mu_hlt_noreco_pt",prefix.Data()),";HLT muon p_{T} [GeV]",100,0.,100.);
  h_em_el_hlt_noreco_eta = new TH1F(Form("%s_em_el_hlt_noreco_eta",prefix.Data()),";HLT electron #eta",100,-3.,3.);
  h_em_mu_hlt_noreco_eta = new TH1F(Form("%s_em_mu_hlt_noreco_eta",prefix.Data()),";HLT muon #eta",100,-3.,3.);

  h_me_mu_hlt_noreco_pt = new TH1F(Form("%s_me_mu_hlt_noreco_pt",prefix.Data()),";HLT muon p_{T} [GeV]",100,0.,100.);
  h_me_el_hlt_noreco_pt = new TH1F(Form("%s_me_el_hlt_noreco_pt",prefix.Data()),";HLT electron p_{T} [GeV]",100,0.,100.);
  h_me_mu_hlt_noreco_eta = new TH1F(Form("%s_me_mu_hlt_noreco_eta",prefix.Data()),";HLT muon #eta",100,-3.,3.);
  h_me_el_hlt_noreco_eta = new TH1F(Form("%s_me_el_hlt_noreco_eta",prefix.Data()),";HLT electron #eta",100,-3.,3.);

  h_mm_overlap = new TH1F(Form("%s_mm_overlap",prefix.Data()),";Triggers Passed",3,-0.5,2.5);
  h_mm_overlap->GetXaxis()->SetBinLabel(1,"#mu#mu only");
  h_mm_overlap->GetXaxis()->SetBinLabel(2,"#mu#mutk only");
  h_mm_overlap->GetXaxis()->SetBinLabel(3,"both");

  h_ee_events = new TH2F(Form("%s_ee_events",prefix.Data()),";leading e;subleading e",nbins_el,-0.5,nbins_el-0.5,nbins_el,-0.5,nbins_el-0.5);
  h_em_events = new TH2F(Form("%s_em_events",prefix.Data()),";leading e;subleading #mu",nbins_el,-0.5,nbins_el-0.5,nbins_mu,-0.5,nbins_mu-0.5);
  h_me_events = new TH2F(Form("%s_me_events",prefix.Data()),";leading #mu;subleading e",nbins_mu,-0.5,nbins_mu-0.5,nbins_el,-0.5,nbins_el-0.5);
  h_mm_events = new TH2F(Form("%s_mm_events",prefix.Data()),";leading #mu;subleading #mu",nbins_mu,-0.5,nbins_mu-0.5,nbins_mu,-0.5,nbins_mu-0.5);
  h_mmtk_events = new TH2F(Form("%s_mmtk_events",prefix.Data()),";leading #mu;subleading #mu",nbins_mu,-0.5,nbins_mu-0.5,nbins_mu,-0.5,nbins_mu-0.5);
  h_mmor_events = new TH2F(Form("%s_mmor_events",prefix.Data()),";leading #mu;subleading #mu",nbins_mu,-0.5,nbins_mu-0.5,nbins_mu,-0.5,nbins_mu-0.5);

  labelAxis(h_ee_events,0,11); labelAxis(h_ee_events,1,11); 
  labelAxis(h_em_events,0,11); labelAxis(h_em_events,1,13); 
  labelAxis(h_me_events,0,13); labelAxis(h_me_events,1,11); 
  labelAxis(h_mm_events,0,13); labelAxis(h_mm_events,1,13); 
  labelAxis(h_mmtk_events,0,13); labelAxis(h_mmtk_events,1,13); 
  labelAxis(h_mmor_events,0,13); labelAxis(h_mmor_events,1,13); 

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()


void dilepStudyLooper::makeOutput(const TString& prefix){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  TString tpsuffix = "";
  if( doTenPercent ) tpsuffix = "_tenPercent";

  outFile   = new TFile(Form("output/%s%s.root",prefix.Data(),tpsuffix.Data()), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix.Data(),frsuffix,tpsuffix), "RECREATE");
  //  outFile   = new TFile("baby.root","RECREATE");
  outFile->cd();
  //  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in dilepStudyLooper.h
}

//--------------------------------------------------------------------

float dilepStudyLooper::electronPFiso(const unsigned int index, const bool cor) {
    
  float pt     = cms2.els_p4().at(index).pt();
  float etaAbs = fabs(cms2.els_etaSC().at(index));

  // get effective area
  float AEff = 0.;
  if (etaAbs <= 1.0) AEff = 0.10;
  else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
  else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
  else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
  else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
  else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
  else if (etaAbs > 2.4) AEff = 0.13;

  float pfiso_ch = cms2.els_iso03_pf2012ext_ch().at(index);
  float pfiso_em = cms2.els_iso03_pf2012ext_em().at(index);
  float pfiso_nh = cms2.els_iso03_pf2012ext_nh().at(index);
    
  // rho
  float rhoPrime = std::max(cms2.evt_kt6pf_foregiso_rho(), float(0.0));
  float pfiso_n = pfiso_em + pfiso_nh;
  if (cor) pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
  float pfiso = (pfiso_ch + pfiso_n) / pt;

  return pfiso;
}

//--------------------------------------------------------------------

float dilepStudyLooper::muonPFiso(const unsigned int imu, const bool cor) {
  float chiso = cms2.mus_isoR03_pf_ChargedHadronPt().at(imu);
  float nhiso = cms2.mus_isoR03_pf_NeutralHadronEt().at(imu);
  float emiso = cms2.mus_isoR03_pf_PhotonEt().at(imu);
  float deltaBeta = cms2.mus_isoR03_pf_PUPt().at(imu);
  float pt = cms2.mus_p4().at(imu).pt();

  float absiso = chiso + nhiso + emiso;
  if (cor) absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return (absiso / pt);

}

//--------------------------------------------------------------------

void dilepStudyLooper::labelAxis(TH2F* h, int axis, int lep) {

  TAxis* h_axis;
  if (axis == 0) h_axis = h->GetXaxis();
  else if (axis == 1) h_axis = h->GetYaxis();
  else {
    std::cout << "labelAxis: didn't recognize axis: " << axis << ", returning.." << std::endl;
    return;
  }

  // electrons
  if (lep == 11) {
    h_axis->SetBinLabel(1,"All");
    h_axis->SetBinLabel(2,"Cand");
    h_axis->SetBinLabel(3,"WWDenom");
    h_axis->SetBinLabel(4,"Loose");
    h_axis->SetBinLabel(5,"Med");
    h_axis->SetBinLabel(6,"Iso");
    h_axis->SetBinLabel(7,"L+Iso");
    h_axis->SetBinLabel(8,"M+Iso");
  }

  // muons
  else if (lep == 13) {
    h_axis->SetBinLabel(1,"All");
    h_axis->SetBinLabel(2,"Cand");
    h_axis->SetBinLabel(3,"Tight");
    h_axis->SetBinLabel(4,"Iso < 1.0");
    h_axis->SetBinLabel(5,"Iso < 0.7");
    h_axis->SetBinLabel(6,"Iso < 0.4");
    h_axis->SetBinLabel(7,"Iso < 0.15");
    h_axis->SetBinLabel(8,"T+Iso");
  }

  else {
    std::cout << "labelAxis: didn't recognize lep: " << lep << ", returning.." << std::endl;
    return;
  }

  return;
}

//--------------------------------------------------------------------
// WW definitions - taken from HWW2012CORE
//--------------------------------------------------------------------

bool dilepStudyLooper::ElectronFOIdV4(unsigned int i) {

	float pt = cms2.els_p4().at(i).pt();
	float etaSC = cms2.els_etaSC().at(i);

	if (fabs(etaSC)<1.479) {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.01		||
			fabs(cms2.els_dEtaIn().at(i))>0.007 	||
			fabs(cms2.els_dPhiIn().at(i))>0.15 		||
			cms2.els_hOverE().at(i)>0.12 			||
			cms2.els_tkIso().at(i)/pt>0.2 			||
			(cms2.els_ecalIso().at(i) - 1.0)/pt>0.2 ||
			cms2.els_hcalIso().at(i)/pt>0.2 ) return false;
	} else {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.03		|| 
			fabs(cms2.els_dEtaIn().at(i))>0.009 	||
			fabs(cms2.els_dPhiIn().at(i))>0.10 		|| 
			cms2.els_hOverE().at(i)>0.10 			||
			cms2.els_tkIso().at(i)/pt>0.2 			||
		    	cms2.els_ecalIso().at(i)/pt>0.2 		||
			cms2.els_hcalIso().at(i)/pt>0.2 ) return false;
	}

    // MIT conversion
	if ( isFromConversionMIT(i) ) return false;
	// conversion rejection - hit based
	if ( cms2.els_exp_innerlayers().at(i) > 0 ) return false;
	
	return true;
} 

bool dilepStudyLooper::ElectronFOV4(unsigned int i){
    return ww_elBase(i) && ElectronFOIdV4(i) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool dilepStudyLooper::ww_elBase(unsigned int index){
    if (cms2.els_p4().at(index).pt() < ptthresh_low) return false;
    if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
    return true;
}

double dilepStudyLooper::dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

bool dilepStudyLooper::ww_eld0PV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    double dxyPV = cms2.els_d0()[index]-
        cms2.vtxs_position()[vtxIndex].x()*sin(cms2.els_trk_p4()[index].phi())+
        cms2.vtxs_position()[vtxIndex].y()*cos(cms2.els_trk_p4()[index].phi());
    return fabs(dxyPV) < 0.02;
}

bool dilepStudyLooper::ww_eldZPV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    // double dzPV = cms2.els_z0corr()[index]-cms2.vtxs_position()[iMax].z();
    double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
    return fabs(dzpv)<0.1;
}

int dilepStudyLooper::primaryVertex() {
  return 0;
}

                                                               
//--------------------------------------------------------------------
// object-trigger matching: taken from SingleLepton2012/looper/stopUtils.cc
//--------------------------------------------------------------------

int dilepStudyLooper::findTriggerIndex(const TString& trigName)
{
  vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
  if(found_it != end_it) return found_it - begin_it;
  return -1;
}

//--------------------------------------------------------------------

TString dilepStudyLooper::triggerName(const TString& triggerPattern){

  //-------------------------------------------------------
  // get exact trigger name corresponding to given pattern
  //-------------------------------------------------------

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger) return "TRIGGER_NOT_FOUND";

  return exact_hltname;

}


bool dilepStudyLooper::objectPassTrigger(const LorentzVector &obj, const TString& trigname, int type, float drmax ){

  if (type != 82 && type != 83) {
    cout << __FILE__ << " " << __LINE__ << " Error! invalid HLT object type: " << type << endl;
    return false;
  }

  TString exact_trigname = triggerName( trigname );

  if( exact_trigname.Contains("TRIGGER_NOT_FOUND") ){
    cout << __FILE__ << " " << __LINE__ << " Error! couldn't find trigger name " << trigname << endl;
    return false;
  }

  int trigidx = findTriggerIndex(exact_trigname);
  std::vector<LorentzVector> trigp4 = hlt_trigObjs_p4()[trigidx];
  std::vector<int> trigid = hlt_trigObjs_id()[trigidx];

  if( trigp4.size() == 0 ) return false;

  for (unsigned int i = 0; i < trigp4.size(); ++i){
    if (trigid.at(i) != type) continue;
    float dr = dRbetweenVectors(trigp4[i], obj);
    if( dr < drmax ) return true;
  }

  return false;
}

//----------------------------

float dilepStudyLooper::getdphi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//-------------------------------------------------------------------

// based on muonIsoValuePF2012 from muonSelections.cc

isovals dilepStudyLooper::muonChIsoValuePF2012 (const unsigned int imu, const float R, const int ivtx)
{

  // isolation sums
  isovals vals;
  vals.chiso00 = 0.;
  vals.chiso04 = 0.;
  vals.chiso07 = 0.;
  vals.chiso10 = 0.;
       
  // loop on pfcandidates
  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ++ipf) {
            
    // skip electrons and muons
    const int particleId = abs(cms2.pfcands_particleId()[ipf]);
    if (particleId == 11)    continue;
    if (particleId == 13)    continue;

    // deltaR between electron and cadidate
    const float dR = dRbetweenVectors(cms2.pfcands_p4()[ipf], cms2.mus_p4()[imu]);
    if (dR > R)              continue;

    // charged hadrons closest vertex
    // should be the primary vertex
    if (cms2.pfcands_charge().at(ipf) != 0) {
      //        if (particleId == 211 || particleId == 321 || particleId == 2212 || particleId == 999211) {
      if (cms2.pfcands_vtxidx().at(ipf) != ivtx) continue;
      if (dR < 0.0001)
	continue;

      float pt = cms2.pfcands_p4()[ipf].pt();
      vals.chiso00 += pt;
      if (pt > 0.4) vals.chiso04 += pt;
      if (pt > 0.7) vals.chiso07 += pt;
      if (pt > 1.0) vals.chiso10 += pt;
    }

  } // loop on pf cands

  return vals;
}

