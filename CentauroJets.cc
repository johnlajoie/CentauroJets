
#include "CentauroJets.h"

#include <phool/phool.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4VtxPointv1.h>

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>

#include <HepMC/GenEvent.h>
#include <HepMC/SimpleVector.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <phgeom/PHGeomUtility.h>

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/Centauro.hh>

#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

// GEANT
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>

#include <iostream>
#include <limits>

#include <utility>

#define LogDebug(exp)		std::cout<<"DEBUG: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogError(exp)		std::cout<<"ERROR: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWarning(exp)	std::cout<<"WARNING: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"

using namespace std;
using namespace fastjet; 

// Tower/cluster cutoffs
#define TOWER_E_CUT 0.200 
#define CLUSTER_E_CUTOFF 0.100

// Cluster/track matching cuts
#define BECAL_CLUST_TRACKMATCH 0.10
#define IHCAL_CLUST_TRACKMATCH 0.35
#define OHCAL_CLUST_TRACKMATCH 0.35
#define FEMC_CLUST_TRACKMATCH 0.15
#define LFHCAL_CLUST_TRACKMATCH 0.35

// HCAL neutral energy scales
#define BARREL_HCAL_NEUT_SCALE (1.0/0.46)
#define FWD_HCAL_NEUT_SCALE (1.0/0.77)

double XYtoPhi(double x, double y)
{
  // -Pi to +Pi
  Double_t phi = atan2(y,x);
  if(phi<-TMath::Pi()) phi += TMath::TwoPi();
  if(phi>=TMath::Pi()) phi -= TMath::TwoPi();
  return phi;
}

double XYtoPhi_02PI(double x, double y)
{
  // 0 to 2pi
  Double_t phi = atan2(y,x);
  if(phi<0.0) phi += TMath::TwoPi();
  if(phi>=TMath::TwoPi()) phi -= TMath::TwoPi();
  return phi;
}

double getEta(double pt, double pz){
  double theta = XYtoPhi(pz,pt);
  double eta = -log(tan(theta/2.0));
  return eta; 
} 

double DeltaPhi(double phi1, double phi2){
  Double_t dphi;
  dphi = phi1 - phi2;
  if(dphi<-TMath::Pi()) dphi += TMath::TwoPi();
  if(dphi>=TMath::Pi()) dphi -= TMath::TwoPi();
  return dphi;
}

int EncodeUserIndex(int charge, bool emPart)
{

  int uidx = 0;
  if(charge>=1) 
    uidx = 1; 
  else if(charge<=-1) 
    uidx = 2; 

  if(emPart) uidx += 10; 

  return uidx; 
 
}

int GetChargeFromUserIndex(int uidx)
{

  int  charge = uidx; 
  if(charge >=10) charge -= 10;
  if(charge == 2) charge = -1; 

  return charge; 
 
}

bool isEMParticle(int uidx)
{

  if(uidx>=10) 
    return true;
  else
    return false; 

}

double JetCharge( std::vector<fastjet::PseudoJet> *tconstit, double jetP ){

  // Calculate jet charge - arXiv:1209.2421v2
  // since the jet is mostly along the z-axis for the current jets, 
  // use p instead of pT

  const double kappa = 0.5; 
  double jc = 0; 

  for(unsigned int i=0; i<tconstit->size(); i++){

    int charge = GetChargeFromUserIndex(tconstit->at(i).user_index()); 

    if(charge==0) continue;
    double cptot = sqrt( pow(tconstit->at(i).px(),2) + 
			 pow(tconstit->at(i).py(),2) +
			 pow(tconstit->at(i).pz(),2) ); 
    jc += (charge)*pow(cptot, kappa); 
  }

  jc /= pow(jetP,kappa); 

  return jc; 

}

double JetChargedFraction( std::vector<fastjet::PseudoJet> *tconstit, double jetP ){

  double cpx = 0.0; 
  double cpy = 0.0; 
  double cpz = 0.0; 

  for(unsigned int i=0; i<tconstit->size(); i++){
    // eliminate neutrals
    if(GetChargeFromUserIndex(tconstit->at(i).user_index())==0) continue;   
    // eliminate electrons
    if( (tconstit->at(i).user_index()==11) || (tconstit->at(i).user_index()==12) ) continue; 
    cpx += tconstit->at(i).px(); 
    cpy += tconstit->at(i).py(); 
    cpz += tconstit->at(i).pz(); 
   }

  return sqrt(cpx*cpx + cpy*cpy + cpz*cpz)/jetP; 

}

double JetNeutralMomentum( std::vector<fastjet::PseudoJet> *tconstit ){

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    // eliminate charged
    if(GetChargeFromUserIndex(tconstit->at(i).user_index())!=0) continue;
    // eliminate photons
    if(tconstit->at(i).user_index()==10) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot.Mag(); 

}

double JetChargedMomentum( std::vector<fastjet::PseudoJet> *tconstit ){

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    if(GetChargeFromUserIndex(tconstit->at(i).user_index())==0) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot.Mag(); 

}

double JetEMMomentum( std::vector<fastjet::PseudoJet> *tconstit ){

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    if(isEMParticle(tconstit->at(i).user_index())) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot.Mag(); 

}


//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
CentauroJets::CentauroJets(const string &name) : SubsysReco(name) {
	//initialize
	_event = 0;
	_outfile_name = name;
	rand = (TRandom *)new TRandom3(0); 

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int CentauroJets::Init(PHCompositeNode *topNode) {
	cout << PHWHERE << " Opening file " << _outfile_name << endl;
	PHTFileServer::get().open(_outfile_name, "RECREATE");

	// Jet evaluation

	_eval_tree_event = new TTree("event", "CentauroJets");
	_eval_tree_event->Branch("event", &event, "event/I");
	_eval_tree_event->Branch("x1", &_hepmcp_x1, "_hepmcp_x1/D");
	_eval_tree_event->Branch("x2", &_hepmcp_x2, "_hepmcp_x2/D");
	_eval_tree_event->Branch("Q2", &_hepmcp_Q2, "_hepmcp_Q2/D");
	_eval_tree_event->Branch("procid", &_hepmcp_procid, "_hepmcp_procid/I");
	_eval_tree_event->Branch("id1", &_hepmcp_id1, "_hepmcp_id1/I");
	_eval_tree_event->Branch("id2", &_hepmcp_id2, "_hepmcp_id2/I");
	_eval_tree_event->Branch("vtx_x", &vtx_x, "vtx_x/D");
	_eval_tree_event->Branch("vtx_y", &vtx_y, "vtx_y/D");
	_eval_tree_event->Branch("vtx_z", &vtx_z, "vtx_z/D");
	_eval_tree_event->Branch("vtx_t", &vtx_t, "vtx_t/D");	
	_eval_tree_event->Branch("measQ2", &measQ2, "measQ2/D");
	_eval_tree_event->Branch("meas_x", &meas_x, "meas_x/D");
	_eval_tree_event->Branch("meas_E_p", &meas_E_p, "meas_E_p/D"); 
	_eval_tree_event->Branch("electron_eta",&electron_eta,"electron_eta/D"); 
	_eval_tree_event->Branch("electron_cluster_dR",&electron_cluster_dR,"electron_cluster_dR/D"); 

	_eval_tree_event->Branch("breit_vphot_e",&breit_vphot_e,"breit_vphot_e/D"); 
	_eval_tree_event->Branch("breit_vphot_pz",&breit_vphot_pz,"breit_vphot_pz/D"); 
	_eval_tree_event->Branch("breit_vphot_pt",&breit_vphot_pt,"breit_vphot_pt/D"); 
	_eval_tree_event->Branch("breit_initial_proton_e",&breit_initial_proton_e,"breit_initial_proton_e/D"); 
	_eval_tree_event->Branch("breit_initial_proton_pz",&breit_initial_proton_pz,"breit_initial_proton_pz/D"); 
	_eval_tree_event->Branch("breit_initial_proton_pt",&breit_initial_proton_pt,"breit_initial_proton_pt/D"); 

	_eval_tree_event->Branch("jet_pT",&jet_pT); 
	_eval_tree_event->Branch("jet_p",&jet_p); 
	_eval_tree_event->Branch("jet_E",&jet_E); 
	_eval_tree_event->Branch("jet_eta",&jet_eta); 
	_eval_tree_event->Branch("jet_phi",&jet_phi); 
	_eval_tree_event->Branch("jet_z",&jet_z);
	_eval_tree_event->Branch("jet_qperp",&jet_qperp);
	_eval_tree_event->Branch("jet_largest",&jet_largest);
	_eval_tree_event->Branch("jet_pidx",&jet_pidx);
	_eval_tree_event->Branch("jet_pdR",&jet_pdR);
	_eval_tree_event->Branch("jet_lab_eta",&jet_lab_eta); 
	_eval_tree_event->Branch("jet_lab_phi",&jet_lab_phi); 
	_eval_tree_event->Branch("jet_nc",&jet_nc);

	_eval_tree_event->Branch("tjet_pT",&tjet_pT); 
	_eval_tree_event->Branch("tjet_p",&tjet_p); 
	_eval_tree_event->Branch("tjet_E",&tjet_E); 
	_eval_tree_event->Branch("tjet_eta",&tjet_eta); 
	_eval_tree_event->Branch("tjet_phi",&tjet_phi); 
	_eval_tree_event->Branch("tjet_z",&tjet_z);
	_eval_tree_event->Branch("tjet_qperp",&tjet_qperp);
	_eval_tree_event->Branch("tjet_largest",&tjet_largest);
	_eval_tree_event->Branch("tjet_pidx",&tjet_pidx);
	_eval_tree_event->Branch("tjet_pdR",&tjet_pdR);
	_eval_tree_event->Branch("tjet_lab_eta",&tjet_lab_eta); 
	_eval_tree_event->Branch("tjet_lab_phi",&tjet_lab_phi); 
	_eval_tree_event->Branch("tjet_nc",&tjet_nc);
	_eval_tree_event->Branch("tjet_Q",&tjet_Q);

	_eval_tree_event->Branch("tcjet_pT",&tcjet_pT); 
	_eval_tree_event->Branch("tcjet_p",&tcjet_p); 
	_eval_tree_event->Branch("tcjet_E",&tcjet_E); 
	_eval_tree_event->Branch("tcjet_eta",&tcjet_eta); 
	_eval_tree_event->Branch("tcjet_phi",&tcjet_phi); 
	_eval_tree_event->Branch("tcjet_z",&tcjet_z);
	_eval_tree_event->Branch("tcjet_qperp",&tcjet_qperp);
	_eval_tree_event->Branch("tcjet_largest",&tcjet_largest);
	_eval_tree_event->Branch("tcjet_pidx",&tcjet_pidx);
	_eval_tree_event->Branch("tcjet_pdR",&tcjet_pdR);
	_eval_tree_event->Branch("tcjet_lab_eta",&tcjet_lab_eta); 
	_eval_tree_event->Branch("tcjet_lab_phi",&tcjet_lab_phi); 
	_eval_tree_event->Branch("tcjet_nc",&tcjet_nc);
	_eval_tree_event->Branch("tcjet_Q",&tcjet_Q);
	_eval_tree_event->Branch("tcjet_cf",&tcjet_cf);
	_eval_tree_event->Branch("tcjet_neut_p",&tcjet_neut_p); 
	_eval_tree_event->Branch("tcjet_chgd_p",&tcjet_chgd_p); 
	_eval_tree_event->Branch("tcjet_em_p",&tcjet_em_p); 

	_eval_tree_event->Branch("pjet_pT",&pjet_pT); 
	_eval_tree_event->Branch("pjet_p",&pjet_p); 
	_eval_tree_event->Branch("pjet_E",&pjet_E); 
	_eval_tree_event->Branch("pjet_eta",&pjet_eta); 
	_eval_tree_event->Branch("pjet_phi",&pjet_phi); 
	_eval_tree_event->Branch("pjet_z",&pjet_z);
	_eval_tree_event->Branch("pjet_qperp",&pjet_qperp);
	_eval_tree_event->Branch("pjet_largest",&pjet_largest);
	_eval_tree_event->Branch("pjet_nc",&pjet_nc);
	_eval_tree_event->Branch("pjet_Q",&pjet_Q);
	_eval_tree_event->Branch("pjet_cf",&pjet_cf);
	_eval_tree_event->Branch("pjet_tcidx",&pjet_tcidx);
	_eval_tree_event->Branch("pjet_tcdR",&pjet_tcdR);
	_eval_tree_event->Branch("pjet_neut_p",&pjet_neut_p); 
	_eval_tree_event->Branch("pjet_chgd_p",&pjet_chgd_p); 
	_eval_tree_event->Branch("pjet_em_p",&pjet_em_p); 
 
	_eval_tree_event->Branch("tfpjet_pT",&tfpjet_pT); 
	_eval_tree_event->Branch("tfpjet_p",&tfpjet_p); 
	_eval_tree_event->Branch("tfpjet_E",&tfpjet_E); 
	_eval_tree_event->Branch("tfpjet_eta",&tfpjet_eta); 
	_eval_tree_event->Branch("tfpjet_phi",&tfpjet_phi); 
	_eval_tree_event->Branch("tfpjet_z",&tfpjet_z);
	_eval_tree_event->Branch("tfpjet_qperp",&tfpjet_qperp);
	_eval_tree_event->Branch("tfpjet_largest",&tfpjet_largest);
 	_eval_tree_event->Branch("tfpjet_nc",&tfpjet_nc);
	_eval_tree_event->Branch("tfpjet_Q",&tfpjet_Q);
	_eval_tree_event->Branch("tfpjet_cf",&tfpjet_cf);
	_eval_tree_event->Branch("tfpjet_neut_p",&tfpjet_neut_p); 
	_eval_tree_event->Branch("tfpjet_chgd_p",&tfpjet_chgd_p); 
	_eval_tree_event->Branch("tfpjet_em_p",&tfpjet_em_p); 

	// cluster evaluation for charged tracks
	_eval_charged_tracks_cent = new TTree("clusteval_cent", "Charged track clusters (central)");
	_eval_charged_tracks_cent->Branch("event", &event, "event/I");
	_eval_charged_tracks_cent->Branch("pid",&ct_pid); 
	_eval_charged_tracks_cent->Branch("p_meas",&ct_p_meas); 
	_eval_charged_tracks_cent->Branch("p_true",&ct_p_true); 
	_eval_charged_tracks_cent->Branch("eta_meas",&ct_eta_meas); 
	_eval_charged_tracks_cent->Branch("eta_true",&ct_eta_true); 
	_eval_charged_tracks_cent->Branch("dist",&ct_dist); 
	_eval_charged_tracks_cent->Branch("e_bemc",&ct_e_bemc); 
	_eval_charged_tracks_cent->Branch("e_ihcal",&ct_e_ihcal); 
	_eval_charged_tracks_cent->Branch("e_ohcal",&ct_e_ohcal); 
	_eval_charged_tracks_cent->Branch("e_tot",&ct_e_tot); 

	_eval_charged_tracks_fwd = new TTree("clusteval_fwd", "Charged track clusters (fwd)");
	_eval_charged_tracks_fwd->Branch("event", &event, "event/I");
	_eval_charged_tracks_fwd->Branch("pid",&ct_pid); 
	_eval_charged_tracks_fwd->Branch("p_meas",&ct_p_meas); 
	_eval_charged_tracks_fwd->Branch("p_true",&ct_p_true); 
	_eval_charged_tracks_fwd->Branch("eta_meas",&ct_eta_meas); 
	_eval_charged_tracks_fwd->Branch("eta_true",&ct_eta_true); 
	_eval_charged_tracks_fwd->Branch("dist",&ct_dist); 
	_eval_charged_tracks_fwd->Branch("e_femc",&ct_e_femc); 
	_eval_charged_tracks_fwd->Branch("e_lfhcal",&ct_e_lfhcal); 
	_eval_charged_tracks_fwd->Branch("e_tot",&ct_e_tot); 

	// evaluation for calo tracks
	_eval_calo_tracks_cent = new TTree("calotrackeval_cent", "Calotrack Evaluation (central)");
	_eval_calo_tracks_cent->Branch("event", &event, "event/I");
	_eval_calo_tracks_cent->Branch("pid",&cat_pid); 
	_eval_calo_tracks_cent->Branch("p_true",&cat_p_true); 
	_eval_calo_tracks_cent->Branch("eta_meas",&cat_eta_meas); 
	_eval_calo_tracks_cent->Branch("eta_true",&cat_eta_true); 
	_eval_calo_tracks_cent->Branch("match",&cat_match); 
	_eval_calo_tracks_cent->Branch("e_tot",&cat_e_tot); 

	_eval_calo_tracks_fwd = new TTree("calotrackeval_fwd", "Calotrack Evaluation (fwd)");
	_eval_calo_tracks_fwd->Branch("event", &event, "event/I");
	_eval_calo_tracks_fwd->Branch("pid",&cat_pid); 
	_eval_calo_tracks_fwd->Branch("p_true",&cat_p_true); 
	_eval_calo_tracks_fwd->Branch("eta_meas",&cat_eta_meas); 
	_eval_calo_tracks_fwd->Branch("eta_true",&cat_eta_true); 
	_eval_calo_tracks_fwd->Branch("match",&cat_match); 
	_eval_calo_tracks_fwd->Branch("e_tot",&cat_e_tot); 

	// Diagnostic histograms

	_h_track_cluster_match = new TH1D("_h_track_cluster_match","",200,0.0,1.0); 
	_h_track_cluster_match_becal = new TH1D("_h_track_cluster_match_becal","",200,0.0,1.0); 
	_h_track_cluster_match_ihcal = new TH1D("_h_track_cluster_match_ihcal","",200,0.0,1.0); 
	_h_track_cluster_match_ohcal = new TH1D("_h_track_cluster_match_ohcal","",200,0.0,1.0); 
	_h_track_cluster_match_femc = new TH1D("_h_track_cluster_match_femc","",200,0.0,1.0); 
	_h_track_cluster_match_lfhcal = new TH1D("_h_track_cluster_match_lfhcal","",200,0.0,1.0); 

	_h_becal_ihcal_match = new TH1D("_h_becal_ihcal_match","",200,0.0,5.0); 
	_h_becal_ihcal_match_eta = new TH1D("_h_becal_ihcal_match_eta","",200,0.0,6.5); 
	_h_becal_ihcal_match_phi = new TH1D("_h_becal_ihcal_match_phi","",200,0.0,6.5); 
	_h_becal_ohcal_match = new TH1D("_h_becal_ohcal_match","",200,0.0,5.0); 
	_h_becal_ohcal_match_eta = new TH1D("_h_becal_ohcal_match_eta","",200,0.0,5.0); 
	_h_becal_ohcal_match_phi = new TH1D("_h_becal_ohcal_match_phi","",200,0.0,5.0); 
	_h_ihcal_ohcal_match = new TH1D("_h_ihcal_ohcal_match","",200,0.0,5.0); 
	_h_ihcal_ohcal_match_eta = new TH1D("_h_ihcal_ohcal_match_eta","",200,0.0,5.0); 
	_h_ihcal_ohcal_match_phi = new TH1D("_h_ihcal_ohcal_match_phi","",200,0.0,5.0); 
	_h_femc_lfhcal_match = new TH1D("_h_femc_lfhcal_match","",200,0.0,5.0); 
	_h_femc_lfhcal_match_eta = new TH1D("_h_femc_lfhcal_match_eta","",200,0.0,5.0); 
	_h_femc_lfhcal_match_phi = new TH1D("_h_femc_lfhcal_match_phi","",200,0.0,5.0); 

	_h_calotrack_prim_match_cent = new TH1D("_h_calotrack_prim_match_cent","",400,0.0,2.0); 
	_h_calotrack_prim_match_cent_gamma = new TH1D("_h_calotrack_prim_match_cent_gamma","",400,0.0,2.0); 
	_h_calotrack_prim_match_cent_neutron = new TH1D("_h_calotrack_prim_match_cent_neutron","",400,0.0,2.0); 
	_h_calotrack_pid_prim_match_cent = new TH1D("_h_calotrack_pid_prim_match_cent","",6001,-3000.5,3000.5);

	_h_calotrack_prim_match_fwd = new TH1D("_h_calotrack_prim_match_fwd","",400,0.0,2.0); 
	_h_calotrack_prim_match_fwd_gamma = new TH1D("_h_calotrack_prim_match_fwd_gamma","",400,0.0,2.0); 
	_h_calotrack_prim_match_fwd_neutron = new TH1D("_h_calotrack_prim_match_fwd_neutron","",400,0.0,2.0); 
	_h_calotrack_pid_prim_match_fwd = new TH1D("_h_calotrack_pid_prim_match_fwd","",6001,-3000.5,3000.5);

	_h_calotrack_cent_NES = new TH1D("_h_calotrack_cent_NES","",200,0.0,2.0); 
	_h_calotrack_fwd_NES = new TH1D("_h_calotrack_fwd_NES","",200,0.0,2.0); 

	_h_calotrack_cent_NES_2D = new TH2D("_h_calotrack_cent_NES_2D","",80,0.0,40.0,200,0.0,2.0); 
	_h_calotrack_fwd_NES_2D = new TH2D("_h_calotrack_fwd_NES_2D","",80,0.0,40.0,200,0.0,2.0); 

	_h_calotrack_cent_gamma_NES_2D = new TH2D("_h_calotrack_cent_gamma_NES_2D","",80,0.0,40.0,200,0.0,2.0); 
	_h_calotrack_fwd_gamma_NES_2D = new TH2D("_h_calotrack_fwd_gamma_NES_2D","",80,0.0,40.0,200,0.0,2.0); 

	_h_calotrack_cent_neutron_NES_2D = new TH2D("_h_calotrack_cent_neutron_NES_2D","",80,0.0,40.0,200,0.0,2.0); 
	_h_calotrack_fwd_neutron_NES_2D = new TH2D("_h_calotrack_fwd_neutron_NES_2D","",80,0.0,40.0,200,0.0,2.0); 

	_h_nclusters_becal = new TH1D("_h_nclusters_becal","",51,-0.5,50.5); 
	_h_nclusters_ihcal = new TH1D("_h_nclusters_ihcal","",51,-0.5,50.5); 
	_h_nclusters_ohcal = new TH1D("_h_nclusters_ohcal","",51,-0.5,50.5); 
	_h_nclusters_femc = new TH1D("_h_nclusters_femc","",51,-0.5,50.5); 
	_h_nclusters_lfhcal = new TH1D("_h_nclusters_lfhcal","",51,-0.5,50.5); 

	return Fun4AllReturnCodes::EVENT_OK;
}

int CentauroJets::InitRun(PHCompositeNode *topNode) {

	return Fun4AllReturnCodes::EVENT_OK;

}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int CentauroJets::process_event(PHCompositeNode *topNode) {
	_event++;
	if (_event % 1000 == 0)
		cout << PHWHERE << "Events processed: " << _event << endl;

	GetNodes(topNode);

	fill_tree(topNode);

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int CentauroJets::End(PHCompositeNode *topNode) {

	PHTFileServer::get().cd(_outfile_name);

	_eval_tree_event->Write();

	_eval_charged_tracks_cent->Write(); 
	_eval_charged_tracks_fwd->Write(); 

	_eval_calo_tracks_cent->Write(); 
	_eval_calo_tracks_fwd->Write(); 

	_h_track_cluster_match->Write(); 
	_h_track_cluster_match_becal->Write(); 
	_h_track_cluster_match_ihcal->Write(); 
	_h_track_cluster_match_ohcal->Write(); 
	_h_track_cluster_match_femc->Write(); 
	_h_track_cluster_match_lfhcal->Write(); 

	_h_becal_ihcal_match->Write(); 
	_h_becal_ihcal_match_eta->Write(); 
	_h_becal_ihcal_match_phi->Write(); 
	_h_becal_ohcal_match->Write(); 
	_h_becal_ohcal_match_eta->Write(); 
	_h_becal_ohcal_match_phi->Write(); 
	_h_ihcal_ohcal_match->Write(); 
	_h_ihcal_ohcal_match_eta->Write(); 
	_h_ihcal_ohcal_match_phi->Write(); 
	_h_femc_lfhcal_match->Write(); 
	_h_femc_lfhcal_match_eta->Write(); 
	_h_femc_lfhcal_match_phi->Write(); 

	_h_calotrack_prim_match_cent->Write(); 
	_h_calotrack_prim_match_cent_gamma->Write(); 
	_h_calotrack_prim_match_cent_neutron->Write(); 
	_h_calotrack_pid_prim_match_cent->Write(); 

	_h_calotrack_prim_match_fwd->Write(); 
	_h_calotrack_prim_match_fwd_gamma->Write(); 
	_h_calotrack_prim_match_fwd_neutron->Write(); 
	_h_calotrack_pid_prim_match_fwd->Write(); 

	_h_calotrack_cent_NES->Write(); 
	_h_calotrack_fwd_NES->Write(); 

	_h_calotrack_cent_NES_2D->Write(); 
	_h_calotrack_fwd_NES_2D->Write(); 

	_h_calotrack_cent_gamma_NES_2D->Write(); 
	_h_calotrack_fwd_gamma_NES_2D->Write(); 

	_h_calotrack_cent_neutron_NES_2D->Write(); 
	_h_calotrack_fwd_neutron_NES_2D->Write(); 

	_h_nclusters_becal->Write(); 
	_h_nclusters_ihcal->Write(); 
	_h_nclusters_ohcal->Write(); 
	_h_nclusters_femc->Write(); 
	_h_nclusters_lfhcal->Write(); 

	delete rand; 

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void CentauroJets::fill_tree(PHCompositeNode *topNode) {

  if (!_truth_container) {
    LogError("_truth_container not found!");
    return;
  }
	
  if (!_trackmap) {
    LogError("_trackmap not found!");
    return;
  }

  event = _event;

  // DIS kinematics and jet reconstruction

  _hepmcp_x1 = -9999.0; 
  _hepmcp_x2 = -9999.0; 
  _hepmcp_Q2 = -9999.0; 
  _hepmcp_procid = -9999.0; 
  _hepmcp_id1 = -9999; 
  _hepmcp_id2 = -9999; 
	
  e_p_initial = -9999.0; 
  p_p_initial = -9999.0; 

  true_electron_headon = NULL; 
  TLorentzVector true_e_initial; 
  TLorentzVector true_p_initial; 

  bool firstQuark = true; 
	
  PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (hepmceventmap) {
    for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
	 eventIter != hepmceventmap->end();
	 ++eventIter){

      PHHepMCGenEvent* hepmcevent = eventIter->second;

      if (hepmcevent) {

	// Get the transformation from the event generator to the lab frame

	EventToLab = hepmcevent->get_LorentzRotation_EvtGen2Lab(); 

	// Get the truth event
	      
	HepMC::GenEvent* truthevent = hepmcevent->getEvent();
	if (!truthevent){
	  cout << PHWHERE
	       << "no evt pointer under phhepmvgeneventmap found "
	       << endl;
	}
	else{

	  HepMC::PdfInfo* pdfinfo = truthevent->pdf_info();

	  _hepmcp_x1 = pdfinfo->x1();
	  _hepmcp_x2 = pdfinfo->x2();
	  _hepmcp_Q2 = pdfinfo->scalePDF();
	  _hepmcp_procid = truthevent->signal_process_id();

	  // No good - not stored by eicsmear
	  //_hepmcp_id1 = pdfinfo->id1();
	  //_hepmcp_id2 = pdfinfo->id2();

	  std::pair< HepMC::GenParticle*, HepMC::GenParticle * > bpart = truthevent->beam_particles(); 

	  HepMC::FourVector b1 =  bpart.first->momentum(); 
	  HepMC::FourVector b2 =  bpart.second->momentum();

	  // Record the initial beam momentum - needed for the lab kinematics
	  e_p_initial = fabs(b1.pz()); 
	  p_p_initial = fabs(b2.pz()); 

	  // Transform to the lab frame
	  CLHEP::HepLorentzVector hb1(b1.px(), b1.py(), b1.pz(), b1.e());
	  CLHEP::HepLorentzVector hb2(b2.px(), b2.py(), b2.pz(), b2.e());
	  hb1 = EventToLab * hb1; 
	  hb2 = EventToLab * hb2; 

	  true_e_initial = TLorentzVector(hb1.px(), hb1.py(), hb1.pz(), hb1.e()); 
	  true_p_initial = TLorentzVector(hb2.px(), hb2.py(), hb2.pz(), hb2.e());

	}

	// look for the struck quark

	_hepmcp_id1 = 11; // electron

	HepMC::GenEvent::particle_iterator pitr = truthevent->particles_begin();
	HepMC::GenEvent::particle_iterator evtStart; 
	for( ; pitr!=truthevent->particles_end(); pitr++){
	  HepMC::GenParticle *ptcle = *pitr; 

	  if((abs(ptcle->pdg_id())<9)&&firstQuark){
	    _hepmcp_id2 = ptcle->pdg_id(); 
	    firstQuark = false; 
	    //std::cout << "found quark = " << _hepmcp_id2  << std::endl; 
	  }

	  if(!firstQuark) break; 

	}

      }

    }
  }

  // get the vertex

  PHG4VtxPoint *vtxPoint = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
  vtx_x =  vtxPoint->get_x(); 
  vtx_y =  vtxPoint->get_y(); 
  vtx_z =  vtxPoint->get_z(); 
  vtx_t =  vtxPoint->get_t(); 

  // Find the scattered electron 

  PHG4TruthInfoContainer::ConstRange range = _truth_container->GetPrimaryParticleRange();

  double minQ2diff = 1.0; 

  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr) {

    PHG4Particle* g4particle = truth_itr->second;
    if(!g4particle) {
      LogDebug("");
      continue;
    }

    if (g4particle->get_pid() != 11) continue;
 
    TLorentzVector e_i(0.0,0.0,-e_p_initial,e_p_initial); 
    TLorentzVector e_f(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e()); 
    TLorentzVector v_phot = (e_i - e_f); 
    double checkQ2 = -v_phot.Mag2();

    if(fabs(checkQ2-_hepmcp_Q2)<minQ2diff){
      // NOTE: This is in the head-on frame
      true_electron_headon = g4particle;
      minQ2diff = fabs(checkQ2-_hepmcp_Q2); 
    }

  }

  if(!true_electron_headon){
    LogWarning("scattered electron not found in primaries!"); 

    // Diagnostic for single particle events - cluster/track evaluation

    BuildChargedCaloTracks(topNode,"CENT"); 
    BuildChargedCaloTracks(topNode,"FWD"); 

    // Do the CaloTracks - diagnostic for single particle events
    // note that these will be in the lab frame

    std::vector<fastjet::PseudoJet> tcpseudojets;

    TLorentzRotation breit; // initialized as identity
    TRotation breitRotInv;  // initialized as identity

    BuildCaloTracks(topNode, "CENT", tcpseudojets, breit, breitRotInv, "", -1); 
    BuildCaloTracks(topNode, "FWD", tcpseudojets, breit, breitRotInv, "", -1); 

    return; 
  }

  // Find the scattered electron in the reconstructed data

  SvtxTrack* electron = NULL;	
  for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
       track_itr != _trackmap->end(); track_itr++) {

    SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

    if ((temp->get_truth_track_id() - true_electron_headon->get_track_id()) == 0) {
      electron = temp;
      break; 
    }
  }
	  
  // What did we find? 

  meas_x = -9999.0; 
  measQ2 = -9999.0; 
  meas_E_p = -9999.0;
  electron_eta = -9999.0; 
  electron_cluster_dR = -9999.0; 
	
  breit_vphot_e = -9999.0; 
  breit_initial_proton_e = -9999.0;
  breit_vphot_pz = -9999.0; 
  breit_vphot_pt = -9999.0;
  breit_initial_proton_pz = -9999.0; 
  breit_initial_proton_pt = -9999.0; 

  jet_pT.clear(); 
  jet_p.clear(); 
  jet_E.clear(); 
  jet_eta.clear(); 
  jet_phi.clear(); 
  jet_z.clear();
  jet_qperp.clear();
  jet_largest.clear(); 
  jet_lab_eta.clear(); 
  jet_lab_phi.clear(); 
  jet_nc.clear();
  jet_pidx.clear(); 
  jet_pdR.clear(); 

  tjet_pT.clear(); 
  tjet_p.clear(); 
  tjet_E.clear(); 
  tjet_eta.clear(); 
  tjet_phi.clear(); 
  tjet_z.clear();
  tjet_qperp.clear();
  tjet_largest.clear(); 
  tjet_lab_eta.clear(); 
  tjet_lab_phi.clear(); 
  tjet_nc.clear();
  tjet_pidx.clear(); 
  tjet_pdR.clear(); 
  tjet_Q.clear(); 

  tcjet_pT.clear(); 
  tcjet_p.clear(); 
  tcjet_E.clear(); 
  tcjet_eta.clear(); 
  tcjet_phi.clear(); 
  tcjet_z.clear();
  tcjet_qperp.clear();
  tcjet_largest.clear(); 
  tcjet_lab_eta.clear(); 
  tcjet_lab_phi.clear(); 
  tcjet_nc.clear();
  tcjet_pidx.clear(); 
  tcjet_pdR.clear(); 
  tcjet_Q.clear(); 
  tcjet_cf.clear(); 
  tcjet_neut_p.clear(); 
  tcjet_chgd_p.clear(); 
  tcjet_em_p.clear(); 

  pjet_pT.clear(); 
  pjet_p.clear(); 
  pjet_E.clear(); 
  pjet_eta.clear(); 
  pjet_phi.clear(); 
  pjet_z.clear();
  pjet_qperp.clear();
  pjet_largest.clear(); 
  pjet_nc.clear();
  pjet_Q.clear(); 
  pjet_cf.clear(); 
  pjet_tcidx.clear(); 
  pjet_tcdR.clear(); 
  pjet_neut_p.clear(); 
  pjet_chgd_p.clear(); 
  pjet_em_p.clear(); 
  
  tfpjet_pT.clear(); 
  tfpjet_p.clear(); 
  tfpjet_E.clear(); 
  tfpjet_eta.clear(); 
  tfpjet_phi.clear(); 
  tfpjet_z.clear();
  tfpjet_qperp.clear();
  tfpjet_largest.clear(); 
  tfpjet_nc.clear();
  tfpjet_Q.clear(); 
  tfpjet_cf.clear(); 
  tfpjet_neut_p.clear(); 
  tfpjet_neut_p.clear(); 
  tfpjet_em_p.clear(); 

  if(electron!=NULL) {

    //std::cout << " electron found eta = " << electron->get_eta() 
    //	      << " p = " << electron->get_p() << " GeV" << std::endl; 

    double ElectronClusterEnergy = 0.0;
    double ElectronClusterDistance = 9999.0;
    std::string ECDetName = "";
    int ECIdx = -1; 
    RawCluster *ElectronClusterPointer = NULL; 

    for (SvtxTrack::ConstStateIter state_itr = electron->begin_states();
	 state_itr != electron->end_states(); state_itr++) {
      SvtxTrackState *temp = dynamic_cast<SvtxTrackState*>(state_itr->second);			  
      //std::cout << " State found at pathlength = " << temp->get_pathlength() 
      //		<< " " << temp->get_name() << std::endl; 
      if( (temp->get_pathlength()>0.0) && ((temp->get_name()=="BECAL") || (temp->get_name()=="EEMC")) ) {
	double clustE = 9999.0; 
	double clustdR = 9999.0; 
	int clustIdx = -1; 
	RawCluster *testCluster = getCluster(topNode, temp->get_name(), temp->get_eta(), temp->get_phi(), clustE, clustIdx, clustdR); 
	if((clustE>0.0)&&(clustdR<ElectronClusterDistance)) {
	  ElectronClusterEnergy = clustE;
	  ElectronClusterDistance = clustdR; 
	  ECDetName = temp->get_name(); 
	  ECIdx = clustIdx; 
	  ElectronClusterPointer = testCluster; 
	}
      }

    }

    meas_E_p = ElectronClusterEnergy/electron->get_p(); 
    electron_eta = electron->get_eta();
    electron_cluster_dR = ElectronClusterDistance; 

    // Calculate the event kinematics
	  
    // Ideally should be taken from event record? 
    double crossing_angle = 0.025; 

    TLorentzVector e_initial(0.0,0.0,-e_p_initial,e_p_initial); 
    TLorentzVector e_final(electron->get_px(), electron->get_py(), electron->get_pz(),electron->get_p()); 
    TLorentzVector virtual_photon = (e_initial - e_final); 
    measQ2 = -virtual_photon.Mag2(); 

    TLorentzVector p_initial(-p_p_initial*sin(crossing_angle), 0.0, p_p_initial*cos(crossing_angle), 
			     sqrt(pow(0.938,2) + pow(p_p_initial,2)));
    meas_x = measQ2/(2*virtual_photon*p_initial);

    // Set up the transformation (boost) to the Breit frame
    TVector3 P3 = p_initial.Vect(); 
    TVector3 q3 = virtual_photon.Vect(); 
    TVector3 boost = -(2.0*meas_x*P3 + q3)*(1.0/(2.0*meas_x*p_initial.E() + virtual_photon.E())); 
    TLorentzRotation breit = TLorentzRotation().Boost(boost); 
    TLorentzRotation breitInv = TLorentzRotation().Boost(-boost); 

    TLorentzVector p_initial_breit = (breit * p_initial); 
    TLorentzVector e_initial_breit = (breit * e_initial); 
    TLorentzVector e_final_breit = (breit * e_final); 
    TLorentzVector virtual_photon_breit = (breit * virtual_photon); 

    // DEBUG
    // std::cout << std::endl; 
    // std::cout << _hepmcp_Q2 << " " << measQ2 << " " << _hepmcp_x2 << " " << meas_x << std::endl; 
    // boost.Print(); 
    // std::cout << boost.Mag() << std::endl; 
    // virtual_photon_breit.Print(); 
    // p_initial_breit.Print(); 
    // TVector3 P3b = p_initial_breit.Vect(); 
    // TVector3 q3b = virtual_photon_breit.Vect();
    // // This should be zero!
    // (2.0*meas_x*P3b + q3b).Print();  
    // std::cout << std::endl; 

    // Now rotate so the virtual photon momentum is all along the negative z-axis

    TVector3 zhat =  -virtual_photon_breit.Vect().Unit(); 
    TVector3 yhat = (e_initial_breit.Vect().Cross(e_final_breit.Vect())).Unit(); 
    TVector3 xhat = yhat.Cross(zhat); 

    // std::cout << zhat.Dot(yhat) << " " << zhat.Dot(xhat) << " " << yhat.Dot(xhat) << std::endl; 
    // zhat.Print(); 
    // yhat.Print(); 
    // xhat.Print(); 

    TRotation breitRot; 
    breitRot.RotateAxes(xhat,yhat,zhat); 
    TRotation breitRotInv = breitRot.Inverse(); 

    p_initial_breit.Transform(breitRotInv); 
    virtual_photon_breit.Transform(breitRotInv); 
    e_initial_breit.Transform(breitRotInv); 
    e_final_breit.Transform(breitRotInv); 
	  
    breit_vphot_e = virtual_photon_breit.E(); 
    breit_initial_proton_e = p_initial_breit.E();
    breit_vphot_pz = virtual_photon_breit.Pz(); 
    breit_vphot_pt = virtual_photon_breit.Pt();
    breit_initial_proton_pz = p_initial_breit.Pz(); 
    breit_initial_proton_pt = p_initial_breit.Pt(); 

    // DEBUG
    // std::cout << " Post-Rotation: " << std::endl; 
    // virtual_photon_breit.Pt(); 
    // p_initial_breit.Print(); 
    // std::cout << std::endl;

    // Set up the Centauro algorithm for this event (in Breit frame)

    fastjet::JetDefinition* jetdef = 
      new fastjet::JetDefinition(new fastjet::contrib::CentauroPlugin(0.8));

    // Collect the calorimeter jets

    std::vector<fastjet::PseudoJet> pseudojets;
    std::vector<fastjet::PseudoJet> tcpseudojets;

    bool use_cemc = true; 
    bool use_ihcal = true; 
    bool use_ohcal = true; 
    bool use_femc = true; 
    bool use_fhcal = true; 

    // tower jets or cluster jets
    bool tower_jets = true; 

    if(tower_jets){

      // Tower Jets
      if(use_cemc) FillTowerPseudoJets(topNode, "BECAL", pseudojets, breit, breitRotInv, 1.0, ECDetName, ElectronClusterPointer); 
      if(use_ihcal) FillTowerPseudoJets(topNode, "HCALIN", pseudojets, breit, breitRotInv, 1.0); 
      if(use_ohcal) FillTowerPseudoJets(topNode, "HCALOUT", pseudojets, breit, breitRotInv, 1.0); 
      if(use_femc) FillTowerPseudoJets(topNode, "FEMC", pseudojets, breit, breitRotInv, 1.0); 
      if(use_fhcal) FillTowerPseudoJets(topNode, "LFHCAL", pseudojets, breit, breitRotInv, 1.0); 

    }
    else{

      // Cluster Jets
      if(use_cemc) FillClusterPseudoJets(topNode, "BECAL", pseudojets, breit, breitRotInv, 1.0, ECDetName, ECIdx); 
      if(use_ihcal) FillClusterPseudoJets(topNode, "HCALIN", pseudojets, breit, breitRotInv, 1.0); 
      if(use_ohcal) FillClusterPseudoJets(topNode, "HCALOUT", pseudojets, breit, breitRotInv, 1.0); 
      if(use_femc) FillClusterPseudoJets(topNode, "FEMC", pseudojets, breit, breitRotInv, 1.0); 
      if(use_fhcal) FillClusterPseudoJets(topNode, "LFHCAL", pseudojets, breit, breitRotInv, 1.0); 

    }

    fastjet::ClusterSequence jetFinder(pseudojets, *jetdef);
    // get jets sorted by energy
    std::vector<fastjet::PseudoJet> fastjets = sorted_by_E(jetFinder.inclusive_jets(0.0));
 
    for(unsigned int ijet=0; ijet<fastjets.size(); ijet++){
	    
      // tag highest E jet
      if(ijet==0)
	jet_largest.push_back(1); 
      else
	jet_largest.push_back(0);

      jet_pT.push_back(fastjets[ijet].perp()); 
      jet_p.push_back(sqrt(pow(fastjets[ijet].px(),2) + 
			   pow(fastjets[ijet].py(),2) + 
			   pow(fastjets[ijet].pz(),2)));
      jet_E.push_back(fastjets[ijet].E()); 
      jet_eta.push_back(fastjets[ijet].eta()); 
      jet_phi.push_back(fastjets[ijet].phi_02pi()); 

      TLorentzVector pjet(fastjets[ijet].px(),fastjets[ijet].py(),fastjets[ijet].pz(),fastjets[ijet].E());
      double zjet = (fastjets[ijet].E()-fastjets[ijet].pz())/sqrt(measQ2); 
      jet_z.push_back(zjet);
  
      jet_qperp.push_back(fastjets[ijet].perp()/zjet);

      // Transform back to lab frame
      pjet.Transform(breitRot); 
      pjet = (breitInv * pjet); 

      jet_lab_eta.push_back(pjet.Eta());  
      jet_lab_phi.push_back(pjet.Phi()); 

      jet_nc.push_back(fastjets[ijet].constituents().size());
	   
    }

    // Track Jets

    pseudojets.clear(); 

    for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
	 track_itr != _trackmap->end(); track_itr++) {

      SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

      // Skip scattered electron
      if ((temp->get_truth_track_id() - true_electron_headon->get_track_id()) == 0) continue; 

      double px = temp->get_px();
      double py = temp->get_py();
      double pz = temp->get_pz();

      // Transform to Breit frame
      TLorentzVector track(px,py,pz,temp->get_p()); 
      TLorentzVector breit_track = (breit*track); 

      breit_track.Transform(breitRotInv); 

      fastjet::PseudoJet pseudojet (breit_track.Px(),breit_track.Py(),breit_track.Pz(),breit_track.E());

      // Need to flag electrons/positrons to seperate EM fraction
      bool em_part = isThisAnElectron(temp); 

      pseudojet.set_user_index(EncodeUserIndex(temp->get_charge(),em_part)); 
      pseudojets.push_back(pseudojet);
      tcpseudojets.push_back(pseudojet);
  	   
    }

    fastjet::ClusterSequence tjetFinder(pseudojets, *jetdef);
    // get jets sorted by energy
    std::vector<fastjet::PseudoJet> tfastjets = sorted_by_E(tjetFinder.inclusive_jets(0.0));
 
    for(unsigned int ijet=0; ijet<tfastjets.size(); ijet++){
	    
      // tag highest E jet
      if(ijet==0)
	tjet_largest.push_back(1); 
      else
	tjet_largest.push_back(0);

      tjet_pT.push_back(tfastjets[ijet].perp()); 
      tjet_p.push_back(sqrt(pow(tfastjets[ijet].px(),2) + 
			    pow(tfastjets[ijet].py(),2) + 
			    pow(tfastjets[ijet].pz(),2)));
      tjet_E.push_back(tfastjets[ijet].E()); 
      tjet_eta.push_back(tfastjets[ijet].eta()); 
      tjet_phi.push_back(tfastjets[ijet].phi_02pi()); 

      TLorentzVector pjet(tfastjets[ijet].px(),tfastjets[ijet].py(),tfastjets[ijet].pz(),tfastjets[ijet].E());
      double zjet = (tfastjets[ijet].E()-tfastjets[ijet].pz())/sqrt(measQ2); 
      tjet_z.push_back(zjet);  

      tjet_qperp.push_back(tfastjets[ijet].perp()/zjet);

      // Transform back to lab frame
      pjet.Transform(breitRot); 
      pjet = (breitInv * pjet); 

      tjet_lab_eta.push_back(pjet.Eta());  
      tjet_lab_phi.push_back(pjet.Phi()); 

      tjet_nc.push_back(tfastjets[ijet].constituents().size());

      // Jet Charge 
      std::vector<fastjet::PseudoJet> tconstit = tfastjets[ijet].constituents();
      tjet_Q.push_back(JetCharge(&tconstit,sqrt(pow(tfastjets[ijet].px(),2) + 
						pow(tfastjets[ijet].py(),2) + 
						pow(tfastjets[ijet].pz(),2))));

    }

    // Add unassociated clusters to the track constituents and find jets again

    //if(use_cemc) FillClusterPseudoJets(topNode, "BECAL", tcpseudojets, breit, breitRotInv, 1.0, ECDetName, ECIdx, true); 
    //if(use_ihcal) FillClusterPseudoJets(topNode, "HCALIN", tcpseudojets, breit, breitRotInv, 1.0, "", -1, true); 
    //if(use_ohcal) FillClusterPseudoJets(topNode, "HCALOUT", tcpseudojets, breit, breitRotInv, 1.0, "", -1, true); 
    //if(use_femc) FillClusterPseudoJets(topNode, "FEMC", tcpseudojets, breit, breitRotInv, 1.0, "", -1, true); 
    //if(use_fhcal) FillClusterPseudoJets(topNode, "LFHCAL", tcpseudojets, breit, breitRotInv, 1.0, "", -1, true); 

    // Add Calo tracks to charged track constituents and find jets again
    
    BuildCaloTracks(topNode, "CENT", tcpseudojets, breit, breitRotInv, ECDetName, ECIdx); 
    BuildCaloTracks(topNode, "FWD", tcpseudojets, breit, breitRotInv, ECDetName, ECIdx); 

    fastjet::ClusterSequence tcjetFinder(tcpseudojets, *jetdef);
    // get jets sorted by energy
    std::vector<fastjet::PseudoJet> tcfastjets = sorted_by_E(tcjetFinder.inclusive_jets(0.0));
 
    for(unsigned int ijet=0; ijet<tcfastjets.size(); ijet++){
	    
      // tag highest E jet
      if(ijet==0)
	tcjet_largest.push_back(1); 
      else
	tcjet_largest.push_back(0);

      tcjet_pT.push_back(tcfastjets[ijet].perp()); 
      tcjet_p.push_back(sqrt(pow(tcfastjets[ijet].px(),2) + 
			     pow(tcfastjets[ijet].py(),2) + 
			     pow(tcfastjets[ijet].pz(),2)));
      tcjet_E.push_back(tcfastjets[ijet].E()); 
      tcjet_eta.push_back(tcfastjets[ijet].eta()); 
      tcjet_phi.push_back(tcfastjets[ijet].phi_02pi()); 

      TLorentzVector pjet(tcfastjets[ijet].px(),tcfastjets[ijet].py(),tcfastjets[ijet].pz(),tcfastjets[ijet].E());
      double zjet = (tcfastjets[ijet].E()-tcfastjets[ijet].pz())/sqrt(measQ2); 
      tcjet_z.push_back(zjet);  

      tcjet_qperp.push_back(tcfastjets[ijet].perp()/zjet);

      // Transform back to lab frame
      pjet.Transform(breitRot); 
      pjet = (breitInv * pjet); 

      tcjet_lab_eta.push_back(pjet.Eta());  
      tcjet_lab_phi.push_back(pjet.Phi()); 

      tcjet_nc.push_back(tcfastjets[ijet].constituents().size());
	 
      // Jet Charge
      std::vector<fastjet::PseudoJet> tconstit = tcfastjets[ijet].constituents();
	    
      double jptot = sqrt(pow(tcfastjets[ijet].px(),2) + 
			  pow(tcfastjets[ijet].py(),2) + 
			  pow(tcfastjets[ijet].pz(),2)); 

      tcjet_Q.push_back(JetCharge(&tconstit,jptot));

      tcjet_cf.push_back(JetChargedFraction(&tconstit,jptot));

      tcjet_neut_p.push_back(JetNeutralMomentum(&tconstit)); 
      tcjet_chgd_p.push_back(JetChargedMomentum(&tconstit)); 
      tcjet_em_p.push_back(JetEMMomentum(&tconstit)); 

    }

    // Primary Jets

    std::vector<fastjet::PseudoJet> pfastjets; 
    std::vector<fastjet::PseudoJet> tfpfastjets; 

    // Primary Jets (in the experiment Breit frame)
    GetPrimaryJets(topNode, jetdef, &pfastjets, breit, breitRotInv, p_initial_breit, virtual_photon_breit, false); 

    // Primary Jets in "true" Breit frame
    // Set up the transformation (boost) from the lab to the "true" Breit frame
    // Transformation based on truth kinematics
 
    // Get the final electron in the lab frame
    CLHEP::HepLorentzVector efp(true_electron_headon->get_px(),true_electron_headon->get_py(),
				true_electron_headon->get_pz(),true_electron_headon->get_e());
    efp = EventToLab * efp;  
    TLorentzVector true_e_final(efp.px(), efp.py(), efp.pz(), efp.e()); 
    TLorentzVector true_virtual_photon = (true_e_initial - true_e_final); 

    // Recalculate the kinematics 
    // (clean up round-off errors and keep everything self-consistent)
    double _tf_Q2 = -true_virtual_photon.Mag2(); 
    double _tf_x2 = _tf_Q2/(2*true_virtual_photon*true_p_initial);

    TVector3 P3t = true_p_initial.Vect(); 
    TVector3 q3t = true_virtual_photon.Vect(); 
    TVector3 boost_tf = -(2.0*_tf_x2*P3t + q3t)*(1.0/(2.0*_tf_x2*true_p_initial.E() + true_virtual_photon.E())); 
    TLorentzRotation breit_tf = TLorentzRotation().Boost(boost_tf); 

    TLorentzVector p_initial_breit_tf = (breit_tf * true_p_initial); 
    TLorentzVector e_initial_breit_tf = (breit_tf * true_e_initial); 
    TLorentzVector e_final_breit_tf = (breit_tf * true_e_final); 
    TLorentzVector virtual_photon_breit_tf = (breit_tf * true_virtual_photon); 

    // DEBUG
    // std::cout << std::endl; 
    // std::cout << _hepmcp_Q2 << " " << measQ2 << " " << _hepmcp_x2 << " " << meas_x << " " 
    // 	    << _tf_Q2 << " " << _tf_x2 << std::endl; 
    // boost_tf.Print(); 
    // std::cout << boost_tf.Mag() << std::endl; 
    // virtual_photon_breit_tf.Print(); 
    // p_initial_breit_tf.Print(); 
    // TVector3 P3bt = p_initial_breit_tf.Vect(); 
    // TVector3 q3bt = virtual_photon_breit_tf.Vect();
    // // This should be zero!
    // (2.0*_tf_x2*P3bt + q3bt).Print();  
    // std::cout << std::endl; 

    // Now rotate so the virtual photon momentum is all along the negative z-axis
    TVector3 zhat_tf =  -virtual_photon_breit_tf.Vect().Unit(); 
    TVector3 yhat_tf = (e_initial_breit_tf.Vect().Cross(e_final_breit_tf.Vect())).Unit(); 
    TVector3 xhat_tf = yhat_tf.Cross(zhat_tf); 

    // std::cout << zhat_tf.Dot(yhat_tf) << " " << zhat.Dot(xhat_tf) << " " << yhat.Dot(xhat_tf) << std::endl; 
    // zhat_tf.Print(); 
    // yhat_tf.Print(); 
    // xhat_tf.Print(); 

    TRotation breitRot_tf; 
    breitRot_tf.RotateAxes(xhat_tf,yhat_tf,zhat_tf); 
    TRotation breitRotInv_tf = breitRot_tf.Inverse(); 

    p_initial_breit_tf.Transform(breitRotInv_tf); 
    virtual_photon_breit_tf.Transform(breitRotInv_tf); 
    e_initial_breit_tf.Transform(breitRotInv_tf); 
    e_final_breit_tf.Transform(breitRotInv_tf); 

    // DEBUG
    // std::cout << " Post-Rotation: " << std::endl; 
    // virtual_photon_breit_tf.Print(); 
    // p_initial_breit_tf.Print(); 
    // std::cout << std::endl;

    GetPrimaryJets(topNode, jetdef, &tfpfastjets, breit_tf, breitRotInv_tf, p_initial_breit_tf, virtual_photon_breit_tf, true); 

    delete jetdef;

    // Match the reco jets to the primary jets
	  
    for(unsigned int ijet=0; ijet<fastjets.size(); ijet++){

      int closest_idx = -1; 
      double closest_dR = 9999.0; 

      for(unsigned int pjet=0; pjet<pfastjets.size(); pjet++){

	double dR = sqrt(pow(fastjets[ijet].eta() - pfastjets[pjet].eta(),2) + 
			 pow(DeltaPhi(fastjets[ijet].phi(),pfastjets[pjet].phi()),2) ); 

	if(dR<closest_dR){
	  closest_idx = pjet; 
	  closest_dR = dR; 
	}

      }

      jet_pidx.push_back(closest_idx); 
      jet_pdR.push_back(closest_dR); 

    }

    for(unsigned int ijet=0; ijet<tfastjets.size(); ijet++){

      int closest_idx = -1; 
      double closest_dR = 9999.0; 

      for(unsigned int pjet=0; pjet<pfastjets.size(); pjet++){

	//double dR = sqrt(pow(tfastjets[ijet].eta() - pfastjets[pjet].eta(),2) + 
	//	      pow(DeltaPhi(tfastjets[ijet].phi(),pfastjets[pjet].phi()),2) ); 

	TVector3 v1(tfastjets[ijet].px(),tfastjets[ijet].py(),tfastjets[ijet].pz()); 
	TVector3 v2(pfastjets[pjet].px(),pfastjets[pjet].py(),pfastjets[pjet].pz());

	double tc_p = sqrt(pow(tfastjets[ijet].px(),2) + 
			   pow(tfastjets[ijet].py(),2) +
			   pow(tfastjets[ijet].pz(),2));

	double p_p = sqrt(pow(pfastjets[pjet].px(),2) + 
			  pow(pfastjets[pjet].py(),2) + 
			  pow(pfastjets[pjet].pz(),2)); 

	double dR = acos(v1.Dot(v2)/(tc_p*p_p)); 

	if(dR<closest_dR){
	  closest_idx = pjet; 
	  closest_dR = dR; 
	}

      }

      tjet_pidx.push_back(closest_idx); 
      tjet_pdR.push_back(closest_dR); 

    }

    for(unsigned int ijet=0; ijet<tcfastjets.size(); ijet++){

      int closest_idx = -1; 
      double closest_dR = 9999.0; 

      for(unsigned int pjet=0; pjet<pfastjets.size(); pjet++){

	//double dR = sqrt(pow(tcfastjets[ijet].eta() - pfastjets[pjet].eta(),2) + 
	//	      pow(DeltaPhi(tcfastjets[ijet].phi(),pfastjets[pjet].phi()),2) ); 

	// use the angle between the two jets

	TVector3 v1(tcfastjets[ijet].px(),tcfastjets[ijet].py(),tcfastjets[ijet].pz()); 
	TVector3 v2(pfastjets[pjet].px(),pfastjets[pjet].py(),pfastjets[pjet].pz());

	double tc_p = sqrt(pow(tcfastjets[ijet].px(),2) + 
			   pow(tcfastjets[ijet].py(),2) +
			   pow(tcfastjets[ijet].pz(),2));

	double p_p = sqrt(pow(pfastjets[pjet].px(),2) + 
			  pow(pfastjets[pjet].py(),2) + 
			  pow(pfastjets[pjet].pz(),2)); 

	double dR = acos(v1.Dot(v2)/(tc_p*p_p)); 

	if(dR<closest_dR){
	  closest_idx = pjet; 
	  closest_dR = dR; 
	}

      }

      tcjet_pidx.push_back(closest_idx); 
      tcjet_pdR.push_back(closest_dR); 

    }

    // Now - loop over primary jets and tag the closest track+cluster jet

    for(unsigned int pjet=0; pjet<pfastjets.size(); pjet++){


      int closest_idx = -1; 
      double closest_dR = 9999.0; 

      for(unsigned int ijet=0; ijet<tcfastjets.size(); ijet++){

	TVector3 v1(tcfastjets[ijet].px(),tcfastjets[ijet].py(),tcfastjets[ijet].pz()); 
	TVector3 v2(pfastjets[pjet].px(),pfastjets[pjet].py(),pfastjets[pjet].pz());

	double tc_p = sqrt(pow(tcfastjets[ijet].px(),2) + 
			   pow(tcfastjets[ijet].py(),2) +
			   pow(tcfastjets[ijet].pz(),2));

	double p_p = sqrt(pow(pfastjets[pjet].px(),2) + 
			  pow(pfastjets[pjet].py(),2) + 
			  pow(pfastjets[pjet].pz(),2)); 

	double dR = acos(v1.Dot(v2)/(tc_p*p_p)); 

	if(dR<closest_dR){
	  closest_idx = pjet; 
	  closest_dR = dR; 
	}

      }

      pjet_tcidx.push_back(closest_idx); 
      pjet_tcdR.push_back(closest_dR); 

    }

  }
  else{
    LogWarning("reconstructed electron not found!"); 

    // DEBUG
    /*
    for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
	 track_itr != _trackmap->end(); track_itr++) {

      SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

      std::cout << temp->get_truth_track_id() << " " << true_electron_headon->get_track_id() << std::endl; 

    }
    */

  }

  _eval_tree_event->Fill();	

  return;
}

void CentauroJets::FillTowerPseudoJets( PHCompositeNode *topNode, std::string detName, 
					std::vector<fastjet::PseudoJet> &pseudojets, 
					TLorentzRotation &breit, TRotation &breitRot, 
					const double scale_factor, std::string ecDet, RawCluster *rcluster){

  string towernodename = "TOWER_CALIB_" + detName;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernodename.c_str());
  if (!towers)
    {
      std::cout << PHWHERE << ": Could not find node " << towernodename.c_str() << std::endl;
      return;
    }
  string towergeomnodename = "TOWERGEOM_" + detName;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename.c_str());
  if (! towergeom)
    {
      cout << PHWHERE << ": Could not find node " << towergeomnodename.c_str() << endl;
      return;
    }

  RawTowerContainer::ConstRange begin_end  = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;
  for (; itr != begin_end.second; ++itr){

    RawTowerDefs::keytype towerid = itr->first;

    // Make sure this tower isn't in the electron cluster
    
    if((ecDet==detName) && rcluster){

      bool inList = false; 

      RawCluster::TowerConstRange begin_end  = rcluster->get_towers();
      RawCluster::TowerConstIterator itr = begin_end.first;
      for (; itr != begin_end.second; ++itr){

	RawTowerDefs::keytype rtowerid = itr->first;
	
	if(rtowerid==towerid) {
	  inList = true; 
	  break; 
	}
    
      }

      if(inList) continue; 

    }

    RawTowerGeom *tgeo = towergeom->get_tower_geometry(towerid); 
    RawTower *ctower = towers->getTower(towerid);

    // Tower pedestal energy cut: 
    if(ctower->get_energy()<=TOWER_E_CUT) continue; 

    double x = tgeo->get_center_x(); 
    double y = tgeo->get_center_y(); 
    double z = tgeo->get_center_z(); 

    double eta = getEta(sqrt(pow(x,2)+pow(y,2)),z-vtx_z); 
    double phi = XYtoPhi(x,y); 
 
    double pt = scale_factor*ctower->get_energy() / cosh(eta);
    double px = pt * cos(phi);
    double py = pt * sin(phi);
    double pz = pt * sinh(eta);

    // Transform to Breit frame
    TLorentzVector tower(px,py,pz,scale_factor*ctower->get_energy()); 
    TLorentzVector breit_tower = (breit*tower); 
    breit_tower.Transform(breitRot); 

    fastjet::PseudoJet pseudojet (breit_tower.Px(),breit_tower.Py(),breit_tower.Pz(),breit_tower.E()); 
    pseudojet.set_user_index(0); 
    pseudojets.push_back(pseudojet);
    
  }

  return; 

}

bool CentauroJets::VetoClusterWithTrack(double eta, double phi, std::string detName){

  // Does this cluster have a track pointing to it? 
      
  double minDist = 9999.0; 

  for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
       track_itr != _trackmap->end(); track_itr++) {

    SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

    for (SvtxTrack::ConstStateIter state_itr = temp->begin_states();
	 state_itr != temp->end_states(); state_itr++) {

      SvtxTrackState *tstate = dynamic_cast<SvtxTrackState*>(state_itr->second);			  

      if( (tstate->get_pathlength()>0.0) && (tstate->get_name()==detName) ) {

	double deta = eta -  tstate->get_eta(); 
	double dPhi = DeltaPhi(phi, tstate->get_phi()); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) ); 
	if(dist<minDist){
	  minDist = dist; 
	}

      }

    }

  }

  //_h_track_cluster_match->Fill(minDist);

  double cutDist = 0.15; 

  if(detName=="BECAL") {
    //_h_track_cluster_match_becal->Fill(minDist);
    cutDist = BECAL_CLUST_TRACKMATCH; 
  }

  if(detName=="HCALIN") {
    //_h_track_cluster_match_ihcal->Fill(minDist);
    cutDist = IHCAL_CLUST_TRACKMATCH; 
  }

  if(detName=="HCALOUT") {
    //_h_track_cluster_match_ohcal->Fill(minDist);
    cutDist = OHCAL_CLUST_TRACKMATCH; 
  }

  if(detName=="FEMC") {
    //_h_track_cluster_match_femc->Fill(minDist);
    cutDist = FEMC_CLUST_TRACKMATCH; 
  }

  if(detName=="LFHCAL") {
    //_h_track_cluster_match_lfhcal->Fill(minDist);
    cutDist = LFHCAL_CLUST_TRACKMATCH; 
  }

  // Veto this cluster if a track points to it. 
  if(minDist<cutDist)
    return true;
  else
    return false; 

}

void CentauroJets::BuildCaloTracks(PHCompositeNode *topNode, std::string type, 
				   std::vector<fastjet::PseudoJet> &pseudojets, 
				   TLorentzRotation &breit, TRotation &breitRot, 
				   std::string ecDet, int ecIdx ){

  std::string fwd_calos[3] = {"FEMC","LFHCAL",""}; 
  std::string cent_calos[3] = {"BECAL","HCALIN","HCALOUT"}; 
  std::string detName[3] = {"","",""}; 

  // Store the starting index so we don't double-count 
  // when filling the diagnostic output
  int startIdx = pseudojets.size(); 

  // Collect the required cluster lists

  RawClusterContainer *clusterList[3] = {NULL, NULL, NULL};

  for(int i=0; i<3; i++){

    if(type=="CENT") 
      detName[i] = cent_calos[i]; 
    else if(type=="FWD") 
      detName[i] = fwd_calos[i]; 
    else{
      cout << " Unknown calo track type = " << type << endl; 
      return; 
    }
    
    if(detName[i]=="") continue; 
    
    string clusternodename = "CLUSTER_" + detName[i];
    clusterList[i] = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
    if (!clusterList[i]) {
      cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
      return;
    }    

  }

  // Create the lists to mark a cluster as used
  
  std::vector<bool> cused0(clusterList[0]->size(),false); 

  int size1; 
  if(clusterList[1])
    size1 = clusterList[1]->size(); 
  else
    size1 = 1; 
  std::vector<bool> cused1(size1,false);
  
  int size2; 
  if(clusterList[2])
    size2 = clusterList[2]->size(); 
  else
    size2 = 1; 
  std::vector<bool> cused2(size2,false); 
    
  // Create and fill the track matching lists

  std::vector<bool> tmatched0(clusterList[0]->size(),false); 
  std::vector<bool> tmatched1(size1,false);
  std::vector<bool> tmatched2(size2,false); 

  for (unsigned int i = 0; i < 3; i++) {

    if(!clusterList[i]) continue; 

    for (unsigned int k = 0; k < clusterList[i]->size(); k++) {

      RawCluster *rcluster = clusterList[i]->getCluster(k);

      // eliminate noise clusters
      if(rcluster->get_energy()<CLUSTER_E_CUTOFF) continue; 

      double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
      double phi = rcluster->get_phi(); 

      if(i==0) {tmatched0[k] = VetoClusterWithTrack(eta, phi, detName[i]); cused0[k] = tmatched0[k];}
      if(i==1) {tmatched1[k] = VetoClusterWithTrack(eta, phi, detName[i]); cused1[k] = tmatched1[k];}
      if(i==2) {tmatched2[k] = VetoClusterWithTrack(eta, phi, detName[i]); cused2[k] = tmatched2[k];}

    }

    if(detName[i]=="BECAL") _h_nclusters_becal->Fill(clusterList[i]->size()); 
    if(detName[i]=="HCALIN") _h_nclusters_ihcal->Fill(clusterList[i]->size()); 
    if(detName[i]=="HCALOUT") _h_nclusters_ohcal->Fill(clusterList[i]->size()); 
    if(detName[i]=="FEMC") _h_nclusters_femc->Fill(clusterList[i]->size()); 
    if(detName[i]=="LFHCAL") _h_nclusters_lfhcal->Fill(clusterList[i]->size()); 

  }

  // Mark the electron cluster as used
  // (the EMCal must always be the first in the list for this to work)
    
  if( ecDet==detName[0] ) cused0[(int)ecIdx] = true; 

  // First, seed with the EMCal and look for backing energy in the HCALs

  for (unsigned int k = 0; k < clusterList[0]->size(); k++) {

    if(cused0[k]) continue; 

    RawCluster *rcluster0 = clusterList[0]->getCluster(k);

    // eliminate noise clusters
    if(rcluster0->get_energy()<CLUSTER_E_CUTOFF) continue; 

    double eta = getEta(rcluster0->get_r(),rcluster0->get_z()-vtx_z);
    double phi = rcluster0->get_phi(); 
 
    double pt = rcluster0->get_energy() / cosh(eta);
    double px = pt * cos(phi);
    double py = pt * sin(phi);
    double pz = pt * sinh(eta);

    // Create the cluster
    TLorentzVector cluster(px,py,pz,rcluster0->get_energy()); 
    
    cused0[k] = true; 

    // mark as a photon candidate - will be set to false if we 
    // find hadronic energy further down

    bool photon_candidate = true; 

    // Look for the closest match in the next detector 

    // For the HCALs we need to collect all the clusters
    // that are within the window (collect split clusters)
    
    if(detName[1]!=""){

      TLorentzVector cluster1(0.0,0.0,0.0,0.0);

      for (unsigned int j = 0; j < clusterList[1]->size(); j++) {

	if(cused1[j]) continue; 

	RawCluster *rcluster1 = clusterList[1]->getCluster(j);

	// eliminate noise clusters
	if(rcluster1->get_energy()<CLUSTER_E_CUTOFF) continue; 

	double eta1 = getEta(rcluster1->get_r(),rcluster1->get_z()-vtx_z);
	double phi1 = rcluster1->get_phi(); 

	double deta = eta -  eta1; 
	double dPhi = DeltaPhi(phi, phi1); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) );
     
	double mDist = 0.0;
	double scale = 1.0; 
	if(type=="CENT"){
	  _h_becal_ihcal_match->Fill(dist);
	  _h_becal_ihcal_match_eta->Fill(deta); 
	  _h_becal_ihcal_match_phi->Fill(dPhi); 
	  mDist = IHCAL_CLUST_TRACKMATCH;
	  scale = BARREL_HCAL_NEUT_SCALE; 
	}
	else if(type=="FWD"){
	  _h_femc_lfhcal_match->Fill(dist);
	  _h_femc_lfhcal_match_eta->Fill(deta);
	  _h_femc_lfhcal_match_phi->Fill(dPhi);
	  mDist = LFHCAL_CLUST_TRACKMATCH;
	  scale = FWD_HCAL_NEUT_SCALE; 
	}

	if(dist<mDist){ 

	  double pt1 = rcluster1->get_energy() / cosh(eta1);
	  double px1 = pt1 * cos(phi1);
	  double py1 = pt1 * sin(phi1);
	  double pz1 = pt1 * sinh(eta1);

	  TLorentzVector clusterAdd(px1,py1,pz1,rcluster1->get_energy()); 
	  cluster1 += (clusterAdd*scale); 

	  photon_candidate = false; 

	  cused1[j] = true; 

	}

      }

      // Add what we found to the existing cluster
      cluster += cluster1;

    }

    // And the last detector (if it exists)

    if(detName[2]!="") {

      TLorentzVector cluster2(0.0,0.0,0.0,0.0); 

      for (unsigned int j = 0; j < clusterList[2]->size(); j++) {

	if(cused2[j]) continue; 

	RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	// eliminate noise clusters
	if(rcluster2->get_energy()<CLUSTER_E_CUTOFF) continue; 

	double eta2 = getEta(rcluster2->get_r(),rcluster2->get_z()-vtx_z);
	double phi2 = rcluster2->get_phi(); 

	double deta = eta -  eta2; 
	double dPhi = DeltaPhi(phi, phi2); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) ); 
	
	_h_becal_ohcal_match->Fill(dist); 
	_h_becal_ohcal_match_eta->Fill(deta); 
	_h_becal_ohcal_match_phi->Fill(dPhi); 

	if(dist<OHCAL_CLUST_TRACKMATCH){ 

	  double pt2 = rcluster2->get_energy() / cosh(eta2);
	  double px2 = pt2 * cos(phi2);
	  double py2 = pt2 * sin(phi2);
	  double pz2 = pt2 * sinh(eta2);

	  TLorentzVector clusterAdd(px2,py2,pz2,rcluster2->get_energy()); 
	  cluster2 += (clusterAdd*BARREL_HCAL_NEUT_SCALE); 

	  photon_candidate = false; 

	  cused2[j] = true; 

	}

      }

      // Add what we found to the existing cluster
      cluster += cluster2;

    }

    // Finally - take the combined cluster and add it to the constituents
    // Transform to Breit frame
    TLorentzVector breit_cluster = (breit*cluster); 
    breit_cluster.Transform(breitRot); 

    fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
    pseudojet.set_user_index(EncodeUserIndex(0,photon_candidate)); 
    pseudojets.push_back(pseudojet);

  }

  // Next, seed with the HCAL and follow up with the second HCAL segment
  // necessary to capture hadronic showers w/o cluster in EMCal
  // we need a double-loop to aggregate the split clusters

  if(detName[1]!=""){

    for (unsigned int k = 0; k < clusterList[1]->size(); k++) {

      if(cused1[k]) continue;

      RawCluster *rcluster = clusterList[1]->getCluster(k);

      // eliminate noise clusters
      if(rcluster->get_energy()<CLUSTER_E_CUTOFF) continue; 

      double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
      double phi = rcluster->get_phi(); 

      double pt = rcluster->get_energy() / cosh(eta);
      double px = pt * cos(phi);
      double py = pt * sin(phi);
      double pz = pt * sinh(eta);

      double scale = 1.0; 
      if(type=="CENT"){
	scale = BARREL_HCAL_NEUT_SCALE; 
      }
      else if(type=="FWD"){
	scale = FWD_HCAL_NEUT_SCALE; 
      }

      // Create the cluster
      TLorentzVector cluster(px,py,pz,rcluster->get_energy()); 
      cluster *= scale; 

      cused1[k] = true; 

      // Aggregate split clusters

      for (unsigned int j = 0; j < clusterList[1]->size(); j++) {
      
	if(cused1[j]) continue;

	RawCluster *rcluster2 = clusterList[1]->getCluster(j);

	// eliminate noise clusters
        if(rcluster2->get_energy()<CLUSTER_E_CUTOFF) continue; 

	double eta2 = getEta(rcluster2->get_r(),rcluster2->get_z()-vtx_z);
	double phi2 = rcluster2->get_phi(); 

	double deta = eta -  eta2; 
	double dPhi = DeltaPhi(phi, phi2); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) );

	double mDist = 0.0; 
	if(type=="CENT"){
	  mDist = IHCAL_CLUST_TRACKMATCH;
	}
	else if(type=="FWD"){
	  mDist = LFHCAL_CLUST_TRACKMATCH; 
	}

	if(dist<mDist){ 

	  double pt2 = rcluster2->get_energy() / cosh(eta2);
	  double px2 = pt2 * cos(phi2);
	  double py2 = pt2 * sin(phi2);
	  double pz2 = pt2 * sinh(eta2);

	  TLorentzVector clusterAdd(px2,py2,pz2,rcluster2->get_energy()); 
	  cluster += (clusterAdd*scale); 

	  cused1[j] = true; 

	}

      }

      // And the last detector (if it exists)

      if(detName[2]!="") {

	TLorentzVector cluster2(0.0,0.0,0.0,0.0); 

 	for (unsigned int j = 0; j < clusterList[2]->size(); j++) {

	  if(cused2[j]) continue; 

	  RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	  // eliminate noise clusters
	  if(rcluster2->get_energy()<CLUSTER_E_CUTOFF) continue; 

	  double eta2 = getEta(rcluster2->get_r(),rcluster2->get_z()-vtx_z);
	  double phi2 = rcluster2->get_phi(); 

	  double deta = eta -  eta2; 
	  double dPhi = DeltaPhi(phi, phi2); 

	  double dist = sqrt( pow(deta,2) + pow(dPhi,2) );

	  if(type=="CENT"){
	    _h_ihcal_ohcal_match->Fill(dist); 
	    _h_ihcal_ohcal_match_eta->Fill(deta); 
	    _h_ihcal_ohcal_match_phi->Fill(dPhi); 
	  }

	  if(dist<OHCAL_CLUST_TRACKMATCH){ 

	    double pt2 = rcluster2->get_energy() / cosh(eta2);
	    double px2 = pt2 * cos(phi2);
	    double py2 = pt2 * sin(phi2);
	    double pz2 = pt2 * sinh(eta2);

	    TLorentzVector clusterAdd(px2,py2,pz2,rcluster2->get_energy()); 
	    cluster2 +=(clusterAdd*BARREL_HCAL_NEUT_SCALE); 

	    cused2[j] = true; 

	  }

	}

	// Add what we found to the existing cluster
	cluster += cluster2;

      }

      // Finally - take the combined cluster and add it to the constituents
      // Transform to Breit frame
      TLorentzVector breit_cluster = (breit*cluster); 
      breit_cluster.Transform(breitRot); 

      fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
      pseudojet.set_user_index(0); 
      pseudojets.push_back(pseudojet);

    }

  }

  // Finally - the last HCAL by itself (if it exists)

  if(detName[2]!="") {

    for (unsigned int k = 0; k < clusterList[2]->size(); k++) {

      if(cused2[k]) continue; 

      RawCluster *rcluster1 = clusterList[2]->getCluster(k);

      // eliminate noise clusters
      if(rcluster1->get_energy()<CLUSTER_E_CUTOFF) continue; 

      double eta1 = getEta(rcluster1->get_r(),rcluster1->get_z()-vtx_z);
      double phi1 = rcluster1->get_phi(); 

      double pt1 = rcluster1->get_energy() / cosh(eta1);
      double px1 = pt1 * cos(phi1);
      double py1 = pt1 * sin(phi1);
      double pz1 = pt1 * sinh(eta1);

      // Create the cluster
      TLorentzVector cluster(px1,py1,pz1,rcluster1->get_energy()); 

      cused2[k] = true; 

      // Aggregate the split clusters

      for (unsigned int j = 0; j < clusterList[2]->size(); j++) {
      
	if(cused2[j]) continue;

	RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	// eliminate noise clusters
        if(rcluster2->get_energy()<CLUSTER_E_CUTOFF) continue; 

	double eta2 = getEta(rcluster2->get_r(),rcluster2->get_z()-vtx_z);
	double phi2 = rcluster2->get_phi(); 

	double deta = eta1 -  eta2; 
	double dPhi = DeltaPhi(phi1, phi2); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) );

	double mDist = 0.0; 
	if(type=="CENT"){
	  mDist = OHCAL_CLUST_TRACKMATCH;
	}

	if(dist<mDist){ 

	  double pt2 = rcluster2->get_energy() / cosh(eta2);
	  double px2 = pt2 * cos(phi2);
	  double py2 = pt2 * sin(phi2);
	  double pz2 = pt2 * sinh(eta2);

	  TLorentzVector clusterAdd(px2,py2,pz2,rcluster2->get_energy()); 
	  cluster += (clusterAdd*BARREL_HCAL_NEUT_SCALE); 

	  cused2[j] = true; 

	}

      }

      // Transform to Breit frame
      TLorentzVector breit_cluster = (breit*cluster); 
      breit_cluster.Transform(breitRot); 

      fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
      pseudojet.set_user_index(0); 
      pseudojets.push_back(pseudojet);

    }

  }

  // Diagnostics - connect the calo tracks to neutral primaries

  cat_pid.clear();  
  cat_p_true.clear(); 
  cat_eta_meas.clear(); 
  cat_eta_true.clear(); 
  cat_match.clear(); 
  cat_e_tot.clear(); 

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  for(unsigned int i=startIdx; i<pseudojets.size(); i++){

    if(pseudojets[i].user_index()==0){

      TVector3 ctrack(pseudojets[i].px(),pseudojets[i].py(),pseudojets[i].pz());
      
      double minDist = 9999.0; 
      int pid = -9999; 
      double prim_p = 9999.0; 
      double prim_Eta = 9999.0; 

      // PRIMARIES ONLY
      PHG4TruthInfoContainer::ConstRange range =
	_truth_container->GetPrimaryParticleRange();

      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
	   truth_itr != range.second; ++truth_itr) {

	PHG4Particle* g4particle = truth_itr->second;
	if(!g4particle) {
	  LogDebug("");
	  continue;
	}

	G4ParticleDefinition* particle = particleTable->FindParticle(g4particle->get_name());
	int charge = -9999; 
	if(particle) 
	  charge = particle->GetPDGCharge();
	else 
	  continue; 

	if(charge!=0) continue; 

	// NOTE: stored HepMC kinematics are in event generator (head-on) frame!
	// Transform to the lab frame
	CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),g4particle->get_pz(),g4particle->get_e());
	efp = EventToLab * efp;  

	// Now transform to the Breit frame
	TLorentzVector partMom(efp.px(), efp.py(), efp.pz(), efp.e()); 
	TLorentzVector partMom_breit = (breit*partMom); 
	partMom_breit.Transform(breitRot); 

	double deta = ctrack.Eta() - partMom_breit.Eta(); 
	double dphi = DeltaPhi(ctrack.Phi(), partMom_breit.Phi()); 

	double dist = sqrt(pow(deta,2) + pow(dphi,2)); 
	if(dist<minDist){
	  minDist = dist; 
	  pid = g4particle->get_pid();
	  prim_p = partMom_breit.Vect().Mag(); 
	  prim_Eta = partMom_breit.Vect().Eta(); 
	}

      }

      if(minDist<9999.0){

	cat_pid.push_back(pid);
	cat_p_true.push_back(prim_p); 
        cat_e_tot.push_back(ctrack.Mag());
	cat_match.push_back(minDist); 
	cat_eta_meas.push_back(ctrack.Eta()); 
	cat_eta_true.push_back(prim_Eta); 

	if(type=="CENT"){
	  _h_calotrack_prim_match_cent->Fill(minDist); 
	  _h_calotrack_pid_prim_match_cent->Fill(pid); 
	  _h_calotrack_cent_NES->Fill(ctrack.Mag()/prim_p); 
	  _h_calotrack_cent_NES_2D->Fill(prim_p,ctrack.Mag()/prim_p); 
	  if(pid==22) {
	    _h_calotrack_prim_match_cent_gamma->Fill(minDist); 
	    _h_calotrack_cent_gamma_NES_2D->Fill(prim_p,ctrack.Mag()/prim_p);
	  }
	  if(pid==2112) {
	    _h_calotrack_prim_match_cent_neutron->Fill(minDist); 
	    _h_calotrack_cent_neutron_NES_2D->Fill(prim_p,ctrack.Mag()/prim_p);
	  }
	}

	if(type=="FWD"){
	  _h_calotrack_prim_match_fwd->Fill(minDist); 
	  _h_calotrack_pid_prim_match_fwd->Fill(pid); 
	  _h_calotrack_fwd_NES->Fill(ctrack.Mag()/prim_p); 
	  _h_calotrack_fwd_NES_2D->Fill(prim_p,ctrack.Mag()/prim_p); 
	  if(pid==22) {
	    _h_calotrack_prim_match_fwd_gamma->Fill(minDist); 
	    _h_calotrack_fwd_gamma_NES_2D->Fill(prim_p,ctrack.Mag()/prim_p);
	  }
	  if(pid==2112) {
	    _h_calotrack_prim_match_fwd_neutron->Fill(minDist); 
	    _h_calotrack_fwd_neutron_NES_2D->Fill(prim_p,ctrack.Mag()/prim_p);
	  }
	}
      
      }

    }

  }

  if(type=="CENT") _eval_calo_tracks_cent->Fill(); 
  if(type=="FWD") _eval_calo_tracks_fwd->Fill(); 

  return; 

}

SvtxTrack *CentauroJets::AttachClusterToTrack(double eta, double phi, std::string detName){

  // Does this cluster have a track pointing to it? 
      
  double minDist = 9999.0; 
  SvtxTrack *closest = NULL; 

  for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
       track_itr != _trackmap->end(); track_itr++) {

    SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

    for (SvtxTrack::ConstStateIter state_itr = temp->begin_states();
	 state_itr != temp->end_states(); state_itr++) {

      SvtxTrackState *tstate = dynamic_cast<SvtxTrackState*>(state_itr->second);			  

      if( (tstate->get_pathlength()>0.0) && (tstate->get_name()==detName) ) {

	double deta = eta -  tstate->get_eta(); 
	double dPhi = DeltaPhi(phi, tstate->get_phi()); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) ); 
	if(dist<minDist){
	  minDist = dist; 
	  closest = temp; 
	}

      }

    }

  }

  _h_track_cluster_match->Fill(minDist);

  double cutDist = 0.15; 

  if(detName=="BECAL") {
    _h_track_cluster_match_becal->Fill(minDist);
    cutDist = BECAL_CLUST_TRACKMATCH; 
  }

  if(detName=="HCALIN") {
    _h_track_cluster_match_ihcal->Fill(minDist);
    cutDist = IHCAL_CLUST_TRACKMATCH; 
  }

  if(detName=="HCALOUT") {
    _h_track_cluster_match_ohcal->Fill(minDist);
    cutDist = OHCAL_CLUST_TRACKMATCH; 
  }

  if(detName=="FEMC") {
    _h_track_cluster_match_femc->Fill(minDist);
    cutDist = FEMC_CLUST_TRACKMATCH; 
  }

  if(detName=="LFHCAL") {
    _h_track_cluster_match_lfhcal->Fill(minDist);
    cutDist = LFHCAL_CLUST_TRACKMATCH; 
  }

  // Return the track the cluster is closest to 
  if(minDist<cutDist)
    return closest;
  else
    return NULL; 

}

void CentauroJets::BuildChargedCaloTracks(PHCompositeNode *topNode, std::string type){

  std::string fwd_calos[3] = {"FEMC","LFHCAL",""}; 
  std::string cent_calos[3] = {"BECAL","HCALIN","HCALOUT"}; 
  std::string detName[3] = {"","",""}; 

  // Collect the required cluster lists

  RawClusterContainer *clusterList[3] = {NULL, NULL, NULL};

  for(int i=0; i<3; i++){

    if(type=="CENT") 
      detName[i] = cent_calos[i]; 
    else if(type=="FWD") 
      detName[i] = fwd_calos[i]; 
    else{
      cout << " Unknown calo track type = " << type << endl; 
      return; 
    }
    
    if(detName[i]=="") continue; 
    
    string clusternodename = "CLUSTER_" + detName[i];
    clusterList[i] = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
    if (!clusterList[i]) {
      cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
      return;
    }    

  }

  // Create and fill the track matching lists

  int size1; 
  if(clusterList[1])
    size1 = clusterList[1]->size(); 
  else
    size1 = 1; 
  
  int size2; 
  if(clusterList[2])
    size2 = clusterList[2]->size(); 
  else
    size2 = 1; 
    
  std::vector<SvtxTrack *> tmatched0(clusterList[0]->size(),NULL); 
  std::vector<SvtxTrack *> tmatched1(size1,NULL);
  std::vector<SvtxTrack *> tmatched2(size2,NULL); 

  std::vector<bool> cused0(clusterList[0]->size(),false); 
  std::vector<bool> cused1(size1,false);
  std::vector<bool> cused2(size2,false); 

  for (unsigned int i = 0; i < 3; i++) {

    if(!clusterList[i]) continue; 

    for (unsigned int k = 0; k < clusterList[i]->size(); k++) {

      RawCluster *rcluster = clusterList[i]->getCluster(k);

      // eliminate noise clusters
      if(rcluster->get_energy()<CLUSTER_E_CUTOFF) continue; 

      double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
      double phi = rcluster->get_phi(); 

      if(i==0) {tmatched0[k] = AttachClusterToTrack(eta, phi, detName[i]);}
      if(i==1) {tmatched1[k] = AttachClusterToTrack(eta, phi, detName[i]);}
      if(i==2) {tmatched2[k] = AttachClusterToTrack(eta, phi, detName[i]);}

    }

  }

  // Clear the vectors

  ct_pid.clear(); 
  ct_p_meas.clear();  
  ct_p_true.clear(); 
  ct_eta_meas.clear(); 
  ct_eta_true.clear(); 
  ct_dist.clear(); 
  ct_e_bemc.clear(); 
  ct_e_ihcal.clear(); 
  ct_e_ohcal.clear(); 
  ct_e_femc.clear(); 
  ct_e_lfhcal.clear(); 
  ct_e_tot.clear(); 

  // Attach the clusters to the tracks and fill the tree
  // start seeded with the first calorimeter (EMCAL)

  for (unsigned int k = 0; k < clusterList[0]->size(); k++) {

    if(!tmatched0[k]) continue; 
    if(cused0[k]) continue; 

    RawCluster *rcluster = clusterList[0]->getCluster(k);

    // Get the truth information for this track

    bool found = false; 
    PHG4TruthInfoContainer::ConstRange range = _truth_container->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
	 truth_itr != range.second; ++truth_itr) {

      PHG4Particle* g4particle = truth_itr->second;
      if(!g4particle) {
	LogDebug("");
	continue;
      }

      if ((tmatched0[k]->get_truth_track_id() - g4particle->get_track_id()) == 0) {

	found = true;

	ct_pid.push_back(g4particle->get_pid()); 

	// Get the final particle in the lab frame
	CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),
				g4particle->get_pz(),g4particle->get_e());
	efp = EventToLab * efp;  
	TLorentzVector match_lf(efp.px(), efp.py(), efp.pz(), efp.e());

	ct_p_meas.push_back(tmatched0[k]->get_p()); 
	ct_p_true.push_back(match_lf.Vect().Mag()); 
	
	ct_eta_meas.push_back(tmatched0[k]->get_eta()); 
        ct_eta_true.push_back(match_lf.Vect().Eta()); 

	double deta = tmatched0[k]->get_eta() - match_lf.Vect().Eta(); 
	double dphi = DeltaPhi(tmatched0[k]->get_phi(),match_lf.Vect().Phi()); 

	ct_dist.push_back(sqrt(pow(dphi,2)+ pow(deta,2))); 

	cused0[k] = true; 

	break; 

      }
	
    }
    
    if(!found) continue; 

    double caloTot = rcluster->get_energy(); 
 
    if(type=="CENT"){
      ct_e_bemc.push_back(rcluster->get_energy()); 
    }
    else if(type=="FWD"){
      ct_e_femc.push_back(rcluster->get_energy()); 
    }

    double stage2_energy = 0.0; 

    if(detName[1]!=""){

      for (unsigned int j = 0; j < clusterList[1]->size(); j++) {
	
	if(!tmatched1[j]) continue;
	if(cused1[j]) continue; 

	// Look for same track match
	if(tmatched1[j]==tmatched0[k]){

	  RawCluster *rcluster1 = clusterList[1]->getCluster(j);

	  caloTot += rcluster1->get_energy(); 
	  stage2_energy = rcluster1->get_energy(); 

	  cused1[j] = true; 

	  //break; 

	}

      }

    }

    if(type=="CENT"){
      ct_e_ihcal.push_back(stage2_energy); 
    }
    else if(type=="FWD"){
      ct_e_lfhcal.push_back(stage2_energy); 
    }

    // Last calorimeter

    double stage3_energy = 0.0; 

    if(detName[2]!=""){

      for (unsigned int j = 0; j < clusterList[2]->size(); j++) {

	if(!tmatched2[j]) continue;
	if(cused2[j]) continue; 

	// Look for same track match
	if(tmatched2[j]==tmatched0[k]){

	  RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	  caloTot += rcluster2->get_energy();
	  stage3_energy += rcluster2->get_energy();
	  
	  cused2[j] = true; 

	  //break; 

	}

      }

    }

    if(type=="CENT"){
      ct_e_ohcal.push_back(stage3_energy); 
    }

    // record the combined energy
    ct_e_tot.push_back(caloTot); 
   
  }

  // Now look for tracks seeded by the first HCAL (no EMCAL cluster)

  if(detName[1]!=""){

    for (unsigned int k = 0; k < clusterList[1]->size(); k++) {

      if(!tmatched1[k]) continue;
      if(cused1[k]) continue; 

      RawCluster *rcluster = clusterList[1]->getCluster(k);

      // Get the truth information for this track

      bool found = false; 
      double stage1_energy = 0.0; 

      PHG4TruthInfoContainer::ConstRange range = _truth_container->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
	   truth_itr != range.second; ++truth_itr) {

	PHG4Particle* g4particle = truth_itr->second;
	if(!g4particle) {
	  LogDebug("");
	  continue;
	}

	if ((tmatched1[k]->get_truth_track_id() - g4particle->get_track_id()) == 0) {

	  if(!found){

	    ct_pid.push_back(g4particle->get_pid()); 

	    // Get the final particle in the lab frame
	    CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),
					g4particle->get_pz(),g4particle->get_e());
	    efp = EventToLab * efp;  
	    TLorentzVector match_lf(efp.px(), efp.py(), efp.pz(), efp.e());

	    ct_p_meas.push_back(tmatched1[k]->get_p()); 
	    ct_p_true.push_back(match_lf.Vect().Mag()); 
	
	    ct_eta_meas.push_back(tmatched1[k]->get_eta()); 
	    ct_eta_true.push_back(match_lf.Vect().Eta()); 

	    double deta = tmatched1[k]->get_eta() - match_lf.Vect().Eta(); 
	    double dphi = DeltaPhi(tmatched1[k]->get_phi(),match_lf.Vect().Phi());

	    ct_dist.push_back(sqrt(pow(dphi,2)+ pow(deta,2))); 

	  }

	  found = true;

	  cused1[k] = true; 

	  stage1_energy += rcluster->get_energy(); 

	  //break; 

	}
	
      }
    
      if(!found) continue; 

      double caloTot = stage1_energy; 
 
      if(type=="CENT"){
	ct_e_bemc.push_back(0.0); 
	ct_e_ihcal.push_back(stage1_energy); 
      }
      else if(type=="FWD"){
	ct_e_femc.push_back(0.0); 
	ct_e_lfhcal.push_back(stage1_energy); 
      }

      // Last calorimeter

      double stage3_energy = 0.0; 

      if(detName[2]!=""){

	for (unsigned int j = 0; j < clusterList[2]->size(); j++) {

	  if(!tmatched2[k]) continue;
	  if(cused2[j]) continue; 

	  // Look for same track match
	  if(tmatched2[j]==tmatched1[k]){

	    RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	    caloTot += rcluster2->get_energy();
	    stage3_energy += rcluster2->get_energy();
	  
	    cused2[j] = true; 

	    //break; 

	  }

	}

      }

      if(type=="CENT"){
	ct_e_ohcal.push_back(stage3_energy); 
      }

      // record the combined energy
      ct_e_tot.push_back(caloTot); 

    }

  }
 
  // Now look for tracks seeded by the last HCAL (no EMCAL/first HCAL cluster)

  if(detName[2]!=""){

    for (unsigned int k = 0; k < clusterList[2]->size(); k++) {

      if(!tmatched2[k]) continue;
      if(cused2[k]) continue; 

      RawCluster *rcluster = clusterList[2]->getCluster(k);

      // Get the truth information for this track

      bool found = false; 
      PHG4TruthInfoContainer::ConstRange range = _truth_container->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
	   truth_itr != range.second; ++truth_itr) {

	PHG4Particle* g4particle = truth_itr->second;
	if(!g4particle) {
	  LogDebug("");
	  continue;
	}

	if ((tmatched2[k]->get_truth_track_id() - g4particle->get_track_id()) == 0) {

	  found = true;

	  ct_pid.push_back(g4particle->get_pid()); 

	  // Get the final particle in the lab frame
	  CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),
				      g4particle->get_pz(),g4particle->get_e());
	  efp = EventToLab * efp;  
	  TLorentzVector match_lf(efp.px(), efp.py(), efp.pz(), efp.e());

	  ct_p_meas.push_back(tmatched2[k]->get_p()); 
	  ct_p_true.push_back(match_lf.Vect().Mag()); 
	
	  ct_eta_meas.push_back(tmatched2[k]->get_eta()); 
	  ct_eta_true.push_back(match_lf.Vect().Eta()); 

	  double deta = tmatched2[k]->get_eta() - match_lf.Vect().Eta(); 
	  double dphi = DeltaPhi(tmatched2[k]->get_phi(),match_lf.Vect().Phi()); 

	  ct_dist.push_back(sqrt(pow(dphi,2)+ pow(deta,2))); 

	  cused2[k] = true; 

	  break; 

	}
	
      }
    
      if(!found) continue; 

      double caloTot = rcluster->get_energy(); 
 
      if(type=="CENT"){
	ct_e_bemc.push_back(0.0); 
	ct_e_ihcal.push_back(0.0); 
	ct_e_ohcal.push_back(rcluster->get_energy()); 
      }

      // record the combined energy
      ct_e_tot.push_back(caloTot); 

    }

  }

  if(type=="CENT"){
    _eval_charged_tracks_cent->Fill(); 
  }
  else if(type=="FWD"){
    _eval_charged_tracks_fwd->Fill();  
  }

  return; 

}

void CentauroJets::FillClusterPseudoJets( PHCompositeNode *topNode, std::string detName, 
					  std::vector<fastjet::PseudoJet> &pseudojets, 
					  TLorentzRotation &breit, TRotation &breitRot, 
					  const double scale_factor, std::string ecDet, int ecIdx, bool TrackVeto){

  string clusternodename = "CLUSTER_" + detName;
  RawClusterContainer *clusterList = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusterList) {
    cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
    return;
  }    

  for (unsigned int k = 0; k < clusterList->size(); ++k) {

    // skip the electron cluster
    if( (ecDet==detName) && (ecIdx==(int)k) ) continue; 

    RawCluster *rcluster = clusterList->getCluster(k);

    //double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
    // match is to eta relative to z=0
    double eta = getEta(rcluster->get_r(),rcluster->get_z());
    double phi = rcluster->get_phi(); 

    // eliminate noise clusters
    if(rcluster->get_energy()<CLUSTER_E_CUTOFF) continue; 

    if(TrackVeto){
      if(VetoClusterWithTrack(eta, phi, detName)) continue; 
    }

    double pt = scale_factor*rcluster->get_energy() / cosh(eta);
    double px = pt * cos(phi);
    double py = pt * sin(phi);
    double pz = pt * sinh(eta);

    // Transform to Breit frame
    TLorentzVector cluster(px,py,pz,scale_factor*rcluster->get_energy()); 
    TLorentzVector breit_cluster = (breit*cluster); 
    breit_cluster.Transform(breitRot); 

    fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
    pseudojet.set_user_index(0); 
    pseudojets.push_back(pseudojet);

  }

  return; 

}

RawCluster *CentauroJets:: getCluster( PHCompositeNode *topNode, std::string detName, 
				       double eta, double phi, double &clustE, int &clIdx, double &dR){

  RawCluster *retCluster = NULL; 

  // pull the clusters
  string clusternodename = "CLUSTER_" + detName;
  RawClusterContainer *clusterList = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusterList) {
    cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
    clustE = -9999.0;
    clIdx = -1;
    return NULL; 
  }    

  // loop over all clusters and find nearest
  double min_r = 9999.0;
  double min_e = -9999.0;
  int min_idx = -1; 
  for (unsigned int k = 0; k < clusterList->size(); ++k) {

    RawCluster *cluster = clusterList->getCluster(k);

    double dphi = DeltaPhi(phi,cluster->get_phi());
    double cluster_eta = getEta(cluster->get_r(),cluster->get_z()); 
    double deta = eta-cluster_eta;
    double r = sqrt(pow(dphi,2)+pow(deta,2));

    if (r < min_r) {
      min_r = r; 
      min_e = cluster->get_energy();
      min_idx = (int)k;
      retCluster = cluster; 
   }

  }

  clustE = min_e; 
  clIdx = min_idx; 
  dR = min_r; 
  
  return retCluster; 

}

double CentauroJets::getClusterByIndex( PHCompositeNode *topNode, std::string detName, int clIdx){

  // pull the clusters
  string clusternodename = "CLUSTER_" + detName;
  RawClusterContainer *clusterList = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusterList) {
    cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
    return -9999.0;
  }    

  double min_e = -9999.0;
  if((clIdx>=0) && (clIdx<(int)clusterList->size())){
    RawCluster *cluster = clusterList->getCluster(clIdx);
    min_e = cluster->get_energy();
  }
  
  return min_e; 


}


//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
int CentauroJets::GetNodes(PHCompositeNode * topNode) {
	//DST objects
	//Truth container
	_truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
	if (!_truth_container && _event < 2) {
	  LogError(" PHG4TruthInfoContainer node not found on node tree");
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

	_trackmap = findNode::getClass<SvtxTrackMap>(topNode,"TrackMap");
	if (!_trackmap && _event < 2) {
	  LogError("TrackMap node not found on node tree");
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

void CentauroJets::GetPrimaryJets(PHCompositeNode *topNode, fastjet::JetDefinition* jetdef, std::vector<fastjet::PseudoJet> *fastjets, 
				  TLorentzRotation &breit, TRotation &breitRot, 
				  TLorentzVector p_initial_breit, TLorentzVector virtual_photon_breit, bool true_frame){

  std::vector<fastjet::PseudoJet> primary_jets;

  // get the list of jets from the primary particles

  if (!_truth_container) {
    LogError("_truth_container not found!");
    return;
  }

  std::vector<fastjet::PseudoJet> pseudojets;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // PRIMARIES ONLY
  PHG4TruthInfoContainer::ConstRange range =
   		_truth_container->GetPrimaryParticleRange();

  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr) {

    PHG4Particle* g4particle = truth_itr->second;
    if(!g4particle) {
      LogDebug("");
      continue;
    }

    // Skip the scattered electron:
    if ((g4particle->get_track_id() - true_electron_headon->get_track_id()) == 0) continue;  

    // remove some particles (muons, taus, neutrinos)...
    // 12 == nu_e
    // 13 == muons
    // 14 == nu_mu
    // 15 == taus
    // 16 == nu_tau
    if ((abs(g4particle->get_pid()) >= 12) && (abs( g4particle->get_pid()) <= 16)) continue;
        
    // NOTE: stored HepMC kinematics are in event generator (head-on) frame!
    // Transform to the lab frame
    CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),g4particle->get_pz(),g4particle->get_e());
    efp = EventToLab * efp;  

    TLorentzVector partMom(efp.px(), efp.py(), efp.pz(), efp.e()); 

    // lab frame cuts
    if((partMom.Pt()<0.200)||(fabs(partMom.Eta())>4.0)) continue;

    TLorentzVector partMom_breit = (breit*partMom); 
    partMom_breit.Transform(breitRot); 

    // add this track to the list of tracks for jets

    fastjet::PseudoJet pseudojet (partMom_breit.Px(),
				  partMom_breit.Py(),
				  partMom_breit.Pz(),
				  partMom_breit.E());

    G4ParticleDefinition* particle = particleTable->FindParticle(g4particle->get_name());
    int charge = -9999.0; 
    if(particle) 
      charge = particle->GetPDGCharge();
    else
      continue; 

    // build the user index

    bool em_part = false; 
    if((abs(g4particle->get_pid())==11) || (g4particle->get_pid()==22)) em_part = true; 

    pseudojet.set_user_index(EncodeUserIndex(charge,em_part));
    pseudojets.push_back(pseudojet);

  }

  if(pseudojets.size()>0) { 

    // Call FastJet
    fastjet::ClusterSequence jetFinder(pseudojets,*jetdef);
    *fastjets = sorted_by_E(jetFinder.inclusive_jets(0.0));

    for (unsigned int ijet = 0; ijet < fastjets->size(); ++ijet) {

      if(true_frame){

	// tag higest E jet
	if(ijet==0)
	  tfpjet_largest.push_back(1); 
	else
	  tfpjet_largest.push_back(0);
      
	tfpjet_pT.push_back(fastjets->at(ijet).perp()); 
	tfpjet_p.push_back(sqrt(pow(fastjets->at(ijet).px(),2) + 
			      pow(fastjets->at(ijet).py(),2) + 
			      pow(fastjets->at(ijet).pz(),2)));
	tfpjet_E.push_back(fastjets->at(ijet).E()); 
	tfpjet_eta.push_back(fastjets->at(ijet).eta()); 
	tfpjet_phi.push_back(fastjets->at(ijet).phi_02pi()); 
  
	TLorentzVector pjet(fastjets->at(ijet).px(),fastjets->at(ijet).py(),fastjets->at(ijet).pz(),fastjets->at(ijet).E());
	double zjet = (fastjets->at(ijet).E()-fastjets->at(ijet).pz())/sqrt(_hepmcp_Q2);
	tfpjet_z.push_back(zjet);

	tfpjet_qperp.push_back(fastjets->at(ijet).perp()/zjet); 

	tfpjet_nc.push_back(fastjets->at(ijet).constituents().size());

	// Jet Charge
	std::vector<fastjet::PseudoJet> tconstit = fastjets->at(ijet).constituents();

	double jptot = sqrt(pow(fastjets->at(ijet).px(),2) + 
			    pow(fastjets->at(ijet).py(),2) + 
			    pow(fastjets->at(ijet).pz(),2)); 

	tfpjet_Q.push_back(JetCharge(&tconstit,jptot));

	tfpjet_cf.push_back(JetChargedFraction(&tconstit,jptot));

	tfpjet_neut_p.push_back(JetNeutralMomentum(&tconstit)); 
	tfpjet_chgd_p.push_back(JetChargedMomentum(&tconstit)); 
	tfpjet_em_p.push_back(JetEMMomentum(&tconstit)); 

      }
      else{

	// tag higest E jet
	if(ijet==0)
	  pjet_largest.push_back(1); 
	else
	  pjet_largest.push_back(0);
      
	pjet_pT.push_back(fastjets->at(ijet).perp()); 
	pjet_p.push_back(sqrt(pow(fastjets->at(ijet).px(),2) + 
			      pow(fastjets->at(ijet).py(),2) + 
			      pow(fastjets->at(ijet).pz(),2)));
	pjet_E.push_back(fastjets->at(ijet).E()); 
	pjet_eta.push_back(fastjets->at(ijet).eta()); 
	pjet_phi.push_back(fastjets->at(ijet).phi_02pi()); 
  
	TLorentzVector pjet(fastjets->at(ijet).px(),fastjets->at(ijet).py(),fastjets->at(ijet).pz(),fastjets->at(ijet).E());
	double zjet = (fastjets->at(ijet).E()-fastjets->at(ijet).pz())/sqrt(measQ2); 
	pjet_z.push_back(zjet);

	pjet_qperp.push_back(fastjets->at(ijet).perp()/zjet); 
  
	pjet_nc.push_back(fastjets->at(ijet).constituents().size());

	// Jet Charge
	std::vector<fastjet::PseudoJet> tconstit = fastjets->at(ijet).constituents();

	double jptot = sqrt(pow(fastjets->at(ijet).px(),2) + 
			    pow(fastjets->at(ijet).py(),2) + 
			    pow(fastjets->at(ijet).pz(),2)); 

	pjet_Q.push_back(JetCharge(&tconstit,jptot));

	pjet_cf.push_back(JetChargedFraction(&tconstit,jptot));

	pjet_neut_p.push_back(JetNeutralMomentum(&tconstit)); 
	pjet_chgd_p.push_back(JetChargedMomentum(&tconstit)); 
	pjet_em_p.push_back(JetEMMomentum(&tconstit)); 

      }

   }

  }


}

bool CentauroJets::isThisAnElectron(SvtxTrack *temp)
{

  // Get the truth information for this track

  bool isElectron = false; 
 
  PHG4TruthInfoContainer::ConstRange range = _truth_container->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr) {

    PHG4Particle* g4particle = truth_itr->second;
    if(!g4particle) {
      LogDebug("");
      continue;
    }

    if ((temp->get_truth_track_id() - g4particle->get_track_id()) == 0) {

      if(abs(g4particle->get_pid())==11) isElectron = true; 
      break; 

    }

  }

  return isElectron; 

}
