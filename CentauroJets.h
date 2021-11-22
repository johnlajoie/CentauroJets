
#ifndef __CentauroJets_H__
#define __CentauroJets_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle; 
class SvtxTrackMap;
class SvtxTrack; 
class TFile;
class TTree;
class TRandom; 
class TLorentzRotation; 
class TRotation; 
class RawTowerContainer; 
class RawTowerGeomContainer; 
class TLorentzVector; 
class RawCluster; 
class TH1D; 
class TH2D; 

namespace fastjet {
  class PseudoJet;
  class JetDefinition; 
}

#include <CLHEP/Vector/LorentzRotation.h>


//Brief: basic ntuple and histogram creation for sim evaluation
class CentauroJets: public SubsysReco
{
 public: 
  //Default constructor
  CentauroJets(const std::string &name="CentauroJets");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *); 

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile_name = file; }

  //User modules
  void fill_tree(PHCompositeNode*);

 private:
  //output filename
  std::string _outfile_name;
   
  //Event counter
  int _event;

  //Get all the nodes
  int GetNodes(PHCompositeNode *);
  
  // Match track to cluster
  double getClusterByIndex( PHCompositeNode *topNode, std::string detName, int clIdx);
  RawCluster *getCluster( PHCompositeNode *topNode, std::string detName, double eta, double phi, double &clustE, int &clIdx, double &dR); 

  // Fill tower pseudojets
  void FillTowerPseudoJets( PHCompositeNode *topNode, std::string detName, 
			    std::vector<fastjet::PseudoJet> &pseudojets, 
			    TLorentzRotation &breit, TRotation &breitRot, 
			    const double scale_factor, std::string ecDet = "", RawCluster *clustPtr = NULL); 

  // Fill Cluster pseudojets
  void FillClusterPseudoJets( PHCompositeNode *topNode, std::string detName, 
			      std::vector<fastjet::PseudoJet> &pseudojets, 
			      TLorentzRotation &breit, TRotation &breitRot, 
			      const double scale_factor, std::string ecDet = "", 
			      int ecIdx = -1, bool trackVeto = false);

  // Primary jets
  void GetPrimaryJets(PHCompositeNode *topNode,  fastjet::JetDefinition* jetdef,  std::vector<fastjet::PseudoJet> *fastjets,
		      TLorentzRotation &breit, TRotation &breitRot, 
		      TLorentzVector p_initial_breit, TLorentzVector virtual_photon_breit, bool true_frame); 

  // Build Calo Tracks (combine clusters across calorimeters)
  void BuildCaloTracks(PHCompositeNode *topNode, std::string type, 
		       std::vector<fastjet::PseudoJet> &pseudojets, 
		       TLorentzRotation &breit, TRotation &breitRot, 
		       std::string ecDet, int ecIdx );

  void BuildChargedCaloTracks(PHCompositeNode *topNode, std::string type);

  bool VetoClusterWithTrack(double eta, double phi, std::string detName);

  SvtxTrack *AttachClusterToTrack(double eta, double phi, std::string detName); 

  bool isThisAnElectron(SvtxTrack *temp); 

  // Event generator transform
  CLHEP::HepLorentzRotation EventToLab; 

  TTree *_eval_tree_event; 

  double e_p_initial; 
  double p_p_initial; 

  double _hepmcp_x1; 
  double _hepmcp_x2; 
  double _hepmcp_Q2; 
  double _hepmcp_procid; 
  int _hepmcp_id1; 
  int _hepmcp_id2; 

  int event; 
  double measQ2; 
  double meas_x; 
  double meas_E_p; 
  double electron_eta; 
  double electron_cluster_dR; 

  double vtx_x; 
  double vtx_y; 
  double vtx_z; 
  double vtx_t; 

  double breit_vphot_e; 
  double breit_initial_proton_e;
  double breit_vphot_pz; 
  double breit_vphot_pt;
  double breit_initial_proton_pz; 
  double breit_initial_proton_pt;

  std::vector<double> jet_pT; 
  std::vector<double> jet_p; 
  std::vector<double> jet_E; 
  std::vector<double> jet_eta; 
  std::vector<double> jet_phi; 
  std::vector<double> jet_z; 
  std::vector<double> jet_qperp; 
  std::vector<int> jet_largest; 
  std::vector<int> jet_pidx; 
  std::vector<double> jet_lab_eta; 
  std::vector<double> jet_lab_phi; 
  std::vector<double> jet_lab_p; 
  std::vector<double> jet_pdR; 
  std::vector<int> jet_nc; 

  std::vector<double> tjet_pT; 
  std::vector<double> tjet_p; 
  std::vector<double> tjet_E; 
  std::vector<double> tjet_eta; 
  std::vector<double> tjet_phi; 
  std::vector<double> tjet_z; 
  std::vector<double> tjet_qperp; 
  std::vector<int> tjet_largest; 
  std::vector<int> tjet_pidx; 
  std::vector<double> tjet_lab_eta; 
  std::vector<double> tjet_lab_phi; 
  std::vector<double> tjet_lab_p; 
  std::vector<double> tjet_pdR; 
  std::vector<int> tjet_nc; 
  std::vector<double> tjet_Q; 

  std::vector<double> tcjet_pT; 
  std::vector<double> tcjet_p; 
  std::vector<double> tcjet_E; 
  std::vector<double> tcjet_eta; 
  std::vector<double> tcjet_phi; 
  std::vector<double> tcjet_z; 
  std::vector<double> tcjet_qperp; 
  std::vector<int> tcjet_largest; 
  std::vector<int> tcjet_pidx; 
  std::vector<double> tcjet_lab_eta; 
  std::vector<double> tcjet_lab_phi; 
  std::vector<double> tcjet_lab_p; 
  std::vector<double> tcjet_pdR; 
  std::vector<int> tcjet_nc; 
  std::vector<double> tcjet_Q; 
  std::vector<double> tcjet_cf; 
  std::vector<double> tcjet_neut_p; 
  std::vector<double> tcjet_chgd_p; 
  std::vector<double> tcjet_em_p; 

  std::vector<double> pjet_pT; 
  std::vector<double> pjet_p; 
  std::vector<double> pjet_E; 
  std::vector<double> pjet_eta; 
  std::vector<double> pjet_phi; 
  std::vector<double> pjet_z; 
  std::vector<double> pjet_qperp; 
  std::vector<int> pjet_largest; 
  std::vector<int> pjet_nc; 
  std::vector<double> pjet_Q; 
  std::vector<double> pjet_cf; 
  std::vector<int> pjet_tcidx; 
  std::vector<double> pjet_tcdR; 
  std::vector<double> pjet_neut_p; 
  std::vector<double> pjet_chgd_p; 
  std::vector<double> pjet_em_p; 

  std::vector<double> tfpjet_pT; 
  std::vector<double> tfpjet_p; 
  std::vector<double> tfpjet_E; 
  std::vector<double> tfpjet_eta; 
  std::vector<double> tfpjet_phi; 
  std::vector<double> tfpjet_z; 
  std::vector<double> tfpjet_qperp; 
  std::vector<int> tfpjet_largest; 
  std::vector<int> tfpjet_nc;   
  std::vector<double> tfpjet_Q; 
  std::vector<double> tfpjet_cf; 
  std::vector<double> tfpjet_neut_p; 
  std::vector<double> tfpjet_chgd_p; 
  std::vector<double> tfpjet_em_p; 

  TTree *_eval_charged_tracks_cent; 
  TTree *_eval_charged_tracks_fwd; 

  std::vector<int> ct_pid; 
  std::vector<double> ct_p_meas; 
  std::vector<double> ct_p_true; 
  std::vector<double> ct_p_true_lab; 
  std::vector<double> ct_eta_meas; 
  std::vector<double> ct_eta_true; 
  std::vector<double> ct_dist; 
  std::vector<double> ct_e_bemc; 
  std::vector<double> ct_e_ihcal; 
  std::vector<double> ct_e_ohcal; 
  std::vector<double> ct_e_femc; 
  std::vector<double> ct_e_lfhcal; 
  std::vector<double> ct_e_tot; 

  TTree *_eval_calo_tracks_cent; 
  TTree *_eval_calo_tracks_fwd; 
  std::vector<int> cat_pid;  
  std::vector<double> cat_p_true; 
  std::vector<double> cat_p_true_lab; 
  std::vector<double> cat_eta_meas; 
  std::vector<double> cat_eta_true; 
  std::vector<double> cat_e_tot; 
  std::vector<double> cat_match; 

  // Scattered electron in event record:
  PHG4Particle* true_electron_headon; 

  //Node pointers
  PHG4TruthInfoContainer* _truth_container;
  SvtxTrackMap* _trackmap;

  //Random number generator
  TRandom *rand; 

  //Diagnostic histograms
  TH1D *_h_track_cluster_match; 
  TH1D *_h_track_cluster_match_becal; 
  TH1D *_h_track_cluster_match_ihcal; 
  TH1D *_h_track_cluster_match_ohcal; 
  TH1D *_h_track_cluster_match_femc; 
  TH1D *_h_track_cluster_match_lfhcal; 

  TH1D *_h_becal_ihcal_match; 
  TH1D *_h_becal_ihcal_match_eta; 
  TH1D *_h_becal_ihcal_match_phi; 
  TH1D *_h_becal_ohcal_match; 
  TH1D *_h_becal_ohcal_match_eta; 
  TH1D *_h_becal_ohcal_match_phi; 
  TH1D *_h_ihcal_ohcal_match; 
  TH1D *_h_ihcal_ohcal_match_eta; 
  TH1D *_h_ihcal_ohcal_match_phi; 
  TH1D *_h_femc_lfhcal_match;
  TH1D *_h_femc_lfhcal_match_eta;
  TH1D *_h_femc_lfhcal_match_phi;

  TH1D *_h_calotrack_prim_match_cent; 
  TH1D *_h_calotrack_prim_match_cent_gamma; 
  TH1D *_h_calotrack_prim_match_cent_neutron; 
  TH1D *_h_calotrack_pid_prim_match_cent; 

  TH1D *_h_calotrack_prim_match_fwd; 
  TH1D *_h_calotrack_prim_match_fwd_gamma; 
  TH1D *_h_calotrack_prim_match_fwd_neutron; 
  TH1D *_h_calotrack_pid_prim_match_fwd; 

  TH1D *_h_calotrack_cent_NES; 
  TH1D *_h_calotrack_fwd_NES; 

  TH2D *_h_calotrack_cent_NES_2D; 
  TH2D *_h_calotrack_fwd_NES_2D; 

  TH2D *_h_calotrack_cent_gamma_NES_2D; 
  TH2D *_h_calotrack_fwd_gamma_NES_2D; 

  TH2D *_h_calotrack_cent_neutron_NES_2D; 
  TH2D *_h_calotrack_fwd_neutron_NES_2D; 

  TH1D *_h_nclusters_becal; 
  TH1D *_h_nclusters_ihcal; 
  TH1D *_h_nclusters_ohcal; 
  TH1D *_h_nclusters_femc; 
  TH1D *_h_nclusters_lfhcal; 
  

};

#endif //* __CentauroJets_H__ *//
