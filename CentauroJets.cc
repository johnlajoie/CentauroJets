
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

// PID includes
#include <eicpidbase/EICPIDParticle.h>
#include <eicpidbase/EICPIDParticleContainer.h>
#include <g4eval/JetEvalStack.h>
#include <g4eval/SvtxEvalStack.h>

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
#define EM_CLUSTER_E_CUTOFF 0.100
#define HAD_CLUSTER_E_CUTOFF 0.200

// Use PID?
#define USE_PID 0

// Use parametrized photon/neutral hadron response?
//#define PARAM_PHOTONS 1
#define PARAM_HAD_NEUTRALS 1

double getMatchingCut(double p, int q, std::string detName){


  double retCut = 0.15; 

  if(detName=="BECAL"){

    if(q<0){

      if(p<=0.75)
	retCut = 0.08; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.065; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.065; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.04; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.03; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.025; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.02; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.02; 
      else
	retCut = 0.02;  

    }
    else if(q>0){
  
      if(p<=0.75)
	retCut = 0.1; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.09; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.08; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.06; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.05; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.035; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.03; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.025; 
      else
	retCut = 0.02;  

    }
    else{


    }

  }

  if(detName=="EEMC"){

    if(q<0){

     if(p<=0.75)
	retCut = 0.15; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.125; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.08; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.06; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.05; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.05; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.045; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.04; 
      else
	retCut = 0.035;  

    }
    else if(q>0){
  
      if(p<=0.75)
	retCut = 0.15; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.125; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.08; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.06; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.05; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.05; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.04; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.035; 
      else
	retCut = 0.035;  

    }
    else{


    }

  }

  if(detName=="FEMC"){

    if(q<0){

      if(p<=0.75)
	retCut = 0.1; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.08; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.08; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.05; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.04; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.035; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.035; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.03; 
      else
	retCut = 0.03;  

    }
    else if(q>0){
  
      if(p<=0.75)
	retCut = 0.125; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.10; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.08; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.05; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.04; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.035; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.03; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.025; 
      else
	retCut = 0.02;  

    }
    else{




    }

  }

  if(detName=="HCALIN"){

    if(q<0){

      if(p<=0.75)
	retCut = 0.35; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.35; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.35; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.35; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.35; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.35; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.35; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.35; 
      else
	retCut = 0.35;  

    }
    else if(q>0){
  
      if(p<=0.75)
	retCut = 0.5; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.5; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.4; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.4; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.4; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.4; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.4; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.4; 
      else
	retCut = 0.4;  

    }
    else{




    }

  }

  if(detName=="HCALOUT"){

    if(q<0){

      if(p<=0.75)
	retCut = 0.5; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.5; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.4; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.4; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.4; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.4; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.4; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.4; 
      else
	retCut = 0.4;  

    }
    else if(q>0){
  
      if(p<=0.75)
	retCut = 0.5; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.5; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.4; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.4; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.4; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.4; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.4; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.4; 
      else
	retCut = 0.4;  

    }
    else{




    }

  }

  if(detName=="LFHCAL"){

    if(q<0){

      if(p<=0.75)
	retCut = 0.5; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.5; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.5; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.5; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.5; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.5; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.5; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.5; 
      else
	retCut = 0.6;  

    }
    else if(q>0){
 
      if(p<=0.75)
	retCut = 0.5; 
      else if(p>0.75 && p<=1.0)
	retCut = 0.5; 
      else if(p>1.0 && p<=2.0)
	retCut = 0.5; 
      else if(p>2.0 && p<=3.0)
	retCut = 0.5; 
      else if(p>3.0 && p<=4.0)
	retCut = 0.5; 
      else if(p>4.0 && p<=5.0)
	retCut = 0.5; 
      else if(p>5.0 && p<=10.0)
	retCut = 0.5; 
      else if(p>10.0 && p<=15.0)
	retCut = 0.5; 
      else
	retCut = 0.6;   

    }
    else{




    }

  }

  // return cut value 

  return retCut; 


}


double GetCaloTrackHadronicEnergyScale(double e[3], std::string type){

  double scale = 1.0; 
  double e_tot = e[0] + e[1] + e[2]; 

  if(type=="CENT"){

    if((e[0]>0.0 || e[1]>0.0) && e[2]>0.0){
      
      if(e_tot>0.215){
        double p_true = (1.0/0.448)*(e_tot - 0.215); 
        scale = p_true/e_tot;
      }
      else{
	scale = 0.0; 
      }
      
    }
    else if(e[0]==0.0 && e[1] == 0.0 && e[2]>0.0){

      if(e_tot>1.47854){
	double p_true = (1.0/0.458695)*(e_tot - 1.47852); 
	scale = p_true/e_tot;
      }
      else{
	scale = 0.0; 
      }

    }
    
  }
  else if(type=="FWD"){
      
    if(e_tot>0.5706){
      double p_true = (1.0/0.38555)*(e_tot - 0.5706); 
      scale = p_true/e_tot;
    }
    else{
      scale = 0.0; 
    }
      
  }

  return scale; 

}



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

TVector3 JetNeutralMomentum( std::vector<fastjet::PseudoJet> *tconstit ){

  // neutral hadrons

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    // eliminate charged
    if(GetChargeFromUserIndex(tconstit->at(i).user_index())!=0) continue;
    // eliminate photons
    if(tconstit->at(i).user_index()==10) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot; 

}

TVector3 JetChargedMomentum( std::vector<fastjet::PseudoJet> *tconstit ){

  // charged hadrons

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    if(GetChargeFromUserIndex(tconstit->at(i).user_index())==0) continue;
    // eliminate electrons/positrons
    if( (tconstit->at(i).user_index()==11) || (tconstit->at(i).user_index()==12) ) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot; 

}

double JetChargedFraction( std::vector<fastjet::PseudoJet> *tconstit, TVector3 pjet ){

  // photons. electrons, positrons

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    if(GetChargeFromUserIndex(tconstit->at(i).user_index())==0) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot.Mag()/pjet.Mag(); 

}

TVector3 JetEMMomentum( std::vector<fastjet::PseudoJet> *tconstit ){

  TVector3 ptot(0.0,0.0,0.0); 

  for(unsigned int i=0; i<tconstit->size(); i++){
    if(!isEMParticle(tconstit->at(i).user_index())) continue;
    TVector3 constit(tconstit->at(i).px(),tconstit->at(i).py(),tconstit->at(i).pz());
    ptot += constit; 
  }

  return ptot; 

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
	_eval_tree_event->Branch("y", &_hepmcp_y, "_hepmcp_y/D");
	_eval_tree_event->Branch("W2", &_hepmcp_W2, "_hepmcp_W2/D");
	_eval_tree_event->Branch("procid", &_hepmcp_procid, "_hepmcp_procid/I");
	_eval_tree_event->Branch("id1", &_hepmcp_id1, "_hepmcp_id1/I");
	_eval_tree_event->Branch("id2", &_hepmcp_id2, "_hepmcp_id2/I");
	_eval_tree_event->Branch("vtx_x", &vtx_x, "vtx_x/D");
	_eval_tree_event->Branch("vtx_y", &vtx_y, "vtx_y/D");
	_eval_tree_event->Branch("vtx_z", &vtx_z, "vtx_z/D");
	_eval_tree_event->Branch("vtx_t", &vtx_t, "vtx_t/D");	
	_eval_tree_event->Branch("measQ2", &measQ2, "measQ2/D");
	_eval_tree_event->Branch("meas_x", &meas_x, "meas_x/D");
	_eval_tree_event->Branch("meas_y", &meas_y, "meas_y/D");
	_eval_tree_event->Branch("measW2", &measW2, "measW2/D");
	_eval_tree_event->Branch("meas_E_p", &meas_E_p, "meas_E_p/D"); 
	_eval_tree_event->Branch("electron_eta",&electron_eta,"electron_eta/D"); 
	_eval_tree_event->Branch("electron_phi",&electron_eta,"electron_phi/D"); 
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
	_eval_tree_event->Branch("jet_lab_p",&jet_lab_p); 
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
	_eval_tree_event->Branch("tjet_lab_p",&tjet_lab_p); 
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
	_eval_tree_event->Branch("tcjet_lab_p",&tcjet_lab_p); 
	_eval_tree_event->Branch("tcjet_nc",&tcjet_nc);
	_eval_tree_event->Branch("tcjet_Q",&tcjet_Q);
	_eval_tree_event->Branch("tcjet_cf",&tcjet_cf);
	_eval_tree_event->Branch("tcjet_neut_p",&tcjet_neut_p); 
	_eval_tree_event->Branch("tcjet_chgd_p",&tcjet_chgd_p); 
	_eval_tree_event->Branch("tcjet_em_p",&tcjet_em_p); 
	_eval_tree_event->Branch("tcjet_neut_pm",&tcjet_neut_pm); 
	_eval_tree_event->Branch("tcjet_chgd_pm",&tcjet_chgd_pm); 
	_eval_tree_event->Branch("tcjet_em_pm",&tcjet_em_pm); 

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
	_eval_tree_event->Branch("pjet_neut_pm",&pjet_neut_pm); 
	_eval_tree_event->Branch("pjet_chgd_pm",&pjet_chgd_pm); 
	_eval_tree_event->Branch("pjet_em_pm",&pjet_em_pm); 
	_eval_tree_event->Branch("pjet_lab_eta",&pjet_lab_eta); 
	_eval_tree_event->Branch("pjet_lab_phi",&pjet_lab_phi); 
	_eval_tree_event->Branch("pjet_lab_p",&pjet_lab_p); 
 
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
	_eval_tree_event->Branch("tfpjet_tcidx",&tfpjet_tcidx);
	_eval_tree_event->Branch("tfpjet_tcdR",&tfpjet_tcdR);
	_eval_tree_event->Branch("tfpjet_neut_p",&tfpjet_neut_p); 
	_eval_tree_event->Branch("tfpjet_chgd_p",&tfpjet_chgd_p); 
	_eval_tree_event->Branch("tfpjet_em_p",&tfpjet_em_p); 
	_eval_tree_event->Branch("tfpjet_neut_pm",&tfpjet_neut_pm); 
	_eval_tree_event->Branch("tfpjet_chgd_pm",&tfpjet_chgd_pm); 
	_eval_tree_event->Branch("tfpjet_em_pm",&tfpjet_em_pm); 
	_eval_tree_event->Branch("tfpjet_lab_eta",&tfpjet_lab_eta); 
	_eval_tree_event->Branch("tfpjet_lab_phi",&tfpjet_lab_phi); 
	_eval_tree_event->Branch("tfpjet_lab_p",&tfpjet_lab_p); 

	// cluster evaluation for charged tracks
	_eval_charged_tracks_cent = new TTree("clusteval_cent", "Charged track clusters (central)");
	_eval_charged_tracks_cent->Branch("event", &event, "event/I");
	_eval_charged_tracks_cent->Branch("pid",&ct_pid,"ct_pid/I"); 
	_eval_charged_tracks_cent->Branch("p_meas",&ct_p_meas,"ct_p_meas/D"); 
	_eval_charged_tracks_cent->Branch("p_true",&ct_p_true,"ct_p_true/D"); 
	_eval_charged_tracks_cent->Branch("eta_meas",&ct_eta_meas,"ct_eta_meas/D"); 
	_eval_charged_tracks_cent->Branch("eta_true",&ct_eta_true,"ct_eta_true/D"); 
	_eval_charged_tracks_cent->Branch("phi_meas",&ct_phi_meas,"ct_phi_meas/D"); 
	_eval_charged_tracks_cent->Branch("phi_true",&ct_phi_true,"ct_phi_true/D"); 
	_eval_charged_tracks_cent->Branch("e_bemc",&ct_e_bemc,"ct_e_bemc/D"); 
	_eval_charged_tracks_cent->Branch("e_ihcal",&ct_e_ihcal,"ct_e_ihcal/D"); 
	_eval_charged_tracks_cent->Branch("e_ohcal",&ct_e_ohcal,"ct_e_ohcal/D"); 
	_eval_charged_tracks_cent->Branch("e_tot",&ct_e_tot,"ct_e_tot/D"); 

	_eval_charged_tracks_fwd = new TTree("clusteval_fwd", "Charged track clusters (fwd)");
	_eval_charged_tracks_fwd->Branch("event", &event, "event/I");
	_eval_charged_tracks_fwd->Branch("pid",&ct_pid,"c_pid/I"); 
	_eval_charged_tracks_fwd->Branch("p_meas",&ct_p_meas,"ct_p_meas/D"); 
	_eval_charged_tracks_fwd->Branch("p_true",&ct_p_true,"ct_p_true/D"); 
	_eval_charged_tracks_fwd->Branch("eta_meas",&ct_eta_meas,"ct_ema_meas/D"); 
	_eval_charged_tracks_fwd->Branch("eta_true",&ct_eta_true,"ct_eta_true/D"); 
	_eval_charged_tracks_fwd->Branch("phi_meas",&ct_phi_meas,"ct_phi_meas/D"); 
	_eval_charged_tracks_fwd->Branch("phi_true",&ct_phi_true,"ct_phi_true/D"); 
	_eval_charged_tracks_fwd->Branch("e_femc",&ct_e_femc,"ct_e_femc/D"); 
	_eval_charged_tracks_fwd->Branch("e_lfhcal",&ct_e_lfhcal,"ct_e_lfhcal/D"); 
	_eval_charged_tracks_fwd->Branch("e_tot",&ct_e_tot,"ct_e_tot/D"); 

	_eval_charged_tracks_bkwd = new TTree("clusteval_bkwd", "Charged track clusters (fwd)");
	_eval_charged_tracks_bkwd->Branch("event", &event, "event/I");
	_eval_charged_tracks_bkwd->Branch("pid",&ct_pid,"c_pid/I"); 
	_eval_charged_tracks_bkwd->Branch("p_meas",&ct_p_meas,"ct_p_meas/D"); 
	_eval_charged_tracks_bkwd->Branch("p_true",&ct_p_true,"ct_p_true/D"); 
	_eval_charged_tracks_bkwd->Branch("eta_meas",&ct_eta_meas,"ct_ema_meas/D"); 
	_eval_charged_tracks_bkwd->Branch("eta_true",&ct_eta_true,"ct_eta_true/D"); 
	_eval_charged_tracks_bkwd->Branch("phi_meas",&ct_phi_meas,"ct_phi_meas/D"); 
	_eval_charged_tracks_bkwd->Branch("phi_true",&ct_phi_true,"ct_phi_true/D"); 
	_eval_charged_tracks_bkwd->Branch("e_eemc",&ct_e_eemc,"ct_e_eemc/D"); 

	// evaluation for calo tracks
	_eval_calo_tracks_cent = new TTree("calotrackeval_cent", "Calotrack Evaluation (central)");
	_eval_calo_tracks_cent->Branch("event", &event, "event/I");
	_eval_calo_tracks_cent->Branch("pid",&cat_pid,"cat_pid/I"); 
	_eval_calo_tracks_cent->Branch("p_true",&cat_p_true, "cat_p_true/D"); 
	_eval_calo_tracks_cent->Branch("eta_meas",&cat_eta_meas, "cat_eta_meas/D"); 
	_eval_calo_tracks_cent->Branch("eta_true",&cat_eta_true, "cat_eta_true/D"); 
	_eval_calo_tracks_cent->Branch("phi_meas",&cat_phi_meas, "cat_phi_meas/D"); 
	_eval_calo_tracks_cent->Branch("phi_true",&cat_phi_true, "cat_phi_true/D"); 
	_eval_calo_tracks_cent->Branch("match",&cat_match, "cat_match/D"); 
	_eval_calo_tracks_cent->Branch("e_tot",&cat_e_tot, "cat_e_tot/D"); 
	_eval_calo_tracks_cent->Branch("e_bemc",&cat_e_bemc,"cat_e_bemc/D"); 
	_eval_calo_tracks_cent->Branch("e_ihcal",&cat_e_ihcal,"cat_e_ihcal/D"); 
	_eval_calo_tracks_cent->Branch("e_ohcal",&cat_e_ohcal,"cat_e_ohcal/D"); 

	_eval_calo_tracks_fwd = new TTree("calotrackeval_fwd", "Calotrack Evaluation (fwd)");
	_eval_calo_tracks_fwd->Branch("event", &event, "event/I");
	_eval_calo_tracks_fwd->Branch("pid",&cat_pid,"cat_pid/I"); 
	_eval_calo_tracks_fwd->Branch("p_true",&cat_p_true, "cat_p_true/D"); 
	_eval_calo_tracks_fwd->Branch("eta_meas",&cat_eta_meas, "cat_eta_meas/D"); 
	_eval_calo_tracks_fwd->Branch("eta_true",&cat_eta_true, "cat_eta_true/D"); 
	_eval_calo_tracks_fwd->Branch("phi_meas",&cat_phi_meas, "cat_phi_meas/D"); 
	_eval_calo_tracks_fwd->Branch("phi_true",&cat_phi_true, "cat_phi_true/D"); 
	_eval_calo_tracks_fwd->Branch("match",&cat_match, "cat_match/D"); 
	_eval_calo_tracks_fwd->Branch("e_tot",&cat_e_tot, "cat_e_tot/D"); 
	_eval_calo_tracks_fwd->Branch("e_femc",&cat_e_femc,"cat_e_femc/D"); 
	_eval_calo_tracks_fwd->Branch("e_lfhcal",&cat_e_lfhcal,"cat_e_lfhcal/D"); 

	_eval_calo_tracks_bkwd = new TTree("calotrackeval_bkwd", "Calotrack Evaluation (bkwd)");
	_eval_calo_tracks_bkwd->Branch("event", &event, "event/I");
	_eval_calo_tracks_bkwd->Branch("pid",&cat_pid,"cat_pid/I"); 
	_eval_calo_tracks_bkwd->Branch("p_true",&cat_p_true, "cat_p_true/D"); 
	_eval_calo_tracks_bkwd->Branch("eta_meas",&cat_eta_meas, "cat_eta_meas/D"); 
	_eval_calo_tracks_bkwd->Branch("eta_true",&cat_eta_true, "cat_eta_true/D"); 
	_eval_calo_tracks_bkwd->Branch("phi_meas",&cat_phi_meas, "cat_phi_meas/D"); 
	_eval_calo_tracks_bkwd->Branch("phi_true",&cat_phi_true, "cat_phi_true/D"); 
	_eval_calo_tracks_bkwd->Branch("match",&cat_match, "cat_match/D"); 
	_eval_calo_tracks_bkwd->Branch("e_tot",&cat_e_tot, "cat_e_tot/D"); 
	_eval_calo_tracks_bkwd->Branch("e_eemc",&cat_e_eemc,"cat_e_eemc/D"); 

	// charged track matching 
	_eval_tmatch_eemc = new TTree("tmatch_eemc", "EEMC Track Match");
	_eval_tmatch_eemc->Branch("deta", &_tm_deta, "_tm_deta/D");
	_eval_tmatch_eemc->Branch("dphi", &_tm_dphi, "_tm_dphi/D");
	_eval_tmatch_eemc->Branch("dist", &_tm_dist, "_tm_dist/D");
	_eval_tmatch_eemc->Branch("eta", &_tm_eta, "_tm_eta/D");
	_eval_tmatch_eemc->Branch("phi", &_tm_phi, "_tm_phi/D");
	_eval_tmatch_eemc->Branch("p", &_tm_p, "_tm_p/D");
	_eval_tmatch_eemc->Branch("ceta", &_tm_ceta, "_tm_ceta/D");
	_eval_tmatch_eemc->Branch("cphi", &_tm_cphi, "_tm_cphi/D");
	_eval_tmatch_eemc->Branch("q", &_tm_q, "_tm_q/I");
	_eval_tmatch_eemc->Branch("e", &_tm_e, "_tm_e/D");

	_eval_tmatch_becal = new TTree("tmatch_becal", "BECAL Track Match");
	_eval_tmatch_becal->Branch("deta", &_tm_deta, "_tm_deta/D");
	_eval_tmatch_becal->Branch("dphi", &_tm_dphi, "_tm_dphi/D");
	_eval_tmatch_becal->Branch("dist", &_tm_dist, "_tm_dist/D");
	_eval_tmatch_becal->Branch("eta", &_tm_eta, "_tm_eta/D");
	_eval_tmatch_becal->Branch("phi", &_tm_phi, "_tm_phi/D");
	_eval_tmatch_becal->Branch("p", &_tm_p, "_tm_p/D");
	_eval_tmatch_becal->Branch("ceta", &_tm_ceta, "_tm_ceta/D");
	_eval_tmatch_becal->Branch("cphi", &_tm_cphi, "_tm_cphi/D");
	_eval_tmatch_becal->Branch("q", &_tm_q, "_tm_q/I");
	_eval_tmatch_becal->Branch("e", &_tm_e, "_tm_e/D");

	_eval_tmatch_ihcal = new TTree("tmatch_ihcal", "IHCAL Track Match");
	_eval_tmatch_ihcal->Branch("deta", &_tm_deta, "_tm_deta/D");
	_eval_tmatch_ihcal->Branch("dphi", &_tm_dphi, "_tm_dphi/D");
	_eval_tmatch_ihcal->Branch("dist", &_tm_dist, "_tm_dist/D");
	_eval_tmatch_ihcal->Branch("eta", &_tm_eta, "_tm_eta/D");
	_eval_tmatch_ihcal->Branch("phi", &_tm_phi, "_tm_phi/D");
	_eval_tmatch_ihcal->Branch("p", &_tm_p, "_tm_p/D");
	_eval_tmatch_ihcal->Branch("ceta", &_tm_ceta, "_tm_ceta/D");
	_eval_tmatch_ihcal->Branch("cphi", &_tm_cphi, "_tm_cphi/D");
	_eval_tmatch_ihcal->Branch("q", &_tm_q, "_tm_q/I");
	_eval_tmatch_ihcal->Branch("e", &_tm_e, "_tm_e/D");

	_eval_tmatch_ohcal = new TTree("tmatch_ohcal", "IHCAL Track Match");
	_eval_tmatch_ohcal->Branch("deta", &_tm_deta, "_tm_deta/D");
	_eval_tmatch_ohcal->Branch("dphi", &_tm_dphi, "_tm_dphi/D");
	_eval_tmatch_ohcal->Branch("dist", &_tm_dist, "_tm_dist/D");
	_eval_tmatch_ohcal->Branch("eta", &_tm_eta, "_tm_eta/D");
	_eval_tmatch_ohcal->Branch("phi", &_tm_phi, "_tm_phi/D");
	_eval_tmatch_ohcal->Branch("p", &_tm_p, "_tm_p/D");
	_eval_tmatch_ohcal->Branch("ceta", &_tm_ceta, "_tm_ceta/D");
	_eval_tmatch_ohcal->Branch("cphi", &_tm_cphi, "_tm_cphi/D");
	_eval_tmatch_ohcal->Branch("q", &_tm_q, "_tm_q/I");
	_eval_tmatch_ohcal->Branch("e", &_tm_e, "_tm_e/D");

	_eval_tmatch_femc = new TTree("tmatch_femc", "FEMC Track Match");
	_eval_tmatch_femc->Branch("deta", &_tm_deta, "_tm_deta/D");
	_eval_tmatch_femc->Branch("dphi", &_tm_dphi, "_tm_dphi/D");
	_eval_tmatch_femc->Branch("dist", &_tm_dist, "_tm_dist/D");
	_eval_tmatch_femc->Branch("eta", &_tm_eta, "_tm_eta/D");
	_eval_tmatch_femc->Branch("phi", &_tm_phi, "_tm_phi/D");
	_eval_tmatch_femc->Branch("p", &_tm_p, "_tm_p/D");
	_eval_tmatch_femc->Branch("ceta", &_tm_ceta, "_tm_ceta/D");
	_eval_tmatch_femc->Branch("cphi", &_tm_cphi, "_tm_cphi/D");
	_eval_tmatch_femc->Branch("q", &_tm_q, "_tm_q/I");
	_eval_tmatch_femc->Branch("e", &_tm_e, "_tm_e/D");

	_eval_tmatch_lfhcal = new TTree("tmatch_lfhcal", "LFHCAL Track Match");
	_eval_tmatch_lfhcal->Branch("deta", &_tm_deta, "_tm_deta/D");
	_eval_tmatch_lfhcal->Branch("dphi", &_tm_dphi, "_tm_dphi/D");
	_eval_tmatch_lfhcal->Branch("dist", &_tm_dist, "_tm_dist/D");
	_eval_tmatch_lfhcal->Branch("eta", &_tm_eta, "_tm_eta/D");
	_eval_tmatch_lfhcal->Branch("phi", &_tm_phi, "_tm_phi/D");
	_eval_tmatch_lfhcal->Branch("p", &_tm_p, "_tm_p/D");
	_eval_tmatch_lfhcal->Branch("ceta", &_tm_ceta, "_tm_ceta/D");
	_eval_tmatch_lfhcal->Branch("cphi", &_tm_cphi, "_tm_cphi/D");
	_eval_tmatch_lfhcal->Branch("q", &_tm_q, "_tm_q/I");
	_eval_tmatch_lfhcal->Branch("e", &_tm_e, "_tm_e/D");

	// calotrack cluster matching

	_eval_cmatch_becal_ihcal = new TTree("cmatch_becal_ihcal", "BECAL-IHCAL Neutral Cluster Match");
	_eval_cmatch_becal_ihcal->Branch("deta", &_cm_deta, "_cm_deta/D");
	_eval_cmatch_becal_ihcal->Branch("dphi", &_cm_dphi, "_cm_dphi/D");
	_eval_cmatch_becal_ihcal->Branch("dist", &_cm_dist, "_cm_dist/D");
	_eval_cmatch_becal_ihcal->Branch("eta", &_cm_eta, "_cm_eta/D");
	_eval_cmatch_becal_ihcal->Branch("phi", &_cm_phi, "_cm_phi/D");
	_eval_cmatch_becal_ihcal->Branch("E", &_cm_E, "_cm_E/D");

	_eval_cmatch_becal_ohcal = new TTree("cmatch_becal_ohcal", "BECAL-OHCAL Neutral Cluster Match");
	_eval_cmatch_becal_ohcal->Branch("deta", &_cm_deta, "_cm_deta/D");
	_eval_cmatch_becal_ohcal->Branch("dphi", &_cm_dphi, "_cm_dphi/D");
	_eval_cmatch_becal_ohcal->Branch("dist", &_cm_dist, "_cm_dist/D");
	_eval_cmatch_becal_ohcal->Branch("eta", &_cm_eta, "_cm_eta/D");
	_eval_cmatch_becal_ohcal->Branch("phi", &_cm_phi, "_cm_phi/D");
	_eval_cmatch_becal_ohcal->Branch("E", &_cm_E, "_cm_E/D");

	_eval_cmatch_ihcal_ohcal = new TTree("cmatch_ihcal_ohcal", "IHCAL-OHCAL Neutral Cluster Match");
	_eval_cmatch_ihcal_ohcal->Branch("deta", &_cm_deta, "_cm_deta/D");
	_eval_cmatch_ihcal_ohcal->Branch("dphi", &_cm_dphi, "_cm_dphi/D");
	_eval_cmatch_ihcal_ohcal->Branch("dist", &_cm_dist, "_cm_dist/D");
	_eval_cmatch_ihcal_ohcal->Branch("eta", &_cm_eta, "_cm_eta/D");
	_eval_cmatch_ihcal_ohcal->Branch("phi", &_cm_phi, "_cm_phi/D");
	_eval_cmatch_ihcal_ohcal->Branch("E", &_cm_E, "_cm_E/D");

	_eval_cmatch_femc_lfhcal = new TTree("cmatch_femc_lfhcal", "FEMC-LFHCAL Neutral Cluster Match");
	_eval_cmatch_femc_lfhcal->Branch("deta", &_cm_deta, "_cm_deta/D");
	_eval_cmatch_femc_lfhcal->Branch("dphi", &_cm_dphi, "_cm_dphi/D");
	_eval_cmatch_femc_lfhcal->Branch("dist", &_cm_dist, "_cm_dist/D");
	_eval_cmatch_femc_lfhcal->Branch("eta", &_cm_eta, "_cm_eta/D");
	_eval_cmatch_femc_lfhcal->Branch("phi", &_cm_phi, "_cm_phi/D");
	_eval_cmatch_femc_lfhcal->Branch("E", &_cm_E, "_cm_E/D");	

	// Diagnostic histograms

	_h_nclusters_eemc = new TH1D("_h_nclusters_eemc","",51,-0.5,50.5); 
	_h_nclusters_becal = new TH1D("_h_nclusters_becal","",51,-0.5,50.5); 
	_h_nclusters_ihcal = new TH1D("_h_nclusters_ihcal","",51,-0.5,50.5); 
	_h_nclusters_ohcal = new TH1D("_h_nclusters_ohcal","",51,-0.5,50.5); 
	_h_nclusters_femc = new TH1D("_h_nclusters_femc","",51,-0.5,50.5); 
	_h_nclusters_lfhcal = new TH1D("_h_nclusters_lfhcal","",51,-0.5,50.5); 

	_h_clusteta_eemc = new TH1D("_h_clusteta_eemc","",100,-5.0,5.0); 
	_h_clusteta_becal = new TH1D("_h_clusteta_becal","",100,-5.0,5.0); 
	_h_clusteta_ihcal = new TH1D("_h_clusteta_ihcal","",100,-5.0,5.0); 
	_h_clusteta_ohcal = new TH1D("_h_clusteta_ohcal","",100,-5.0,5.0); 
	_h_clusteta_femc = new TH1D("_h_clusteta_femc","",100,-5.0,5.0); 
	_h_clusteta_lfhcal = new TH1D("_h_clusteta_lfhcal","",100,-5.0,5.0); 

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
	_eval_charged_tracks_bkwd->Write(); 

	_eval_tmatch_eemc->Write(); 
	_eval_tmatch_becal->Write(); 
	_eval_tmatch_ihcal->Write(); 
	_eval_tmatch_ohcal->Write(); 
	_eval_tmatch_femc->Write(); 
	_eval_tmatch_lfhcal->Write(); 

	_eval_cmatch_becal_ihcal->Write(); 
	_eval_cmatch_becal_ohcal->Write(); 
	_eval_cmatch_ihcal_ohcal->Write(); 
	_eval_cmatch_femc_lfhcal->Write(); 

	_eval_calo_tracks_cent->Write(); 
	_eval_calo_tracks_fwd->Write(); 
	_eval_calo_tracks_bkwd->Write(); 

	_h_nclusters_eemc->Write(); 
	_h_nclusters_becal->Write(); 
	_h_nclusters_ihcal->Write(); 
	_h_nclusters_ohcal->Write(); 
	_h_nclusters_femc->Write(); 
	_h_nclusters_lfhcal->Write(); 

	_h_clusteta_eemc->Write(); 
	_h_clusteta_becal->Write(); 
	_h_clusteta_ihcal->Write(); 
	_h_clusteta_ohcal->Write(); 
	_h_clusteta_femc->Write(); 
	_h_clusteta_lfhcal->Write(); 

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
  _hepmcp_W2 = -9999.0;
  _hepmcp_y = -9999.0;  
  _hepmcp_procid = -9999.0; 
  _hepmcp_id1 = -9999; 
  _hepmcp_id2 = -9999; 
	
  e_p_initial = -9999.0; 
  p_p_initial = -9999.0; 

  true_electron_headon = NULL; // not really deadon but in lab frame 
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
 
    // Lab frame!
    //TLorentzVector e_i(0.0,0.0,-e_p_initial,e_p_initial);
    TLorentzVector e_i = true_e_initial; 
    TLorentzVector e_f(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e()); 
    TLorentzVector v_phot = (e_i - e_f); 
    double checkQ2 = -v_phot.Mag2();

    if(fabs(checkQ2-_hepmcp_Q2)<minQ2diff){
      true_electron_headon = g4particle;
      minQ2diff = fabs(checkQ2-_hepmcp_Q2); 
    }

  }

  if(!true_electron_headon){
    LogWarning("scattered electron not found in primaries!"); 

    // Diagnostic for single particle events - cluster/track evaluation

    BuildChargedCaloTracks(topNode,"CENT"); 
    BuildChargedCaloTracks(topNode,"FWD"); 
    BuildChargedCaloTracks(topNode,"BKWD"); 

    // Do the CaloTracks - diagnostic for single particle events
    // note that these will be in the lab frame

    std::vector<fastjet::PseudoJet> tcpseudojets;

    TLorentzRotation breit; // initialized as identity
    TRotation breitRotInv;  // initialized as identity
    
    BuildCaloTracks(topNode, "CENT", tcpseudojets, breit, breitRotInv, "", -1); 
    BuildCaloTracks(topNode, "FWD", tcpseudojets, breit, breitRotInv, "", -1); 
    BuildCaloTracks(topNode, "BKWD", tcpseudojets, breit, breitRotInv, "", -1); 

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
  meas_y = -9999.0; 
  measW2 = -9999.0; 
  meas_E_p = -9999.0;
  electron_eta = -9999.0; 
  electron_phi = -9999.0; 
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
  jet_lab_p.clear(); 
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
  tjet_lab_p.clear(); 
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
  tcjet_lab_p.clear(); 
  tcjet_nc.clear();
  tcjet_pidx.clear(); 
  tcjet_pdR.clear(); 
  tcjet_Q.clear(); 
  tcjet_cf.clear(); 
  tcjet_neut_p.clear(); 
  tcjet_chgd_p.clear(); 
  tcjet_em_p.clear(); 
  tcjet_neut_pm.clear(); 
  tcjet_chgd_pm.clear(); 
  tcjet_em_pm.clear(); 

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
  pjet_neut_pm.clear(); 
  pjet_chgd_pm.clear(); 
  pjet_em_pm.clear(); 
  pjet_lab_eta.clear(); 
  pjet_lab_phi.clear(); 
  pjet_lab_p.clear(); 
  
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
  tfpjet_tcidx.clear(); 
  tfpjet_tcdR.clear(); 
  tfpjet_neut_p.clear(); 
  tfpjet_neut_p.clear(); 
  tfpjet_em_p.clear(); 
  tfpjet_neut_pm.clear(); 
  tfpjet_neut_pm.clear(); 
  tfpjet_em_pm.clear(); 
  tfpjet_lab_eta.clear(); 
  tfpjet_lab_phi.clear(); 
  tfpjet_lab_p.clear(); 

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
	RawCluster *testCluster = getCluster(topNode, temp->get_name(), temp->get_eta(), temp->get_phi(), 
					     -1, electron->get_p(), clustE, clustIdx, clustdR); 
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
    electron_phi = electron->get_phi();
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
    meas_y = (p_initial*virtual_photon)/(p_initial*e_initial);
    measW2 = (p_initial + virtual_photon)*(p_initial + virtual_photon); 

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
      jet_lab_p.push_back(pjet.Vect().Mag()); 

      jet_nc.push_back(fastjets[ijet].constituents().size());
	   
    }

    // Track Jets

    pseudojets.clear(); 
    tcpseudojets.clear(); 

    EICPIDParticleContainer *pidcontainer = findNode::getClass<EICPIDParticleContainer>(topNode, "EICPIDParticleMap");

    for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
	 track_itr != _trackmap->end(); track_itr++) {

      SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

      // Skip scattered electron
      if ((temp->get_truth_track_id() - true_electron_headon->get_track_id()) == 0) continue; 

      double px = temp->get_px();
      double py = temp->get_py();
      double pz = temp->get_pz();

      double mass = 0.13957039; // assume pion

      // Check the particle ID hypothesis

      if (pidcontainer && USE_PID ){

	// EICPIDParticle are index the same as the tracks
	const EICPIDParticle *pid_particle = pidcontainer->findEICPIDParticle(temp->get_id());

	if (pid_particle)
	  {
	    // top level log likelihood sums.
	    // More detailed per-detector information also available at  EICPIDParticle::get_LogLikelyhood(EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector)
	    float m_tr_electron_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::ElectronCandiate);
	    float m_tr_pion_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::PionCandiate);
	    float m_tr_kaon_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::KaonCandiate);
	    float m_tr_proton_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::ProtonCandiate);

	    if((m_tr_electron_loglikelihood-m_tr_pion_loglikelihood)>0.0) mass = 0.000510998950; 
	    if((m_tr_kaon_loglikelihood-m_tr_pion_loglikelihood)>0.0) {
	      mass = 0.493677; 
	      if((m_tr_proton_loglikelihood-m_tr_kaon_loglikelihood)>0.0) mass = 0.938272; 
	    }

	  }

      } 
      else{
	if(USE_PID) cout << "EICPIDParticle missing! Assuming massless." << endl; 
	mass = 0.0; 
      }

      // Transform to Breit frame
      TLorentzVector track(px,py,pz,sqrt(pow(temp->get_p(),2) + pow(mass,2))); 
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
      tjet_lab_p.push_back(pjet.Vect().Mag()); 

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
    
#ifdef PARAM_PHOTONS
    BuildParametrizedPhotonCaloTracks(topNode, "CENT", tcpseudojets, breit, breitRotInv); 
    BuildParametrizedPhotonCaloTracks(topNode, "FWD", tcpseudojets, breit, breitRotInv); 
    BuildParametrizedPhotonCaloTracks(topNode, "BKWD", tcpseudojets, breit, breitRotInv); 
#else
    BuildCaloTracks(topNode, "CENT", tcpseudojets, breit, breitRotInv, ECDetName, ECIdx); 
    BuildCaloTracks(topNode, "FWD", tcpseudojets, breit, breitRotInv, ECDetName, ECIdx); 
    BuildCaloTracks(topNode, "BKWD", tcpseudojets, breit, breitRotInv, ECDetName, ECIdx); 
#endif


#ifdef PARAM_HAD_NEUTRALS
    BuildParametrizedHadCaloTracks(topNode, "CENT", tcpseudojets, breit, breitRotInv); 
    BuildParametrizedHadCaloTracks(topNode, "FWD", tcpseudojets, breit, breitRotInv); 
    BuildParametrizedHadCaloTracks(topNode, "BKWD", tcpseudojets, breit, breitRotInv); 
#endif

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

      tcjet_nc.push_back(tcfastjets[ijet].constituents().size());
	 
      // Jet Charge
      std::vector<fastjet::PseudoJet> tconstit = tcfastjets[ijet].constituents();
	    
      double jptot = sqrt(pow(tcfastjets[ijet].px(),2) + 
			  pow(tcfastjets[ijet].py(),2) + 
			  pow(tcfastjets[ijet].pz(),2)); 

      tcjet_Q.push_back(JetCharge(&tconstit,jptot));
      tcjet_cf.push_back(JetChargedFraction(&tconstit,pjet.Vect()));

      tcjet_neut_p.push_back(JetNeutralMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 
      tcjet_chgd_p.push_back(JetChargedMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 
      tcjet_em_p.push_back(JetEMMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 

      tcjet_neut_pm.push_back(JetNeutralMomentum(&tconstit).Mag()); 
      tcjet_chgd_pm.push_back(JetChargedMomentum(&tconstit).Mag()); 
      tcjet_em_pm.push_back(JetEMMomentum(&tconstit).Mag()); 

      // Transform back to lab frame
      pjet.Transform(breitRot); 
      pjet = (breitInv * pjet); 

      tcjet_lab_eta.push_back(pjet.Eta());  
      tcjet_lab_phi.push_back(pjet.Phi()); 
      tcjet_lab_p.push_back(pjet.Vect().Mag()); 

    }

    // Primary Jets

    std::vector<fastjet::PseudoJet> pfastjets; 
    std::vector<fastjet::PseudoJet> tfpfastjets; 

    // Primary Jets (in the experiment Breit frame)
    GetPrimaryJets(topNode, jetdef, &pfastjets, breit, breitRotInv, breitInv, breitRot, p_initial_breit, virtual_photon_breit, false); 

    // Primary Jets in "true" Breit frame
    // Set up the transformation (boost) from the lab to the "true" Breit frame
    // Transformation based on truth kinematics
 
    // Get the final electron in the lab frame
    CLHEP::HepLorentzVector efp(true_electron_headon->get_px(),true_electron_headon->get_py(),
				true_electron_headon->get_pz(),true_electron_headon->get_e());
    // Already in lab frame! 
    //efp = EventToLab * efp;  
    TLorentzVector true_e_final(efp.px(), efp.py(), efp.pz(), efp.e()); 
    TLorentzVector true_virtual_photon = (true_e_initial - true_e_final); 

    _hepmcp_W2 = (true_p_initial + true_virtual_photon)*(true_p_initial + true_virtual_photon); 
    _hepmcp_y = (true_p_initial*true_virtual_photon)/(true_p_initial*true_e_initial); 

    // Recalculate the kinematics 
    // (clean up round-off errors and keep everything self-consistent)
    double _tf_Q2 = -true_virtual_photon.Mag2(); 
    double _tf_x2 = _tf_Q2/(2*true_virtual_photon*true_p_initial);

    TVector3 P3t = true_p_initial.Vect(); 
    TVector3 q3t = true_virtual_photon.Vect(); 
    TVector3 boost_tf = -(2.0*_tf_x2*P3t + q3t)*(1.0/(2.0*_tf_x2*true_p_initial.E() + true_virtual_photon.E())); 
    TLorentzRotation breit_tf = TLorentzRotation().Boost(boost_tf); 
    TLorentzRotation breitInv_tf = TLorentzRotation().Boost(-boost_tf); 

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

    GetPrimaryJets(topNode, jetdef, &tfpfastjets, breit_tf, breitRotInv_tf, breitInv_tf, breitRot_tf, 
		   p_initial_breit_tf, virtual_photon_breit_tf, true); 

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


    for(unsigned int pjet=0; pjet<tfpfastjets.size(); pjet++){


      int closest_idx = -1; 
      double closest_dR = 9999.0; 

      for(unsigned int ijet=0; ijet<tcfastjets.size(); ijet++){

	TVector3 v1(tcfastjets[ijet].px(),tcfastjets[ijet].py(),tcfastjets[ijet].pz()); 
	TVector3 v2(tfpfastjets[pjet].px(),tfpfastjets[pjet].py(),tfpfastjets[pjet].pz());

	double tc_p = sqrt(pow(tcfastjets[ijet].px(),2) + 
			   pow(tcfastjets[ijet].py(),2) +
			   pow(tcfastjets[ijet].pz(),2));

	double p_p = sqrt(pow(tfpfastjets[pjet].px(),2) + 
			  pow(tfpfastjets[pjet].py(),2) + 
			  pow(tfpfastjets[pjet].pz(),2)); 

	double dR = acos(v1.Dot(v2)/(tc_p*p_p)); 

	if(dR<closest_dR){
	  closest_idx = pjet; 
	  closest_dR = dR; 
	}

      }

      tfpjet_tcidx.push_back(closest_idx); 
      tfpjet_tcdR.push_back(closest_dR); 

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

bool CentauroJets::PassClusterEtaCut(double eta, std::string detName ){
  
  // Cuts are *after* clusters are adjusted to track match locations

  if(detName=="BECAL") {
    if( (eta>-1.8) && (eta<1.3) ) 
      return true; 
    else
      return false; 
  }

  if(detName=="EEMC") {
    if( (eta>-3.4) && (eta<-1.8) ) 
      return true; 
    else
      return false; 
  }

  if(detName=="FEMC") {
    if( (eta>1.3) && (eta<3.4) ) 
      return true; 
    else
      return false; 
  }

  if(detName=="HCALIN") {
    if( (eta>-1.6) && (eta<1.2) ) 
      return true; 
    else
      return false; 
  }

  if(detName=="HCALOUT") {
    if( (eta>-1.2) && (eta<1.2) ) 
      return true; 
    else
      return false; 
  }

  if(detName=="LFHCAL") {
    if( (eta>1.2) && (eta<3.4) ) 
      return true; 
    else
      return false; 
  }

  return false; 

}

void CentauroJets::GetClusterOffsetFit( double &eta, double &phi, int charge, double p, std::string detName ){

  double new_phi_pos = phi; 
  double new_eta_pos = eta; 
  double new_phi_neg = phi; 
  double new_eta_neg = eta; 

  if(detName=="BECAL") {

    if(charge<=0){

      if(p<0.75){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.360005 + 0.0119768*eta + -0.191993*eta*eta) ;
	else if(eta<-1.0)
	  new_phi_neg = phi - (-1.18532 + -0.624726*eta);
	else
	  new_phi_neg = phi - (-0.669893 + 0.138791*eta);

	new_eta_neg = -0.0432953 + 1.0674*eta; 

      }
      else if(p>=0.75 && p<1.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.268165 + 0.0141392*eta + -0.157004*eta*eta) ;
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.865813 + -0.438634*eta);
	else
	  new_phi_neg = phi - (-0.469131 + 0.0846924*eta);

	new_eta_neg = -0.052961 + 1.04777*eta; 

      }
      else if(p>=1.0 && p<2.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.167403 + 0.0119659*eta + -0.108719*eta*eta) ;
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.575884 + -0.292328*eta);
	else
	  new_phi_neg = phi - (-0.407185 + 0.138671*eta);

	new_eta_neg = -0.0523756 + 1.05244*eta; 

      }
      else if(p>=2.0 && p<3.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.0987133 + 0.00757871*eta + -0.0612283*eta*eta) ;
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.353617 + -0.182957*eta);
	else
	  new_phi_neg = phi - ( -0.218702 +  0.0532541*eta);

	new_eta_neg = -0.0513668 + 1.05107*eta;

      }
      else if(p>=3.0 && p<4.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.0653631 + 0.0031072*eta + -0.0466652*eta*eta) ;
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.231377 + -0.115999*eta);
	else
	  new_phi_neg = phi - (-0.200313 +  0.0868494*eta);

	new_eta_neg = -0.0508864 + 1.04994*eta;

      }
      else if(p>=4.0 && p<5.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.0512339 +  0.0046292*eta + -0.0376337*eta*eta);
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.161873 + -0.0788219*eta);
	else
	  new_phi_neg = phi - (-0.153439  + 0.0676106*eta);

	new_eta_neg = -0.0500371 + 1.04697*eta;
	
      }
      else if(p>=5.0 && p<10.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.0313837 +  0.00193104*eta + -0.021448*eta*eta);
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.109924 + -0.0559835*eta);
	else
	  new_phi_neg = phi - (-0.103726 + 0.0493075*eta);

	new_eta_neg = -0.0494103 + 1.04727*eta;

      }
      else if(p>=10.0 && p<15.0){

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.0178984 + 0.00114097*eta + -0.0113011*eta*eta);
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.0623192 + -0.0318633*eta);
	else
	  new_phi_neg = phi - (-0.0419758 +0.0127882*eta);

	new_eta_neg = -0.0484053 + 1.04748*eta;

      }
      else{

	if(abs(eta)<=1.0)
	  new_phi_neg = phi - (-0.0125593 + 0.000676921*eta + -0.00685772*eta*eta);
	else if(eta<-1.0)
	  new_phi_neg = phi - (-0.0423715 + -0.0214144*eta);
	else
	  new_phi_neg = phi - (-0.0288464 + 0.00890873*eta);

	new_eta_neg = -0.0474674 + 1.04464*eta;

      }

    }
    if(charge>=0){

      if(p<0.75){

	if(eta<=-1.0){
	  new_phi_pos = phi - (1.12761 + 0.565174*eta);
	  new_eta_pos = -0.0443036 + 1.09318*eta; 
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.356764 +  -0.00946361*eta +  0.190724*eta*eta);
	  new_eta_pos = -0.056254 + 1.05487*eta; 
	}
	else{
	  new_phi_pos = phi - (0.614763 + -0.0985048*eta);
	  new_eta_pos = -0.505265 + 1.47364*eta; 
	}

      }
      else if(p>=0.75 && p<1.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (1.01863 + 0.528971*eta) ;
	  new_eta_pos = -0.0468911 + 1.07121*eta; 
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.299654 + -0.0121693*eta + 0.17995*eta*eta);
	  new_eta_pos = -0.0541906 + 1.04853*eta; 
	}
	else{
	  new_phi_pos = phi - (0.693192 + -0.252914*eta);
	  new_eta_pos = -0.270632 + 1.25931*eta; 
	}


      }
      else if(p>=1.0 && p<2.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.592628 + 0.299501*eta) ;
	  new_eta_pos = -0.0452062+ 1.05673*eta; 
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.17873 + -0.00604416*eta + 0.112273*eta*eta);
	  new_eta_pos = -0.0516357 + 1.05131*eta; 
	}
	else{
	  new_phi_pos = phi - (0.200865 + 0.0559365*eta);
	  new_eta_pos = -0.111421 + 1.11325*eta; 
	}

      }
      else if(p>=2.0 && p<3.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.38207 + 0.194595*eta) ;
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.108245 + -0.00600031*eta + 0.0685143*eta*eta);
	}
	else{
	  new_phi_pos = phi - (0.188961 + -0.015496*eta);
	}

	new_eta_pos = -0.05072 + 1.04944*eta; 

      }
      else if(p>=3.0 && p<4.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.266011 + 0.133414*eta) ;
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.0781154 + -0.00342547*eta + 0.0534418*eta*eta);
	}
	else{
	  new_phi_pos = phi - ( 0.115435 + 0.00274316*eta);
	}

	new_eta_pos = -0.0491458 + 1.04895*eta; 

      }
      else if(p>=4.0 && p<5.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.193749 + 0.0914972*eta) ;
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.0611015 + -0.00434471*eta + 0.043988*eta*eta);
	}
	else{
	  new_phi_pos = phi - (0.166096 + -0.0617353*eta);
	}

	new_eta_pos = -0.0617353 + 1.04652*eta; 
	
      }
      else if(p>=5.0 && p<10.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.135296 + 0.0675732*eta) ;
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.0375606 + -0.00284185*eta + 0.0290346*eta*eta);
	}
	else{
	  new_phi_pos = phi - (0.101738 + -0.0357625*eta);
	}

	new_eta_pos = -0.0477418 + 1.04549*eta; 

      }
      else if(p>=10.0 && p<15.0){

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.0726225 + 0.0342078*eta) ;
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.0207465 + -0.00137859*eta + 0.0177749*eta*eta);
	}
	else{
	  new_phi_pos = phi - (0.0345588);
	}

	new_eta_pos = -0.0471691 + 1.04462*eta; 

      }
      else{

	if(eta<=-1.0){
	  new_phi_pos = phi - (0.0551176 + 0.0262792*eta) ;
	}
	else if(eta>-1.0 && eta<=1.0){
	  new_phi_pos = phi - (0.0156055+ -0.00102952*eta + 0.0113895*eta*eta);
	}
	else{
	  new_phi_pos = phi - ( 0.036707 + -0.0102788*eta);
	}

	new_eta_pos =  -0.0469741 + 1.04337*eta; 

      }

    }

  }

  if(detName=="EEMC") {

    if(charge<=0){

      if(p<0.75){
	new_phi_neg = phi - ( -0.375567 + 0.000602845*eta); 
	new_eta_neg = 0.00264235 + 0.994413*eta; 	
      }
      else if(p>=0.75 && p<1.0){
	new_phi_neg = phi - ( -0.31117 +  -0.00976997*eta); 
	new_eta_neg = 0.0230682 + 1.00216*eta; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_neg = phi - ( -0.156496 + 0.00188529*eta ); 
	new_eta_neg =  0.0318744 + 1.00521*eta; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_neg = phi - ( -0.0980362 + -0.000439882*eta ); 
	new_eta_neg = 0.0323534 + 1.00485*eta;
      }
      else if(p>=3.0 && p<4.0){
	new_phi_neg = phi - (-0.0707824 + -0.000495876*eta);
	new_eta_neg = 0.0219295 + 1.00167*eta;
      }
      else if(p>=4.0 && p<5.0){
	new_phi_neg = phi - (-0.0504497 + 0.0013824*eta );
	new_eta_neg = 0.0208685 + 1.00119*eta;
      }
      else if(p>=5.0 && p<10.0){
	new_phi_neg = phi - (-0.0346814 +  -0.000599061*eta);
	new_eta_neg = 0.0206542 + 1.00283*eta;
      }
      else if(p>=10.0 && p<15.0){
	new_phi_neg = phi - ( -0.0156642 + 0.00167558*eta);
	new_eta_neg = 0.0239701 + 1.00501*eta;
      }
      else{
	new_phi_neg = phi - (-0.0137015 + 0.000308886*eta); 
	new_eta_neg = 0.0178779 + 1.00294*eta; 
      }

    }
    if(charge>=0){

      if(p<0.75){
	new_phi_pos = phi - (0.405946 + 0.00963873*eta); 
	new_eta_pos = -0.000140854 + 0.993976*eta; 	
      }
      else if(p>=0.75 && p<1.0){
	new_phi_pos = phi - ( 0.298192 + 0.00563709*eta );  
	new_eta_pos = 0.00695631 + 0.99481*eta; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_pos = phi - (0.149684 + -0.00442267*eta);
	new_eta_pos = 0.0287754 + 1.00318*eta; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_pos = phi - (0.097814 + 0.000770295*eta); 
	new_eta_pos =  0.0311729 + 1.00486*eta;
      }
      else if(p>=3.0 && p<4.0){
	new_phi_pos = phi - (0.0705052 + 0.000872545*eta);
	new_eta_pos =  0.0210661 + 1.00137*eta;
      }
      else if(p>=4.0 && p<5.0){
	new_phi_pos = phi - (0.0543699 + -0.00019804*eta);
	new_eta_pos = 0.0217715 + 1.00208*eta;
      }
      else if(p>=5.0 && p<10.0){
	new_phi_pos = phi - (0.0326268 + -0.000204726*eta);
	new_eta_pos =  0.0256571 +  1.00449*eta;
      }
      else if(p>=10.0 && p<15.0){
	new_phi_pos = phi - ( 0.0167583 + -0.00117102*eta);
	new_eta_pos = 0.0203821 + 1.00359*eta;
      }
      else{
	new_phi_pos = phi - (0.0133368 + -0.000432056*eta); 
	new_eta_pos = 0.0211667 + 1.00413*eta; 
      }

    }
    
  }

  if(detName=="FEMC") {

    if(charge<=0){

      if(p<0.75){
	new_phi_neg = phi - (0.273612 + -0.420079*eta + 0.077171*eta*eta) - (-0.222667 + 0.212773*eta + -0.0366431*eta*eta); 
	new_eta_neg =  0.12501 + 0.980493*eta; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_neg = phi - (0.0105355 + -0.108362*eta +  0.0190665*eta*eta) - (0.00877707 + -0.00479967*eta); 
	new_eta_neg = 0.0414085 + 0.987142*eta;
     }
      else if(p>=1.0 && p<2.0){
	new_phi_neg = phi - (0.0202877 + -0.0732956*eta + 0.0121526*eta*eta); 
	new_eta_neg = -0.0163598 + 1.00125*eta; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_neg = phi - (-0.00392383 + -0.00392383*eta + 0.00601271*eta*eta) - (0.000861668 + -0.0225117*eta); 
	new_eta_neg = -0.0228555 + 1.00222*eta;
      }
      else if(p>=3.0 && p<4.0){
	new_phi_neg = phi - (0.00567538 + -0.0320377*eta + 0.00589812*eta*eta);
	new_eta_neg = -0.0214558 + 1.00145*eta;
     }
      else if(p>=4.0 && p<5.0){
	new_phi_neg = phi - (-0.00832566 + -0.0147959*eta + 0.00260797*eta*eta);
	new_eta_neg = -0.0218401 + 1.00223*eta; 
     }
      else if(p>=5.0 && p<10.0){
	new_phi_neg = phi - (-0.0076235 + -0.00697934*eta +  0.00108543*eta*eta);
	new_eta_neg =  -0.020734 + 1.00258*eta;
      }
      else if(p>=10.0 && p<15.0){
	new_phi_neg = phi - (-0.00273903 + -0.00566229*eta +  0.00098751*eta*eta);
	new_eta_neg = -0.017006 + 1.00236*eta;
      }
      else{
	new_phi_neg = phi - (-0.00328948 +  -0.00317229*eta + 0.000517406*eta*eta);
	new_eta_neg =  -0.017006 + 1.00236*eta;
      }

    }
    if(charge>=0){

      if(p<0.75){
	new_phi_pos = phi - (-0.0212205 + 0.157972*eta + -0.0268174*eta*eta); 
	new_eta_pos = 0.045431 + 0.992563*eta; 	
       }
      else if(p>=0.75 && p<1.0){
	new_phi_pos = phi - ( -3.48109 + 5.71202*eta + -3.32969*eta*eta + 0.852357*pow(eta,3) + -0.0807565*pow(eta,4)); 
	new_eta_pos = 0.0145786 + 0.996376*eta; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_pos = phi - (0.0454196 + 0.0233758*eta +  -0.00283064*eta*eta); 
	new_eta_pos = -0.0111863 + 0.998818*eta; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_pos = phi - (0.0101583 + 0.0295474*eta + -0.00507624*eta*eta); 
	new_eta_pos =  -0.0211347 + 1.00097*eta;
      }
      else if(p>=3.0 && p<4.0){
	new_phi_pos = phi - (0.0071257 +  0.0206688*eta + -0.00338545*eta*eta);
	new_eta_pos = -0.0205476 + 1.00104*eta;
      }
      else if(p>=4.0 && p<5.0){
	new_phi_pos = phi - (0.00891172 + 0.0144097*eta + -0.00239306*eta*eta);
	new_eta_pos = -0.0219039 + 1.00205*eta;
      }
      else if(p>=5.0 && p<10.0){
	new_phi_pos = phi - (0.00724539 + 0.006736*eta + -0.00094245*eta*eta);
	new_eta_pos = -0.0205005 + 1.00253*eta;
      }
      else if(p>=10.0 && p<15.0){
	new_phi_pos = phi - (0.00276497 + 0.00590355*eta + -0.00109209*eta*eta);
	new_eta_pos = -0.019365 + 1.00317*eta;
      }
      else{
	new_phi_pos = phi - (0.0020188 + 0.00402878*eta + -0.000666731*eta*eta); 
	new_eta_pos = -0.0169182 + 1.00251*eta; 
       }

    }


  }

  if(detName=="HCALIN") {

    if(charge<=0){

      if(p<0.75){
	new_phi_neg = phi - (-4.74674e-01); 
	new_eta_neg = 0.00715007 + 0.871202*eta + -0.0244734*eta*eta + + 0.164271*eta*eta*eta;  	
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_neg = phi - (-0.346706 +  0.00827334*eta); 
	new_eta_neg = -0.00172789 + 0.965483*eta + -0.00384477*eta*eta + 0.0458575*eta*eta*eta; 
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_neg = phi - ( -0.207237 + 0.0441341*eta); 
	new_eta_neg = 0.00549633 + 1.02316*eta + -0.0170673*eta*eta + -0.00350476*eta*eta*eta; 
	// phase2
	new_phi_neg -= 1.60553e-02; 
	new_eta_neg -= -2.10639e-02; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_neg = phi - (-0.127344 + -0.00619482*eta); 
	new_eta_neg = -0.013462 + 0.919855*eta + 0.0208697*eta*eta + 0.114938*eta*eta*eta;
	// phase2
	new_phi_neg -= 2.64288e-02; 
	new_eta_neg -= -2.70426e-02; 
      }
      else if(p>=3.0 && p<4.0){
	new_phi_neg = phi - (-0.0902535 + -0.0321342*eta);
	new_eta_neg = -0.00697089 + 0.990275*eta;
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
       }
      else if(p>=4.0 && p<5.0){
	new_phi_neg = phi - (-0.0661588 + -0.00668812*eta);
	new_eta_neg = -0.0144468 + 1.02263*eta;
	// phase2
	//new_phi_neg -= ; 
	new_eta_neg -= -1.74516e-02; 
      }
      else if(p>=5.0 && p<10.0){
	new_phi_neg = phi - (-0.0400165 + -0.0134691*eta);
	new_eta_neg = -0.0108213 + 1.00154*eta;
	// phase2
	//new_phi_neg -= ; 
	new_eta_neg -= -7.86095e-03; 
      }
      else if(p>=10.0 && p<15.0){
	new_phi_neg = phi - (-0.0201962 + -0.00145956*eta);
	new_eta_neg =  0.00582504 + 1.00652*eta;
	// phase2
	new_phi_neg -= 4.32095e-03; 
	new_eta_neg -= 2.72789e-03; 
      }
      else{
	new_phi_neg = phi - (-0.0179695 + (6.97606e-05)*eta); 
	new_eta_neg = -0.00109953 + 1.00648*eta; 
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }


    }
    if(charge>=0){

      if(p<0.75){
	new_phi_pos = phi - (0.423163 +  0.00299655*eta); 
	new_eta_pos = 0.0174807 + 0.724978*eta + 0.0232405*eta*eta +  0.31458*eta*eta*eta; 	
	// phase2
	new_phi_pos -= 0.0396; 
	new_eta_pos -= 0.01262; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_pos = phi - (0.357757 + 0.0483568*eta + 0.0358377*eta*eta); 
	new_eta_pos = 0.023188 + 1.00114*eta; 
	// phase2
	new_phi_pos -= -8.07508e-03; 
	new_eta_pos -= 3.25449e-02; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_pos = phi - (0.25608 + 0.00654453*eta); 
	new_eta_pos = -0.008085 + 1.03385*eta; 
	// phase2
	new_phi_pos -= -3.70456e-02; 
	new_eta_pos -=  7.88243e-03; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_pos = phi - (0.145427 + 0.0134292*eta); 
	new_eta_pos =  0.000756643 + 1.01152*eta;
	// phase2
	new_phi_pos -= -3.23846e-02; 
	new_eta_pos -= 9.24605e-03; 
      }
      else if(p>=3.0 && p<4.0){
	new_phi_pos = phi - (0.0992292);
	new_eta_pos = -0.0105483 + 1.02258*eta;
	// phase2
	new_phi_pos -= -6.00624e-03; 
	new_eta_pos -= -1.54023e-02; 
      }
      else if(p>=4.0 && p<5.0){
	new_phi_pos = phi - (0.10065);
	new_eta_pos = -0.0143708 + 0.987576*eta;
	// phase2
	new_phi_pos -= -3.34988e-02; 
	new_eta_pos -= -7.50205e-03; 
      }
      else if(p>=5.0 && p<10.0){
	new_phi_pos = phi - (0.0444226);
	new_eta_pos = -0.00545667 + 1.00571*eta;
	// phase2
	new_phi_pos -=  -2.84635e-03; 
	new_eta_pos -= -4.32661e-03; 
      }
      else if(p>=10.0 && p<15.0){
	new_phi_pos = phi - (0.0263697);
	new_eta_pos = -0.00286946 + 1.02353*eta;
	// phase2
	new_phi_pos -= 1.66569e-03; 
	new_eta_pos -= -8.34484e-03; 
      }
      else{
	new_phi_pos = phi - (0.021827); 
	new_eta_pos = -0.00387954 + 1.01595*eta; 
	// phase2
	new_phi_pos -= 8.93771e-04; 
	new_eta_pos -= 4.79322e-04; 
      }

    }

  }


  if(detName=="HCALOUT") {

    if(charge<=0){

      if(p<0.75){
	new_phi_neg = phi - (-0.376324); 
	new_eta_neg = 0.0343632 + 1.01011*eta; 	
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_neg = phi - (-0.240655); 
	new_eta_neg = 0.0340835 + 0.994548*eta; 
	// phase2
	new_phi_neg -= -7.44996e-03; 
	new_eta_neg -= 8.32341e-03; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_neg = phi - (-0.201036); 
	new_eta_neg = 0.00732727 + 1.10957*eta; 
	// phase2
	new_phi_neg -= 2.64288e-02; 
	new_eta_neg -= -2.10639e-02; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_neg = phi - (-0.150696 + 0.0234669*eta); 
	new_eta_neg = 0.0143835+  0.963022*eta + -0.0429179*eta*eta + 0.0943905*eta*eta*eta;
	// phase2
	//new_phi_neg -= ; 
	new_eta_neg -= -2.70426e-02; 
      }
      else if(p>=3.0 && p<4.0){
	new_phi_neg = phi - (-1.13335e-01);
	new_eta_neg = -0.000437473 + 1.02849*eta;
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }
      else if(p>=4.0 && p<5.0){
	new_phi_neg = phi - (-9.50999e-02);
	new_eta_neg = 0.000166058 + 1.03416*eta;
	// phase2
	//new_phi_neg -= ; 
	new_eta_neg -= -1.74516e-02; 
      }
      else if(p>=5.0 && p<10.0){
	new_phi_neg = phi - (-0.0790791);
	new_eta_neg = -0.00500539 + 1.04528*eta;
	// phase2
	//new_phi_neg -= ; 
	new_eta_neg -= -7.86095e-03; 
      }
      else if(p>=10.0 && p<15.0){
	new_phi_neg = phi - (-0.0684646);
	new_eta_neg = 0.000269227 + 1.04872*eta;
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }
      else{
	new_phi_neg = phi - (-0.0597346); 
	new_eta_neg = -0.00761354 + 1.0511*eta; 
	// phase2
	//new_phi_neg -= ; 
	//new_eta_neg -= ; 
      }

    }
    if(charge>=0){

      if(p<0.75){
	new_phi_pos = phi - ( 0.269714 + 0.0150204*eta + 0.113267*eta*eta +  0.375956*eta*eta*eta); 
	new_eta_pos = 0.0226544 + 0.783425*eta + -0.130746*eta*eta + 0.555223*eta*eta*eta; 	
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -= 0.02308; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_pos = phi - (0.168688 + 0.0551862*eta + -0.0951666*eta*eta); 
	new_eta_pos =  0.00277313 +  0.849803*eta + 0.0637895*eta*eta + 0.420371*eta*eta*eta; 
	// phase2
	new_phi_pos -= 0.01578; 
	//new_eta_pos -= ; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_pos = phi - (0.0763659); 
	new_eta_pos = 0.00467778 + 0.925522*eta + -0.0168405*eta*eta + 0.234916*eta*eta*eta; 
	// phase2
	new_phi_pos -= 1.80147e-02; 
	//new_eta_pos -= ; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_pos = phi - (0.0559355); 
	new_eta_pos = -0.00734858 +  1.01537*eta;
	// phase2
	new_phi_pos -= -3.64490e-03; 
	new_eta_pos -= 3.84732e-03; 
      }
      else if(p>=3.0 && p<4.0){
	new_phi_pos = phi - (2.39834e-02);
	new_eta_pos = -0.00531465 + 1.02202*eta;
	// phase2
	new_phi_pos -= 9.82503e-03; 
	new_eta_pos -= -8.20733e-03; 
      }
      else if(p>=4.0 && p<5.0){
	new_phi_pos = phi - (9.19719e-03);
	new_eta_pos = -0.00631507 + 1.03453*eta;
	// phase2
	new_phi_pos -= 5.45259e-03; 
	new_eta_pos -= -6.58406e-03; 
      }
      else if(p>=5.0 && p<10.0){
	new_phi_pos = phi - (-0.00893872);
	new_eta_pos = 0.000622195 + 1.04355*eta;
	// phase2
	new_phi_pos -= 8.17611e-03; 
	new_eta_pos -= 3.30926e-03; 
      }
      else if(p>=10.0 && p<15.0){
	new_phi_pos = phi - (-0.0229137);
	new_eta_pos =  -0.00528808 + 1.05282*eta;
	// phase2
	new_phi_pos -= 4.86959e-03; 
	new_eta_pos -= -4.14862e-03; 
      }
      else{
	new_phi_pos = phi - (-3.20578e-02); 
	new_eta_pos = -0.00344713 + 1.05322*eta; 
	// phase2
	new_phi_pos -= 1.53409e-02; 
	new_eta_pos -= -2.40236e-03; 
      }

    }

  }


  if(detName=="LFHCAL") {

    if(charge<=0){

      if(p<0.75){
	new_phi_neg = phi - (-0.0686378 + -0.0364384*eta); 
	new_eta_neg = 3.56676 + -4.05714*eta + 2.33484*eta*eta + -0.340504*eta*eta*eta; 	
	// phase2
	new_phi_neg -= -4.32383e-02; 
	new_eta_neg -= 2.53794e-02; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_neg = phi - (-0.0430947 + -0.0391439*eta); 
	new_eta_neg =  1.55612 + -1.02384*eta + 0.861299*eta*eta + -0.117237*eta*eta*eta; 
	// phase2
	new_phi_neg -= 4.33318e-03; 
	new_eta_neg -= 1.22613e-02; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_neg = phi - (-0.0717001 + (1.92653e-05)*eta); 
	new_eta_neg = 0.570644 + 0.0162908*eta + 0.545085*eta*eta + -0.0926713*eta*eta*eta; 
	// phase2
	new_phi_neg -= -3.46997e-03; 
	new_eta_neg -= 1.24461e-02; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_neg = phi - (-0.0489048 + 0.0040359*eta); 
	new_eta_neg =  1.55765 + -1.34762*eta + 1.15313*eta*eta + -0.179604*eta*eta*eta;
	// phase2
	new_phi_neg -= -1.13370e-03; 
	new_eta_neg -= 5.17301e-02; 
      }
      else if(p>=3.0 && p<4.0){
	new_phi_neg = phi - (0.00815418 + -0.0123175*eta);
	new_eta_neg = 1.7756 + -1.69371*eta + 1.3308*eta*eta + -0.207112*eta*eta*eta;
	// phase2
	new_phi_neg -= -4.33432e-03; 
	new_eta_neg -= 4.66335e-02; 
      }
      else if(p>=4.0 && p<5.0){
	new_phi_neg = phi - (0.00530316 + -0.0160722*eta);
	new_eta_neg = 1.51211 + -1.19703*eta + 1.0704*eta*eta + -0.166915*eta*eta*eta;
	// phase2
	new_phi_neg -= -4.33432e-03; 
	new_eta_neg -= 7.19574e-02; 
      }
      else if(p>=5.0 && p<10.0){
	new_phi_neg = phi - (-0.00796392 + -0.000177456*eta);
	new_eta_neg = 0.501163 + 0.0465134*eta +  0.580546*eta*eta + -0.10445*eta*eta*eta;
	// phase2
	new_phi_neg -= -4.40678e-03; 
	new_eta_neg -= 7.27446e-02; 
      }
      else if(p>=10.0 && p<15.0){
	new_phi_neg = phi - (-0.0212337 + 0.011897*eta);
	new_eta_neg =  0.551743 + -0.0508371*eta + 0.643217*eta*eta + -0.116764*eta*eta*eta;
	// phase2
	new_phi_neg -= -6.77285e-03; 
	new_eta_neg -= 7.08640e-02; 
      }
      else{
	new_phi_neg = phi - (0.0122706 + -0.00912465*eta); 
	new_eta_neg =  0.284003 + 0.299047*eta + 0.507267*eta*eta + -0.100554*eta*eta*eta; 
	// phase2
	new_phi_neg -= 3.13490e-03; 
	new_eta_neg -= 7.72306e-02; 
      }

    }
    if(charge>=0){

      if(p<0.75){
	new_phi_pos = phi - (0.0561723 + 0.0557544*eta); 
	new_eta_pos =  3.26025 + -3.22304*eta + 1.80164*eta*eta + -0.247645*eta*eta*eta; 	
	// phase2
	new_phi_pos -= 2.01034e-02; 
	new_eta_pos -= -2.54962e-03; 
      }
      else if(p>=0.75 && p<1.0){
	new_phi_pos = phi - (0.0484116 +  0.0384422*eta); 
	new_eta_pos = 1.40153 + -0.921705*eta + 0.881702*eta*eta + -0.130419*eta*eta*eta; 
	// phase2
	new_phi_pos -= -0.007605; 
	new_eta_pos -= 2.09980e-02; 
      }
      else if(p>=1.0 && p<2.0){
	new_phi_pos = phi - (0.05913 + 0.00231293*eta); 
	new_eta_pos = 1.4977 + -1.205*eta + 1.04887*eta*eta + -0.157632*eta*eta*eta; 
	// phase2
	new_phi_pos -= 2.12404e-03; 
	new_eta_pos -= 3.46216e-02; 
      }
      else if(p>=2.0 && p<3.0){
	new_phi_pos = phi - (0.0373754 + 0.00455856*eta); 
	new_eta_pos =  1.26274 + -0.86757*eta +  0.908855*eta*eta + -0.140645*eta*eta*eta;
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -=  4.06292e-02; 
      }
      else if(p>=3.0 && p<4.0){
	new_phi_pos = phi - (-0.000477254 + 0.0160119*eta);
	new_eta_pos =  1.39004 +  -1.11899*eta +  1.06207*eta*eta + -0.167809*eta*eta*eta;
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -= 6.12419e-02; 
      }
      else if(p>=4.0 && p<5.0){
	new_phi_pos = phi - (-0.0152078 + 0.0176643*eta);
	new_eta_pos = 1.91049 + -1.84218*eta + 1.37816*eta*eta + -0.211968*eta*eta*eta;
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -= 3.99772e-02; 
      }
      else if(p>=5.0 && p<10.0){
	new_phi_pos = phi - (0.0222429 + -0.00503467*eta);
	new_eta_pos = 0.976475 + -0.611702*eta + 0.867474*eta*eta + -0.144315*eta*eta*eta;
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -= 5.48821e-02; 
      }
      else if(p>=10.0 && p<15.0){
	new_phi_pos = phi - (-0.00554241 + 0.00817091*eta);
	new_eta_pos = 0.853172 + -0.506148*eta + 0.859377*eta*eta + -0.149469*eta*eta*eta;
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -= 6.70955e-02; 
      }
      else{
	new_phi_pos = phi - (0.00130492 + 0.00452228*eta); 
	new_eta_pos = 0.422698 + 0.166608*eta + 0.548957*eta*eta + -0.105019*eta*eta*eta; 
	// phase2
	//new_phi_pos -= ; 
	new_eta_pos -= 7.76073e-02; 
      }

    }

  }

  // set final values
  
  if(charge>0){
    phi = new_phi_pos; 
    eta = new_eta_pos;
  }
  else if(charge<0){
    phi = new_phi_neg; 
    eta = new_eta_neg;
  }
  else{
    phi = 0.5*(new_phi_pos + new_phi_neg); 
    eta = 0.5*(new_eta_pos + new_eta_neg);
  }

}

void CentauroJets::ApplyClusterOffsets( double &eta, double &phi, int charge, double p, std::string detName ){

//#define NO_INTERPOLATION

#ifdef NO_INTERPOLATION

  GetClusterOffsetFit(eta, phi, charge, p, detName); 

#else

  if(p>=15.0){
    // No interpolation for p>=15GeV 
    GetClusterOffsetFit(eta, phi, charge, p, detName); 
  }
  else{

    // interpolate results in momentum

    double bins[9] = {0.625,0.875,1.5,2.5,3.5,4.5,7.5,12.5,17.5}; 

    // Adjust low mommentum bins based on detector

    if(detName=="BECAL") {
      bins[0] = 0.6457; 
      bins[1] = 0.8587;
    }
    else if(detName=="EEMC") {
      bins[0] = 0.6505; 
      bins[1] = 0.8468;
    }
    else if(detName=="FEMC") {
      bins[0] = 0.6265; 
      bins[1] = 0.8475;
    }
    else if(detName=="HCALIN") {
      bins[0] = 0.626; 
      bins[1] = 0.8569;
    }
    else if(detName=="HCALOUT") {
      //bins[0] = ; 
      //bins[1] = ;
    }
    else if(detName=="LFHCAL") {
      bins[0] = 0.6462; 
      bins[1] = 0.8719;
    }
    
    double x_low = bins[0]; 
    double x_high = bins[8]; 

    for(int i=0; i<8; i++){
      if(p>=bins[i]) {
	x_low = bins[i];
      }
      else{

	if(i>0){
	  x_high = bins[i]; 
	}
	else{
	  x_high = bins[1]; 
	}
	  
	break; 
      }
    }

    double eta_low = eta; 
    double phi_low = phi; 
    GetClusterOffsetFit(eta_low, phi_low, charge, x_low, detName); 

    double eta_high = eta; 
    double phi_high = phi; 
    GetClusterOffsetFit(eta_high, phi_high, charge, x_high, detName); 

    eta = eta_low + (p - x_low)*( (eta_high - eta_low)/(x_high - x_low));
    phi = phi_low + (p - x_low)*( (phi_high - phi_low)/(x_high - x_low));

  }

#endif

  return; 

}

bool CentauroJets::VetoClusterWithTrack(double eta, double phi, double e, std::string detName){

  // Does this cluster have a track pointing to it? 
      
  double minDist = 9999.0; 
  double minDeta = 9999.0; 
  double minDphi = 9999.0; 
  double minP = 9999.0; 
  double minEta = 9999.0; 
  double minPhi = 9999.0;
  double minQ = 9999.0; 
  double minCEta = 9999.0; 
  double minCPhi = 9999.0;

  for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
       track_itr != _trackmap->end(); track_itr++) {

    SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

    for (SvtxTrack::ConstStateIter state_itr = temp->begin_states();
	 state_itr != temp->end_states(); state_itr++) {

      SvtxTrackState *tstate = dynamic_cast<SvtxTrackState*>(state_itr->second);			  

      if( (tstate->get_pathlength()>0.0) && (tstate->get_name()==detName) ) {

	double ceta = eta; 
	double cphi = phi; 
  
	double e_o = e; 
	if((detName=="HCALIN")||
	   (detName=="HCALOUT")||
	   (detName=="LFHCAL")){
	  e_o = tstate->get_p(); 
	}

	ApplyClusterOffsets(ceta,cphi,temp->get_charge(),e_o,detName); 

	double deta = ceta -  tstate->get_eta(); 
	double dPhi = DeltaPhi(cphi, tstate->get_phi()); 
	double dist = sqrt( pow(deta,2) + pow(dPhi,2) ); 
	if(dist<minDist){
	  minDist = dist; 
	  minDeta = deta; 
	  minDphi = dPhi; 
	  minEta = tstate->get_eta();
	  minPhi = tstate->get_phi(); 
	  minP = tstate->get_p();
	  minQ = temp->get_charge();
	  minCEta = ceta; 
	  minCPhi = cphi; 
	}

      }

    }

  }

  _tm_dist = minDist; 
  _tm_dphi = minDphi; 
  _tm_deta = minDeta;
  _tm_p = minP; 
  _tm_eta = minEta; 
  _tm_phi = minPhi; 
  _tm_ceta = minCEta; 
  _tm_cphi = minCPhi; 
  _tm_q = minQ; 
  _tm_e = e; 

  double cutDist = 0.15; 

  if(detName=="EEMC") {
    _eval_tmatch_eemc->Fill(); 
    cutDist = getMatchingCut(_tm_e, _tm_q, detName); 
  }

  if(detName=="BECAL") {
    _eval_tmatch_becal->Fill(); 
    cutDist = getMatchingCut(_tm_e, _tm_q, detName); 
  }

  if(detName=="HCALIN") {
    _eval_tmatch_ihcal->Fill(); 
    cutDist = getMatchingCut(_tm_p, _tm_q, detName); 
  }

  if(detName=="HCALOUT") {
    _eval_tmatch_ohcal->Fill(); 
    cutDist = getMatchingCut(_tm_p, _tm_q, detName); 
  }

  if(detName=="FEMC") {
    _eval_tmatch_femc->Fill(); 
    cutDist = getMatchingCut(_tm_e, _tm_q, detName); 
  }

  if(detName=="LFHCAL") {
    _eval_tmatch_lfhcal->Fill(); 
    cutDist = getMatchingCut(_tm_p, _tm_q, detName); 
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
  std::string bkwd_calos[3] = {"EEMC","",""}; 
  std::string detName[3] = {"","",""}; 

  // Collect the required cluster lists

  RawClusterContainer *clusterList[3] = {NULL, NULL, NULL};

  for(int i=0; i<3; i++){

    if(type=="CENT") 
      detName[i] = cent_calos[i]; 
    else if(type=="FWD") 
      detName[i] = fwd_calos[i]; 
    else if(type=="BKWD") 
      detName[i] = bkwd_calos[i]; 
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
      if(i==0 && rcluster->get_energy()<EM_CLUSTER_E_CUTOFF) {
	cused0[k] = true; 
	continue;
      }
      if(i>0 && rcluster->get_energy()<HAD_CLUSTER_E_CUTOFF) {
	if(i==1)
	  cused1[k] = true; 
	else
	  cused2[k] = true; 
	continue;
      }

      double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
      double phi = rcluster->get_phi(); 

      if(i==0) {tmatched0[k] = VetoClusterWithTrack(eta, phi, rcluster->get_energy(), detName[i]); cused0[k] = tmatched0[k];}
      if(i==1) {tmatched1[k] = VetoClusterWithTrack(eta, phi, rcluster->get_energy(), detName[i]); cused1[k] = tmatched1[k];}
      if(i==2) {tmatched2[k] = VetoClusterWithTrack(eta, phi, rcluster->get_energy(), detName[i]); cused2[k] = tmatched2[k];}

      if(detName[i]=="EEMC") {
	_h_clusteta_eemc->Fill(eta); 
      }

      if(detName[i]=="BECAL") {
	_h_clusteta_becal->Fill(eta); 
      }

      if(detName[i]=="HCALIN") {
	_h_clusteta_ihcal->Fill(eta); 
      }

      if(detName[i]=="HCALOUT") {
	_h_clusteta_ohcal->Fill(eta); 
      }

      if(detName[i]=="FEMC") {
	_h_clusteta_femc->Fill(eta); 
      }

      if(detName[i]=="LFHCAL") {
	_h_clusteta_lfhcal->Fill(eta); 
      }

    }

    if(detName[i]=="EEMC") {
      _h_nclusters_eemc->Fill(clusterList[i]->size());
    }

    if(detName[i]=="BECAL") {
      _h_nclusters_becal->Fill(clusterList[i]->size()); 
    }

    if(detName[i]=="HCALIN") {
      _h_nclusters_ihcal->Fill(clusterList[i]->size());
    }

    if(detName[i]=="HCALOUT") {
      _h_nclusters_ohcal->Fill(clusterList[i]->size());
    }

    if(detName[i]=="FEMC") {
      _h_nclusters_femc->Fill(clusterList[i]->size()); 
    }

    if(detName[i]=="LFHCAL") {
      _h_nclusters_lfhcal->Fill(clusterList[i]->size()); 
    }

  }

  // Mark the electron cluster as used
  // (the EMCal must always be the first in the list for this to work)
    
  if( ecDet==detName[0] ) cused0[(int)ecIdx] = true; 

  // First, seed with the EMCal and look for backing energy in the HCALs

  for (unsigned int k = 0; k < clusterList[0]->size(); k++) {

    double e_found[3] = {0.0,0.0,0.0}; 

    if(cused0[k]) continue; 

    RawCluster *rcluster0 = clusterList[0]->getCluster(k);

    double eta = getEta(rcluster0->get_r(),rcluster0->get_z()-vtx_z);
    double phi = rcluster0->get_phi(); 
 
    // Apply the track/cluster matching offsets to the eta/phi
    ApplyClusterOffsets( eta, phi, 0, 20.0, detName[0] ); 

    if(!PassClusterEtaCut(eta,detName[0])) continue; 

    double pt = rcluster0->get_energy() / cosh(eta);
    double px = pt * cos(phi);
    double py = pt * sin(phi);
    double pz = pt * sinh(eta);

    // Create the cluster
    TLorentzVector cluster(px,py,pz,rcluster0->get_energy()); 
    
    e_found[0] = rcluster0->get_energy(); 

    cused0[k] = true; 

    // mark as a photon candidate - will be set to false if we 
    // find hadronic energy further down

    bool photon_candidate = true; 

    // Look for the closest match in the next detector 

    if(detName[1]!=""){

      TLorentzVector cluster1(0.0,0.0,0.0,0.0);
      double mDist = 9999.0;
      double c_dist = 9999.0; 
      double c_deta = 9999.0; 
      double c_dPhi = 9999.0; 
      int prev_match = -9999; 

      for (unsigned int j = 0; j < clusterList[1]->size(); j++) {

	if(cused1[j]) continue; 

	RawCluster *rcluster1 = clusterList[1]->getCluster(j);

	double eta1 = getEta(rcluster1->get_r(),rcluster1->get_z()-vtx_z);
	double phi1 = rcluster1->get_phi(); 

	// Apply the track/cluster offsets
	ApplyClusterOffsets( eta1, phi1, 0, 20.0, detName[1] ); 

	if(!PassClusterEtaCut(eta1,detName[1])) continue; 

	double deta = eta -  eta1; 
	double dPhi = DeltaPhi(phi, phi1); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) );
     
	if(dist<mDist){ 

	  double pt1 = rcluster1->get_energy() / cosh(eta1);
	  double px1 = pt1 * cos(phi1);
	  double py1 = pt1 * sin(phi1);
	  double pz1 = pt1 * sinh(eta1);

	  TLorentzVector clusterAdd(px1,py1,pz1,rcluster1->get_energy()); 
	  cluster1 = clusterAdd;
	  
	  mDist = dist; 
	  c_dist = dist;
	  c_deta = deta;
	  c_dPhi = dPhi; 

	  photon_candidate = false; 

	  cused1[j] = true; 

	  if(prev_match>=0) cused1[prev_match] = false; 
	  prev_match = j; 


	}

      }

      if(cluster1.E()>0.0){

	e_found[1] = cluster1.E(); 

	_cm_deta = c_deta;
	_cm_dphi = c_dPhi;
	_cm_dist = c_dist;
	_cm_eta = cluster1.Eta();
	_cm_phi = cluster1.Phi();
	_cm_E = cluster1.E();

	// tuned matching cut
	double dcut = 9999.0; 

	if(type=="CENT"){
	  dcut = 2.5*0.171; 
	  _eval_cmatch_becal_ihcal->Fill();
	}
	else if(type=="FWD"){
	  dcut = 0.7; 
	  _eval_cmatch_femc_lfhcal->Fill();
	}

	if(c_dist<dcut){
	  // Add what we found to the existing cluster
	  cluster += cluster1;
	}
	else{
	  if(prev_match>=0) cused1[prev_match] = false; 
	}

      }

    }

    // And the last detector (if it exists)

    if(detName[2]!="") {

      TLorentzVector cluster2(0.0,0.0,0.0,0.0); 
      double mDist = 9999.0;
      double c_dist = 9999.0; 
      double c_deta = 9999.0; 
      double c_dPhi = 9999.0; 
      int prev_match = -9999; 

      for (unsigned int j = 0; j < clusterList[2]->size(); j++) {

	if(cused2[j]) continue; 

	RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	double eta2 = getEta(rcluster2->get_r(),rcluster2->get_z()-vtx_z);
	double phi2 = rcluster2->get_phi(); 

	// Apply the track/cluster offsets
	ApplyClusterOffsets( eta2, phi2, 0, 20.0, detName[2] ); 

	double deta = eta -  eta2; 
	double dPhi = DeltaPhi(phi, phi2); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) ); 
	
	if(dist<mDist){ 

	  double pt2 = rcluster2->get_energy() / cosh(eta2);
	  double px2 = pt2 * cos(phi2);
	  double py2 = pt2 * sin(phi2);
	  double pz2 = pt2 * sinh(eta2);

	  TLorentzVector clusterAdd(px2,py2,pz2,rcluster2->get_energy()); 
	  cluster2 = clusterAdd; 

	  mDist = dist; 
	  c_dist = dist;
	  c_deta = deta;
	  c_dPhi = dPhi; 

	  photon_candidate = false; 

	  cused2[j] = true; 
	  
	  if(prev_match>=0) cused2[prev_match] = false; 
	  prev_match = j; 

	}

      }

      if(cluster2.E()>0.0){

	e_found[2] = cluster2.E(); 
	
	// tuned matching cut 
	double dcut = 2.5*0.224; 

	_cm_deta = c_deta;
	_cm_dphi = c_dPhi;
	_cm_dist = c_dist;
	_cm_eta = cluster2.Eta();
	_cm_phi = cluster2.Phi();
	_cm_E = cluster2.E();

	_eval_cmatch_becal_ohcal->Fill(); 

	if(c_dist<dcut){
	  // Add what we found to the existing cluster
	  cluster += cluster2;
	}
	else{
	  if(prev_match>=0) cused2[prev_match] = false; 
	}

      }

    }

    // Adjust the energy scale

    if(!photon_candidate){
      double scale = GetCaloTrackHadronicEnergyScale(e_found, type);
      if(scale<=0.0) continue; 
      cluster = cluster*scale;
    }
    
    // Evaluation

    TVector3 ctrack = cluster.Vect();

    cat_e_tot = ctrack.Mag();
    cat_eta_meas = ctrack.Eta(); 
    cat_phi_meas = ctrack.Phi(); 

    GetCaloTrackTruthInfo( ctrack, type ); 

    if(type=="CENT") {
      cat_e_bemc = e_found[0]; 
      cat_e_ihcal = e_found[1]; 
      cat_e_ohcal = e_found[2]; 
      _eval_calo_tracks_cent->Fill();
    }
    if(type=="FWD") {
      cat_e_femc = e_found[0]; 
      cat_e_lfhcal = e_found[1]; 
      _eval_calo_tracks_fwd->Fill();
    }
    if(type=="BKWD") {
      cat_e_eemc = e_found[0]; 
      _eval_calo_tracks_bkwd->Fill(); 
    }
      
    // Finally - take the combined cluster and add it to the constituents
    // Transform to Breit frame
    TLorentzVector breit_cluster = (breit*cluster); 
    breit_cluster.Transform(breitRot); 

    if( PARAM_HAD_NEUTRALS ){ 
    
      if(photon_candidate){
	
	// require a loose match to the nearest primary photon
	if((cat_match<0.1) && (cat_pid==22)){
	  fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
	  pseudojet.set_user_index(EncodeUserIndex(0,photon_candidate)); 
	  pseudojets.push_back(pseudojet);
	}
	  
      }

    }
    else{

      fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
      pseudojet.set_user_index(EncodeUserIndex(0,photon_candidate)); 
      pseudojets.push_back(pseudojet);

    }

  }

  // Next, seed with the HCAL and follow up with the second HCAL segment
  // necessary to capture hadronic showers w/o cluster in EMCal

  if(detName[1]!=""){

    for (unsigned int k = 0; k < clusterList[1]->size(); k++) {

      double e_found[3] = {0.0,0.0,0.0}; 
    
      if(cused1[k]) continue;

      RawCluster *rcluster = clusterList[1]->getCluster(k);

      double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
      double phi = rcluster->get_phi(); 

      // Apply the track/cluster offsets
      ApplyClusterOffsets( eta, phi, 0, 20.0, detName[1] ); 

      if(!PassClusterEtaCut(eta,detName[1])) continue; 

      double pt = rcluster->get_energy() / cosh(eta);
      double px = pt * cos(phi);
      double py = pt * sin(phi);
      double pz = pt * sinh(eta);

      // Create the cluster
      TLorentzVector cluster(px,py,pz,rcluster->get_energy()); 

      e_found[0] = cluster.E(); 

      cused1[k] = true; 

      // And the last detector (if it exists)

      if(detName[2]!="") {

	TLorentzVector cluster2(0.0,0.0,0.0,0.0); 
	double mDist = 9999.0;
	double c_dist = 9999.0; 
	double c_deta = 9999.0; 
	double c_dPhi = 9999.0; 
	int prev_match = -9999; 

 	for (unsigned int j = 0; j < clusterList[2]->size(); j++) {

	  if(cused2[j]) continue; 

	  RawCluster *rcluster2 = clusterList[2]->getCluster(j);

	  double eta2 = getEta(rcluster2->get_r(),rcluster2->get_z()-vtx_z);
	  double phi2 = rcluster2->get_phi(); 

	  // Apply the track/cluster offsets
	  ApplyClusterOffsets( eta2, phi2, 0, 20.0, detName[2] ); 

	  if(!PassClusterEtaCut(eta2,detName[2])) continue; 

	  double deta = eta -  eta2; 
	  double dPhi = DeltaPhi(phi, phi2); 

	  double dist = sqrt( pow(deta,2) + pow(dPhi,2) );

	  if(dist<mDist){ 

	    double pt2 = rcluster2->get_energy() / cosh(eta2);
	    double px2 = pt2 * cos(phi2);
	    double py2 = pt2 * sin(phi2);
	    double pz2 = pt2 * sinh(eta2);

	    TLorentzVector clusterAdd(px2,py2,pz2,rcluster2->get_energy()); 
	    cluster2 = clusterAdd; 

	    mDist = dist; 
	    c_dist = dist;
	    c_deta = deta;
	    c_dPhi = dPhi; 

	    cused2[j] = true; 

	    if(prev_match>=0) cused2[prev_match] = false; 
	    prev_match = j; 

	  }

	}

	if(cluster2.E()>0.0){
	  
	  e_found[1] = cluster2.E(); 

	  if(type=="CENT"){

	    _cm_deta = c_deta;
	    _cm_dphi = c_dPhi;
	    _cm_dist = c_dist;
	    _cm_eta = cluster2.Eta();
	    _cm_phi = cluster2.Phi();
	    _cm_E = cluster2.E();

	    _eval_cmatch_ihcal_ohcal->Fill();
 
	  }

	  double dcut = 0.5; 

	  if(c_dist<dcut){
	    // Add what we found to the existing cluster
	    cluster += cluster2;
	  }
	  else{
	    if(prev_match>=0) cused2[prev_match] = false; 
	  }

	}

      }

      // Adjust the hadronic energy scale
      double scale = GetCaloTrackHadronicEnergyScale(e_found, type);
      if(scale<=0.0) continue; 
      cluster = cluster*scale;

      // Finally - take the combined cluster and add it to the constituents
      // Transform to Breit frame
      TLorentzVector breit_cluster = (breit*cluster); 
      breit_cluster.Transform(breitRot); 

      if( !PARAM_HAD_NEUTRALS ){ 

        fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
        pseudojet.set_user_index(0); 
        pseudojets.push_back(pseudojet);
      
      }

      TVector3 ctrack = cluster.Vect();

      cat_e_tot = ctrack.Mag();
      cat_eta_meas = ctrack.Eta(); 
      cat_phi_meas = ctrack.Phi(); 

      GetCaloTrackTruthInfo( ctrack, type ); 

      if(type=="CENT") {
	cat_e_bemc = 0.0; 
	cat_e_ihcal = e_found[0]; 
	cat_e_ohcal = e_found[1]; 
	_eval_calo_tracks_cent->Fill();
      }

    }

  }

  // Finally - the last HCAL by itself (if it exists)

  if(detName[2]!="") {

    for (unsigned int k = 0; k < clusterList[2]->size(); k++) {

      if(cused2[k]) continue; 

      RawCluster *rcluster1 = clusterList[2]->getCluster(k);

      double eta1 = getEta(rcluster1->get_r(),rcluster1->get_z()-vtx_z);
      double phi1 = rcluster1->get_phi(); 

      // Apply the track/cluster offsets
      ApplyClusterOffsets( eta1, phi1, 0, 20.0, detName[2] ); 

      if(!PassClusterEtaCut(eta1,detName[2])) continue; 

      double pt1 = rcluster1->get_energy() / cosh(eta1);
      double px1 = pt1 * cos(phi1);
      double py1 = pt1 * sin(phi1);
      double pz1 = pt1 * sinh(eta1);

      // Create the cluster
      TLorentzVector cluster(px1,py1,pz1,rcluster1->get_energy()); 

      cused2[k] = true; 

      // Adjust the hadronic energy scale
      double e_found[3] = {0.0,0.0,0.0}; 
      e_found[2] = rcluster1->get_energy(); 
      double scale = GetCaloTrackHadronicEnergyScale(e_found, type);
      if(scale<=0.0) continue; 
      cluster = cluster*scale;

      // Transform to Breit frame
      TLorentzVector breit_cluster = (breit*cluster); 
      breit_cluster.Transform(breitRot); 

      if( !PARAM_HAD_NEUTRALS ){ 

        fastjet::PseudoJet pseudojet (breit_cluster.Px(),breit_cluster.Py(),breit_cluster.Pz(),breit_cluster.E()); 
        pseudojet.set_user_index(0); 
        pseudojets.push_back(pseudojet);

      }

      TVector3 ctrack = cluster.Vect();

      cat_e_tot = ctrack.Mag();
      cat_eta_meas = ctrack.Eta(); 
      cat_phi_meas = ctrack.Phi(); 

      GetCaloTrackTruthInfo( ctrack, type ); 

      if(type=="CENT") {
	cat_e_bemc = 0.0; 
	cat_e_ihcal = 0.0; 
	cat_e_ohcal = cluster.E(); 
	_eval_calo_tracks_cent->Fill();
      }

    }

  }

  return; 

}

void CentauroJets::GetCaloTrackTruthInfo( TVector3 ctrack, std::string type ){

  // Diagnostics - connect the calo tracks to neutral primaries

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  double minDist = 9999.0; 
  int pid = -9999; 
  double prim_p = 9999.0; 
  double prim_Eta = 9999.0; 
  double prim_Phi = 9999.0; 

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

    TLorentzVector partMom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e()); 

    double deta = ctrack.Eta() - partMom.Eta(); 
    double dphi = DeltaPhi(ctrack.Phi(), partMom.Phi()); 

    double dist = sqrt(pow(deta,2) + pow(dphi,2)); 
    if(dist<minDist){
      minDist = dist; 
      pid = g4particle->get_pid();
      prim_p = partMom.Vect().Mag();
      prim_Eta = partMom.Vect().Eta();
      prim_Phi = partMom.Vect().Phi(); 
    }

  }

  cat_pid = pid;
  cat_p_true = prim_p; 
  cat_e_tot = ctrack.Mag();
  cat_match = minDist; 
  cat_eta_meas = ctrack.Eta(); 
  cat_eta_true = prim_Eta; 
  cat_phi_meas = ctrack.Phi(); 
  cat_phi_true = prim_Phi; 

  return; 

}

void CentauroJets::BuildParametrizedPhotonCaloTracks(PHCompositeNode *topNode, std::string type, 
				   std::vector<fastjet::PseudoJet> &pseudojets, 
				   TLorentzRotation &breit, TRotation &breitRot){

  // get the list of jets from the primary particles

  if (!_truth_container) {
    LogError("_truth_container not found!");
    return;
  }

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

    // All we do here is photons
    if(g4particle->get_pid()!=22) continue; 
    int charge = 0; 
 
    // lab frame cuts
    TLorentzVector partMom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e()); 

    double res_a = 0.0; 
    double res_b = 0.0; 

    if(type=="CENT"){
      if (fabs(partMom.Eta())>1.0) continue; 
      res_a = 0.015; 
      res_b = 0.019; 
    }
    else if(type=="FWD"){
      if((partMom.Eta()<1.3) || (partMom.Eta()>3.5)) continue; 
      res_a = 0.072; 
      res_b = 0.00; 
    }
    else if(type=="BKWD"){
      if((partMom.Eta()>-1.8) || (partMom.Eta()<-3.5)) continue; 
      res_a = 0.016; 
      res_b = 0.002; 
    }
    else{
      return; 
    }
					      
    // Parametrize the calorimeter response: 
    
    double c_energy = g4particle->get_e();
    double res = sqrt( pow(res_a/sqrt(c_energy),2) + pow(res_b,2) ); 
    c_energy = rand->Gaus(c_energy, res*c_energy); 
    // minimal cluster energy
    if(c_energy<0.200) continue;

    double ratio = c_energy/partMom.P(); 

    TLorentzVector partMomP(g4particle->get_px()*ratio,g4particle->get_py()*ratio,g4particle->get_pz()*ratio,c_energy); 

    // add this track to the list of tracks for jets

    TLorentzVector partMom_breit = (breit*partMomP); 
    partMom_breit.Transform(breitRot); 

    fastjet::PseudoJet pseudojet (partMom_breit.Px(),
				  partMom_breit.Py(),
				  partMom_breit.Pz(),
				  partMom_breit.E());

    // build the user index and add the particle
    bool em_part = true; 
    pseudojet.set_user_index(EncodeUserIndex(charge,em_part));
    pseudojets.push_back(pseudojet);

  }

  return; 

}

void CentauroJets::BuildParametrizedHadCaloTracks(PHCompositeNode *topNode, std::string type, 
				   std::vector<fastjet::PseudoJet> &pseudojets, 
				   TLorentzRotation &breit, TRotation &breitRot){

  // No HCAL in BKWD direction
  if(type=="BKWD") return; 

  // get the list of jets from the primary particles

  if (!_truth_container) {
    LogError("_truth_container not found!");
    return;
  }

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
    if((abs(g4particle->get_pid())==11) || (g4particle->get_pid()==22)) continue; 

    // eliminate charged 
    G4ParticleDefinition* particle = particleTable->FindParticle(g4particle->get_name());
    int charge = -9999.0; 
    if(particle) 
      charge = particle->GetPDGCharge();
    else
      continue; 
    if(charge!=0) continue;
 
    // lab frame cuts
    TLorentzVector partMom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e()); 

    double res_a = 0.0; 
    double res_b = 0.0; 

    if(type=="CENT"){
      if (fabs(partMom.Eta())>1.0) continue; 
      res_a = 0.74; 
      res_b = 0.16; 
    }
    else if(type=="FWD"){
      if((partMom.Eta()<1.3) || (partMom.Eta()>3.5)) continue; 
      res_a = 0.31; 
      res_b = 0.028; 
    }
    else{
      return; 
    }
					      
    // Parametrize the calorimeter response: 
    
    double c_energy = g4particle->get_e();
    if(g4particle->get_pid()==2112){
      c_energy = partMom.E() - partMom.M(); // kinetic energy for neutrons
    }
    // minimal cluster energy
    if(c_energy<0.500) continue; 

    double res = sqrt( pow(res_a/sqrt(c_energy),2) + pow(res_b,2) ); 
    c_energy = rand->Gaus(c_energy, res*c_energy); 
    if(c_energy<0.0) continue;

    double ratio = c_energy/partMom.P(); 

    TLorentzVector partMomP(g4particle->get_px()*ratio,g4particle->get_py()*ratio,g4particle->get_pz()*ratio,c_energy); 

    // add this track to the list of tracks for jets

    TLorentzVector partMom_breit = (breit*partMomP); 
    partMom_breit.Transform(breitRot); 

    fastjet::PseudoJet pseudojet (partMom_breit.Px(),
				  partMom_breit.Py(),
				  partMom_breit.Pz(),
				  partMom_breit.E());

    // build the user index and add the particle
    bool em_part = false; 
    pseudojet.set_user_index(EncodeUserIndex(charge,em_part));
    pseudojets.push_back(pseudojet);

  }

  return; 

}

SvtxTrack *CentauroJets::AttachClusterToTrack(double eta, double phi, double e, std::string detName){

  // Does this cluster have a track pointing to it? 
      
  double minDist = 9999.0; 
  double minEO = 9999.0;
  int minQ = 9999; 
  SvtxTrack *closest = NULL; 

  for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
       track_itr != _trackmap->end(); track_itr++) {

    SvtxTrack* temp = dynamic_cast<SvtxTrack*>(track_itr->second);

    for (SvtxTrack::ConstStateIter state_itr = temp->begin_states();
	 state_itr != temp->end_states(); state_itr++) {

      SvtxTrackState *tstate = dynamic_cast<SvtxTrackState*>(state_itr->second);			  

      if( (tstate->get_pathlength()>0.0) && (tstate->get_name()==detName) ) {

	double ceta = eta; 
	double cphi = phi; 

	double e_o = e; 
	if((detName=="HCALIN")||
	   (detName=="HCALOUT")||
	   (detName=="LFHCAL")){
	  e_o = tstate->get_p(); 
	}

	ApplyClusterOffsets(ceta,cphi,temp->get_charge(),e_o, detName); 

	double deta = ceta -  tstate->get_eta(); 
	double dPhi = DeltaPhi(cphi, tstate->get_phi()); 

	double dist = sqrt( pow(deta,2) + pow(dPhi,2) ); 
	if(dist<minDist){
	  minDist = dist;
	  closest = temp;
	  minEO = e_o; 
	  minQ = temp->get_charge();
	}

      }

    }

  }

  double cutDist = getMatchingCut(minEO, minQ, detName); 

  // Return the track the cluster is closest to 
  if(minDist<cutDist)
    return closest;
  else
    return NULL; 

}

void CentauroJets::BuildChargedCaloTracks(PHCompositeNode *topNode, std::string type){

  std::string fwd_calos[3] = {"FEMC","LFHCAL",""}; 
  std::string cent_calos[3] = {"BECAL","HCALIN","HCALOUT"}; 
  std::string bkwd_calos[3] = {"EEMC","",""}; 
  std::string detName[3] = {"","",""}; 

  // Collect the required cluster lists

  RawClusterContainer *clusterList[3] = {NULL, NULL, NULL};

  for(int i=0; i<3; i++){

    if(type=="CENT") 
      detName[i] = cent_calos[i]; 
    else if(type=="FWD") 
      detName[i] = fwd_calos[i]; 
    else if(type=="BKWD") 
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
      if(i==0 && rcluster->get_energy()<EM_CLUSTER_E_CUTOFF) {
	cused0[k] = true; 
	continue;
      }
      if(i>0 && rcluster->get_energy()<HAD_CLUSTER_E_CUTOFF) {
	if(i==1)
	  cused1[k] = true; 
	else
	  cused2[k] = true; 
	continue;
      }

      double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
      double phi = rcluster->get_phi();

      if(i==0) {tmatched0[k] = AttachClusterToTrack(eta, phi, rcluster->get_energy(), detName[i]);}
      if(i==1) {tmatched1[k] = AttachClusterToTrack(eta, phi, rcluster->get_energy(), detName[i]);}
      if(i==2) {tmatched2[k] = AttachClusterToTrack(eta, phi, rcluster->get_energy(), detName[i]);}

    }

  }

  // Attach the clusters to the tracks and fill the tree
  // start seeded with the first calorimeter (EMCAL)

  for (unsigned int k = 0; k < clusterList[0]->size(); k++) {

    if(!tmatched0[k]) continue; 
    if(cused0[k]) continue; 

    ct_pid = -9999; 
    ct_p_meas = -9999.0;  
    ct_p_true = -9999.0; 
    ct_eta_meas = -9999.0; 
    ct_eta_true = -9999.0; 
    ct_phi_meas = -9999.0; 
    ct_phi_true = -9999.0; 
    ct_e_eemc = -9999.0; 
    ct_e_bemc = -9999.0; 
    ct_e_ihcal = -9999.0; 
    ct_e_ohcal = -9999.0;  
    ct_e_femc = -9999.0; 
    ct_e_lfhcal = -9999.0; 
    ct_e_tot = -9999.0; 

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

	ct_pid = g4particle->get_pid(); 

	// Get the final particle in the lab frame
	CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),
				g4particle->get_pz(),g4particle->get_e());
	// Not needed per Jin
	//efp = EventToLab * efp;  
	TLorentzVector match_lf(efp.px(), efp.py(), efp.pz(), efp.e());

	ct_p_meas = tmatched0[k]->get_p(); 
	ct_p_true = match_lf.Vect().Mag(); 
		
	ct_eta_meas = tmatched0[k]->get_eta(); 
        ct_eta_true = match_lf.Vect().Eta(); 

	ct_phi_meas = tmatched0[k]->get_phi(); 
        ct_phi_true = match_lf.Vect().Phi(); 

	cused0[k] = true; 

      }
	
    }
    
    if(!found) continue; 

    double caloTot = rcluster->get_energy(); 
 
    if(type=="CENT"){
      ct_e_bemc = rcluster->get_energy(); 
    }
    else if(type=="FWD"){
      ct_e_femc = rcluster->get_energy(); 
    }
    else if(type=="BKWD"){
      ct_e_eemc = rcluster->get_energy(); 
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
	  stage2_energy += rcluster1->get_energy(); 

	  cused1[j] = true; 

	  //break; 

	}

      }

    }

    if(type=="CENT"){
      ct_e_ihcal = stage2_energy; 
    }
    else if(type=="FWD"){
      ct_e_lfhcal = stage2_energy; 
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
      ct_e_ohcal = stage3_energy; 
    }

    // record the combined energy
    ct_e_tot = caloTot; 
  
    if(type=="CENT"){

      double e_found[3] = {0.0,0.0,0.0};
      e_found[0] = ct_e_bemc; 
      e_found[1] = ct_e_ihcal; 
      e_found[2] = ct_e_ohcal;

      double scale =  GetCaloTrackHadronicEnergyScale(e_found, type);
      if(scale<=0.0) continue; 
      
      ct_e_tot *= scale; 
      
      _eval_charged_tracks_cent->Fill(); 
    }
    else if(type=="FWD"){

      double e_found[3] = {0.0,0.0,0.0};
      e_found[0] = ct_e_femc; 
      e_found[1] = ct_e_lfhcal;

      double scale =  GetCaloTrackHadronicEnergyScale(e_found, type);
      if(scale<=0.0) continue; 
      
      ct_e_tot *= scale; 

      _eval_charged_tracks_fwd->Fill();  
    }
    else if(type=="BKWD"){
      _eval_charged_tracks_bkwd->Fill();  
    }

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

	    ct_pid = g4particle->get_pid(); 

	    // Get the final particle in the lab frame
	    CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),
					g4particle->get_pz(),g4particle->get_e());
	    // Not needed per Jin 
	    //efp = EventToLab * efp;  
	    TLorentzVector match_lf(efp.px(), efp.py(), efp.pz(), efp.e());

	    ct_p_meas = tmatched1[k]->get_p(); 
	    ct_p_true = match_lf.Vect().Mag(); 
	
	    ct_eta_meas = tmatched1[k]->get_eta(); 
	    ct_eta_true = match_lf.Vect().Eta(); 

	    found = true;

	    cused1[k] = true; 

	  }

	  stage1_energy += rcluster->get_energy(); 

	  //break; 

	}
	
      }
     
      if(!found) continue;

      double caloTot = stage1_energy; 
 
      if(type=="CENT"){
	ct_e_bemc = 0.0; 
	ct_e_ihcal = stage1_energy; 
      }
      else if(type=="FWD"){
	ct_e_femc = 0.0; 
	ct_e_lfhcal = stage1_energy; 
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
	ct_e_ohcal = stage3_energy; 
      }

      // record the combined energy
      ct_e_tot = caloTot; 

      if(type=="CENT"){

	double e_found[3] = {0.0,0.0,0.0};
	e_found[0] = ct_e_bemc; 
	e_found[1] = ct_e_ihcal; 
	e_found[2] = ct_e_ohcal;

	double scale =  GetCaloTrackHadronicEnergyScale(e_found, type);
	if(scale<=0.0) continue; 
      
	ct_e_tot *= scale; 
 
	_eval_charged_tracks_cent->Fill(); 
      }
      else if(type=="FWD"){

	double e_found[3] = {0.0,0.0,0.0};
	e_found[0] = ct_e_femc; 
	e_found[1] = ct_e_lfhcal;

	double scale =  GetCaloTrackHadronicEnergyScale(e_found, type);
	if(scale<=0.0) continue; 
      
	ct_e_tot *= scale; 

	_eval_charged_tracks_fwd->Fill();  
      }

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

	  ct_pid = g4particle->get_pid(); 

	  // Get the final particle in the lab frame
	  CLHEP::HepLorentzVector efp(g4particle->get_px(),g4particle->get_py(),
				      g4particle->get_pz(),g4particle->get_e());
	  // Not needed per Jin
	  //efp = EventToLab * efp;  
	  TLorentzVector match_lf(efp.px(), efp.py(), efp.pz(), efp.e());

	  ct_p_meas = tmatched2[k]->get_p(); 
	  ct_p_true = match_lf.Vect().Mag(); 
	
	  ct_eta_meas = tmatched2[k]->get_eta(); 
	  ct_eta_true = match_lf.Vect().Eta(); 

	  cused2[k] = true; 

	  break; 

	}
	
      }
    
      if(!found) continue; 

      double caloTot = rcluster->get_energy(); 
 
      if(type=="CENT"){
	ct_e_bemc = 0.0; 
	ct_e_ihcal = 0.0; 
	ct_e_ohcal = rcluster->get_energy(); 
      }

      // record the combined energy
      ct_e_tot = caloTot; 

      if(type=="CENT"){

	double e_found[3] = {0.0,0.0,0.0};
	e_found[0] = ct_e_bemc; 
	e_found[1] = ct_e_ihcal; 
	e_found[2] = ct_e_ohcal;

	double scale =  GetCaloTrackHadronicEnergyScale(e_found, type);
	if(scale<=0.0) continue; 
      
	ct_e_tot *= scale; 
 
	_eval_charged_tracks_cent->Fill(); 
      }

    }

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

    // eliminate noise clusters
    if(rcluster->get_energy()<EM_CLUSTER_E_CUTOFF) continue; 

    double eta = getEta(rcluster->get_r(),rcluster->get_z()-vtx_z);
    double phi = rcluster->get_phi(); 

    ApplyClusterOffsets(eta, phi, 0.0, 20.0, detName); 

    if(TrackVeto){
      if(VetoClusterWithTrack(eta, phi, rcluster->get_energy(), detName)) continue; 
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
				       double eta, double phi, int charge, double ptot, 
				       double &clustE, int &clIdx, double &dR){

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

    double cluster_phi = cluster->get_phi(); 
    double cluster_eta = getEta(cluster->get_r(),cluster->get_z()); 

    ApplyClusterOffsets( cluster_eta, cluster_phi, charge, ptot, detName );
    if(!PassClusterEtaCut(cluster_eta, detName )) continue; 

    double dphi = DeltaPhi(cluster_phi,phi);
    double deta = cluster_eta - eta;
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
				  TLorentzRotation &breit, TRotation &breitRot,TLorentzRotation &breitInv, TRotation &breitRotInv, 
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
    // No transform needed per Jin!
    //efp = EventToLab * efp;  

    TLorentzVector partMom(efp.px(), efp.py(), efp.pz(), efp.e()); 

    G4ParticleDefinition* particle = particleTable->FindParticle(g4particle->get_name());
    int charge = -9999.0; 
    if(particle) 
      charge = particle->GetPDGCharge();
    else
      continue; 

    // lab frame cuts
    if(((charge!=0)&&(partMom.Pt()<0.200))||
        (fabs(partMom.Eta())>3.5)) continue;
    //if(((charge!=0)&&(partMom.Pt()<0.100))||
    //    (fabs(partMom.Eta())>3.5)) continue;
   
    TLorentzVector partMom_breit = (breit*partMom); 
    partMom_breit.Transform(breitRot); 

    // add this track to the list of tracks for jets

    fastjet::PseudoJet pseudojet (partMom_breit.Px(),
				  partMom_breit.Py(),
				  partMom_breit.Pz(),
				  partMom_breit.E());

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

	tfpjet_cf.push_back(JetChargedFraction(&tconstit,pjet.Vect()));

	tfpjet_neut_p.push_back(JetNeutralMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 
	tfpjet_chgd_p.push_back(JetChargedMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 
	tfpjet_em_p.push_back(JetEMMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 

	tfpjet_neut_pm.push_back(JetNeutralMomentum(&tconstit).Mag()); 
	tfpjet_chgd_pm.push_back(JetChargedMomentum(&tconstit).Mag()); 
	tfpjet_em_pm.push_back(JetEMMomentum(&tconstit).Mag()); 

	// Transform back to lab frame
	pjet.Transform(breitRotInv); 
	pjet = (breitInv * pjet); 

        tfpjet_lab_eta.push_back(pjet.Eta());  
	tfpjet_lab_phi.push_back(pjet.Phi()); 
	tfpjet_lab_p.push_back(pjet.Vect().Mag()); 

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

	pjet_cf.push_back(JetChargedFraction(&tconstit,pjet.Vect()));

	pjet_neut_p.push_back(JetNeutralMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 
	pjet_chgd_p.push_back(JetChargedMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 
	pjet_em_p.push_back(JetEMMomentum(&tconstit).Dot(pjet.Vect())/pjet.Vect().Mag()); 

	pjet_neut_pm.push_back(JetNeutralMomentum(&tconstit).Mag()); 
	pjet_chgd_pm.push_back(JetChargedMomentum(&tconstit).Mag()); 
	pjet_em_pm.push_back(JetEMMomentum(&tconstit).Mag()); 

	// Transform back to lab frame
	pjet.Transform(breitRotInv); 
	pjet = (breitInv * pjet); 

        pjet_lab_eta.push_back(pjet.Eta());  
	pjet_lab_phi.push_back(pjet.Phi()); 
	pjet_lab_p.push_back(pjet.Vect().Mag()); 

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
