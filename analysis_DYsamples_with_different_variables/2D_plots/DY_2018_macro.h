//////////////////////////////////////////////////////////
// this class has been automatically generated on
// Wed Sep 22 10:00:04 2021 by ROOT version 6.24/04
// from TTree passedEvents/passedEvents
// found on file: DY_2018.root
//////////////////////////////////////////////////////////


// this program is to generate the ntuples with histograms holding the GEN Z mass and their pt, eta and phi 
// further it will create the histogram to hold the value of GENZ in different pt, eta and phi categories
// further there is a boolean debug_pT which do the purpose of checking the number of entries falling in each pT bins before and after smearing to chek that bins are getting filled right

#ifndef DY_2018_macro_h
#define DY_2018_macro_h
#include "TAxis.h"
#include "TH1.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
const int pT_bins = 13;
const int eta_bins = 14;
class DY_2018_macro {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        GENmass2l;
   Double_t        genLep_pt1;
   Double_t        genLep_pt2;
   Double_t        genLep_eta1;
   Double_t        genLep_eta2;
   Double_t        genLep_phi1;
   Double_t        genLep_phi2;
   Int_t           NInt;
   Int_t           NVtx;
   Double_t        PvX;
   Double_t        PvY;
   Double_t        PvZ;
   Double_t        BsX;
   Double_t        BsY;
   Double_t        BsZ;
   Double_t        massZ;
   Double_t        massZErr;
   Double_t        massZ_FSR;
   Double_t        massZErr_FSR;
   Double_t        massZ_single_BS;
   Double_t        massZ_vtx;
   Double_t        massZ_vtx_FSR;
   Double_t        massZErr_vtx;
   Double_t        massZErr_vtx_FSR;
   Double_t        massZ_vtxChi2;
   Double_t        massZ_vtx_BS;
   Double_t        massZ_vtx_BS_FSR;
   Double_t        massZErr_vtx_BS;
   Double_t        massZErr_vtx_BS_FSR;
   Double_t        massZ_vtxChi2_BS;
   Double_t        PvX1;
   Double_t        PvX2;
   Double_t        PvY1;
   Double_t        PvY2;
   Double_t        PvZ1;
   Double_t        PvZ2;
   Double_t        BsX1;
   Double_t        BsX2;
   Double_t        BsY1;
   Double_t        BsY2;
   Double_t        BsZ1;
   Double_t        BsZ2;
   Double_t        muonPV_x1;
   Double_t        muonPV_y1;
   Double_t        muonPV_z1;
   Double_t        muonPV_x2;
   Double_t        muonPV_y2;
   Double_t        muonPV_z2;
   Double_t        pT1_FromMuonBestTrack;
   Double_t        pT2_FromMuonBestTrack;
   Double_t        eta1_FromMuonBestTrack;
   Double_t        eta2_FromMuonBestTrack;
   Double_t        phi1_FromMuonBestTrack;
   Double_t        phi2_FromMuonBestTrack;
   Double_t        pT1_genFromReco;
   Double_t        pT2_genFromReco;
   Double_t        pT1;
   Double_t        pT2;
   Double_t        eta1;
   Double_t        eta2;
   Double_t        phi1;
   Double_t        phi2;
   Double_t        m1;
   Double_t        m2;
   Double_t        pT_FSR1;
   Double_t        pT_FSR2;
   Double_t        eta_FSR1;
   Double_t        eta_FSR2;
   Double_t        phi_FSR1;
   Double_t        phi_FSR2;
   Double_t        m_FSR1;
   Double_t        m_FSR2;
   Double_t        Tracker1;
   Double_t        Tracker2;
   Double_t        single_BS_pT1;
   Double_t        single_BS_pT2;
   Double_t        single_BS_eta1;
   Double_t        single_BS_eta2;
   Double_t        single_BS_phi1;
   Double_t        single_BS_phi2;
   Double_t        single_BS_m1;
   Double_t        single_BS_m2;
   Double_t        vtx_pT1;
   Double_t        vtx_pT2;
   Double_t        vtx_eta1;
   Double_t        vtx_eta2;
   Double_t        vtx_phi1;
   Double_t        vtx_phi2;
   Double_t        vtx_m1;
   Double_t        vtx_m2;
   Double_t        vtx_BS_pT1;
   Double_t        vtx_BS_pT2;
   Double_t        vtx_BS_eta1;
   Double_t        vtx_BS_eta2;
   Double_t        vtx_BS_phi1;
   Double_t        vtx_BS_phi2;
   Double_t        vtx_BS_m1;
   Double_t        vtx_BS_m2;
   Double_t        vtx_pT_FSR1;
   Double_t        vtx_pT_FSR2;
   Double_t        vtx_eta_FSR1;
   Double_t        vtx_eta_FSR2;
   Double_t        vtx_phi_FSR1;
   Double_t        vtx_phi_FSR2;
   Double_t        vtx_m_FSR1;
   Double_t        vtx_m_FSR2;
   Double_t        vtx_BS_pT_FSR1;
   Double_t        vtx_BS_pT_FSR2;
   Double_t        vtx_BS_eta_FSR1;
   Double_t        vtx_BS_eta_FSR2;
   Double_t        vtx_BS_phi_FSR1;
   Double_t        vtx_BS_phi_FSR2;
   Double_t        vtx_BS_m_FSR1;
   Double_t        vtx_BS_m_FSR2;
   Double_t        d0BS1;
   Double_t        d0BS2;
   Double_t        d0BS_vtx_BS1;
   Double_t        d0BS_vtx_BS2;
   Int_t           Id1;
   Int_t           Id2;
   Double_t        pterr1;
   Double_t        pterr2;
   Double_t        pterr1_VX;
   Double_t        pterr2_VX;
   Double_t        pterr1_VX_BS;
   Double_t        pterr2_VX_BS;
   Double_t        weight;
   Int_t           nFSRPhotons;

   // List of branches
   TBranch        *b_GENmass2l;   //!
   TBranch        *b_genLep_pt1;   //!
   TBranch        *b_genLep_pt2;   //!
   TBranch        *b_genLep_eta1;   //!
   TBranch        *b_genLep_eta2;   //!
   TBranch        *b_genLep_phi1;   //!
   TBranch        *b_genLep_phi2;   //!
   TBranch        *b_NInt;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_PvX;   //!
   TBranch        *b_PvY;   //!
   TBranch        *b_PvZ;   //!
   TBranch        *b_BsX;   //!
   TBranch        *b_BsY;   //!
   TBranch        *b_BsZ;   //!
   TBranch        *b_massZ;   //!
   TBranch        *b_massZErr;   //!
   TBranch        *b_massZ_FSR;   //!
   TBranch        *b_massZErr_FSR;   //!
   TBranch        *b_massZ_single_BS;   //!
   TBranch        *b_massZ_vtx;   //!
   TBranch        *b_massZ_vtx_FSR;   //!
   TBranch        *b_massZErr_vtx;   //!
   TBranch        *b_massZErr_vtx_FSR;   //!
   TBranch        *b_massZ_vtxChi2;   //!
   TBranch        *b_massZ_vtx_BS;   //!
   TBranch        *b_massZ_vtx_BS_FSR;   //!
   TBranch        *b_massZErr_vtx_BS;   //!
   TBranch        *b_massZErr_vtx_BS_FSR;   //!
   TBranch        *b_massZ_vtxChi2_BS;   //!
   TBranch        *b_PvX1;   //!
   TBranch        *b_PvX2;   //!
   TBranch        *b_PvY1;   //!
   TBranch        *b_PvY2;   //!
   TBranch        *b_PvZ1;   //!
   TBranch        *b_PvZ2;   //!
   TBranch        *b_BsX1;   //!
   TBranch        *b_BsX2;   //!
   TBranch        *b_BsY1;   //!
   TBranch        *b_BsY2;   //!
   TBranch        *b_BsZ1;   //!
   TBranch        *b_BsZ2;   //!
   TBranch        *b_muonPV_x1;   //!
   TBranch        *b_muonPV_y1;   //!
   TBranch        *b_muonPV_z1;   //!
   TBranch        *b_muonPV_x2;   //!
   TBranch        *b_muonPV_y2;   //!
   TBranch        *b_muonPV_z2;   //!
   TBranch        *b_pT1_FromMuonBestTrack;   //!
   TBranch        *b_pT2_FromMuonBestTrack;   //!
   TBranch        *b_eta1_FromMuonBestTrack;   //!
   TBranch        *b_eta2_FromMuonBestTrack;   //!
   TBranch        *b_phi1_FromMuonBestTrack;   //!
   TBranch        *b_phi2_FromMuonBestTrack;   //!
   TBranch        *b_pT1_genFromReco;   //!
   TBranch        *b_pT2_genFromReco;   //!
   TBranch        *b_pT1;   //!
   TBranch        *b_pT2;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_m1;   //!
   TBranch        *b_m2;   //!
   TBranch        *b_pT_FSR1;   //!
   TBranch        *b_pT_FSR2;   //!
   TBranch        *b_eta_FSR1;   //!
   TBranch        *b_eta_FSR2;   //!
   TBranch        *b_phi_FSR1;   //!
   TBranch        *b_phi_FSR2;   //!
   TBranch        *b_m_FSR1;   //!
   TBranch        *b_m_FSR2;   //!
   TBranch        *b_Tracker1;   //!
   TBranch        *b_Tracker2;   //!
   TBranch        *b_single_BS_pT1;   //!
   TBranch        *b_single_BS_pT2;   //!
   TBranch        *b_single_BS_eta1;   //!
   TBranch        *b_single_BS_eta2;   //!
   TBranch        *b_single_BS_phi1;   //!
   TBranch        *b_single_BS_phi2;   //!
   TBranch        *b_single_BS_m1;   //!
   TBranch        *b_single_BS_m2;   //!
   TBranch        *b_vtx_pT1;   //!
   TBranch        *b_vtx_pT2;   //!
   TBranch        *b_vtx_eta1;   //!
   TBranch        *b_vtx_eta2;   //!
   TBranch        *b_vtx_phi1;   //!
   TBranch        *b_vtx_phi2;   //!
   TBranch        *b_vtx_m1;   //!
   TBranch        *b_vtx_m2;   //!
   TBranch        *b_vtx_BS_pT1;   //!
   TBranch        *b_vtx_BS_pT2;   //!
   TBranch        *b_vtx_BS_eta1;   //!
   TBranch        *b_vtx_BS_eta2;   //!
   TBranch        *b_vtx_BS_phi1;   //!
   TBranch        *b_vtx_BS_phi2;   //!
   TBranch        *b_vtx_BS_m1;   //!
   TBranch        *b_vtx_BS_m2;   //!
   TBranch        *b_vtx_pT_FSR1;   //!
   TBranch        *b_vtx_pT_FSR2;   //!
   TBranch        *b_vtx_eta_FSR1;   //!
   TBranch        *b_vtx_eta_FSR2;   //!
   TBranch        *b_vtx_phi_FSR1;   //!
   TBranch        *b_vtx_phi_FSR2;   //!
   TBranch        *b_vtx_m_FSR1;   //!
   TBranch        *b_vtx_m_FSR2;   //!
   TBranch        *b_vtx_BS_pT_FSR1;   //!
   TBranch        *b_vtx_BS_pT_FSR2;   //!
   TBranch        *b_vtx_BS_eta_FSR1;   //!
   TBranch        *b_vtx_BS_eta_FSR2;   //!
   TBranch        *b_vtx_BS_phi_FSR1;   //!
   TBranch        *b_vtx_BS_phi_FSR2;   //!
   TBranch        *b_vtx_BS_m_FSR1;   //!
   TBranch        *b_vtx_BS_m_FSR2;   //!
   TBranch        *b_d0BS1;   //!
   TBranch        *b_d0BS2;   //!
   TBranch        *b_d0BS_vtx_BS1;   //!
   TBranch        *b_d0BS_vtx_BS2;   //!
   TBranch        *b_Id1;   //!
   TBranch        *b_Id2;   //!
   TBranch        *b_pterr1;   //!
   TBranch        *b_pterr2;   //!
   TBranch        *b_pterr1_VX;   //!
   TBranch        *b_pterr2_VX;   //!
   TBranch        *b_pterr1_VX_BS;   //!
   TBranch        *b_pterr2_VX_BS;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_nFSRPhotons;   //!

   DY_2018_macro(TTree *tree=0);
   virtual ~DY_2018_macro();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // classes  and varaibles declared further apart from base macro
  

   public:
    // variables
    double pi = TMath::Pi();    

    // to print few information regarding pT of leptons in each event and their range 
    bool debug_base = true;
    // debug_count is to see those events in which lepton pT is exceeding 200 GeV 
    bool debug_count = true;
    // analyze_pT is to just to make sure we are analyzing pT and not eta and phi, can remove it now 
    bool analyze_pT = true;
    // smearing if you are smearing in the program, otherwise it wont smear
    // debug pT is to store the pT values in each bins just as a sanity check
    bool debug_pT = true; 
    bool debug_eta = true; 
    // debug entries is to count and display number of entries in each pT bin with the bollean debug_pT 
    bool debug_entries = true; 

    TLorentzVector genLep1, genLep2, genZ;
    TLorentzVector Lep1, Lep2, recoZ;
    double massZ_reconstruct;
    double massZ_reco;
    double mass_mu = 0.105658;
    TString saving_path = "/afs/cern.ch/user/n/nrawal/work/UL_sample_generation_with_without_rochester/analysis_DYsamples_with_different_variables/2D_plots/";
    //double pt_list_mean[16] = {7.5,22.5,35,45,55,65,82.5,100,120,165};
//    double pt_list[17] = {0,20,30,35,40,42,44,46,48,50,55,60,70,90,120,150,200};
//    TString pt_list_symbol[17] = {"0","20","30","35","40","42","44","46","48","50","55","60","70","90","120","150","200"};
//    TString bin_number_pt[16] = {"1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th","11th","12th","13th","14th","15th","16th"};

    double pt_list[14] = {0,10, 15, 20, 25,31,35,40,43,46,50,60,90,200};
    TString pt_list_symbol[14] ={"0","10","15", "20", "25","31","35","40","43","46","50","60","90","200"};  
    TString bin_number_pt[13] = {"1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th","11th", "12th","13th"};

    double eta_list[15] = {-2.4, -2.1, -1.9, -1.6 , -1.2, -0.8, -0.4,  0.0, 0.4, 0.8, 1.2, 1.6, 1.9, 2.1, 2.4};
    TString eta_list_symbol[15] ={"-2.4", "-2.1", "-1.9", "-1.6", "-1.2", "-0.8", "-0.4", "0.0", "0.4", "0.8", "1.2", "1.6", "1.9", "2.1", "2.4"};  
    TString bin_number_eta[15] = {"1st","2nd","3rd","4th", "5th", "6th", "7th", "8th", "9th", "10th", "11th", "12th", "13th", "14th", "15th"};

    double counter=0;
    double counter_first=0;
    double counter_second=0;
    // to count the number of entries in each pT bin
    double n_positive_bin[pT_bins] = {0};
    double n_negative_bin[pT_bins] ={0};
    // counting total number of entries
    double total_entries_positive, total_entries_negative;   
    
    // to count the number of entries in all events with pT above 200 GeV
    double count_positive_mu = 0 , count_negative_mu = 0; 
    TH1D * hist_gen_pT1 = new TH1D("hist_gen_pT1", " Gen pt of first lepton " , 800,0,200);
    TH1D * hist_gen_pT2 = new TH1D("hist_gen_pT2", " Gen pt of second lepton " ,800,0,200);
    TH1D * hist_gen_eta1 = new TH1D("hist_gen_eta1", " Gen eta of first lepton " , 50,-2.5,2.5);
    TH1D * hist_gen_eta2 = new TH1D("hist_gen_eta2", " Gen eta of second lepton " ,50,-2.5,2.5);
    TH1D * hist_gen_phi1 = new TH1D("hist_gen_phi1", " Gen phi of first lepton " , 50,-pi,pi);
    TH1D * hist_gen_phi2 = new TH1D("hist_gen_phi2", " Gen phi of second lepton " ,50,-pi,pi);

    // histograms for reco pT1 and reco pT2
    TH1D * hist_reco_pT1 = new TH1D("hist_reco_pT1", " Reco pt of first lepton " , 800,0,200);
    TH1D * hist_reco_pT2 = new TH1D("hist_reco_pT2", " Reco pt of second lepton " ,800,0,200);
    TH1D * hist_reco_eta1 = new TH1D("hist_reco_eta1", " Reco eta of first lepton " , 50,-2.5,2.5);
    TH1D * hist_reco_eta2 = new TH1D("hist_reco_eta2", " Reco eta of second lepton " ,50,-2.5,2.5);
    TH1D * hist_reco_phi1 = new TH1D("hist_reco_phi1", " Reco phi of first lepton " , 50,-pi,pi);
    TH1D * hist_reco_phi2 = new TH1D("hist_reco_phi2", " Reco phi of second lepton " ,50,-pi,pi);


    TH1D * hist_id1 = new TH1D("hist_id1", " id  of first lepton " , 30,-15,15);
    TH1D * hist_id2 = new TH1D("hist_id2", " id  of second lepton " , 30,-15,15);
    TH1D * hist_gen_Zmass = new TH1D("hist_gen_Zmass"," Gen Z mass",200,80,100); 
    
    TH1D * hist_gen_positive_mu = new TH1D("hist_gen_positive_mu", " Gen pt (#mu^{+}) " ,800,0,200);
    TH1D * hist_id_positive_mu = new TH1D("hist_id_positive_mu", " id (#mu^{+}) " , 15,-15,0);

    TH1D * hist_gen_positive_eta = new TH1D("hist_gen_positive_eta", " Gen eta of (#mu^{+}) " , 50,-2.5,2.5);
    TH1D * hist_gen_negative_eta = new TH1D("hist_gen_negative_eta", " Gen eta of  (#mu^{-}) " ,50,-2.5,2.5);
    TH1D * hist_gen_positive_phi = new TH1D("hist_gen_positive_phi", " Gen phi of  (#mu^{+})  " , 50,-pi,pi);
    TH1D * hist_gen_negative_phi = new TH1D("hist_gen_negative_phi", " Gen phi of  (#mu^{-})  " ,50,-pi,pi);
 
    TH1D * hist_gen_negative_mu = new TH1D("hist_gen_negative_mu", " Gen pt (#mu^{-}) " , 800,0,200);
    TH1D * hist_id_negative_mu = new TH1D("hist_id_negative_mu", " id (#mu^{-}) " , 15,0,15);
  
    TH1D * hist_reco_Zmass = new TH1D("hist_reco_Zmass", " m(#mu^{-} #mu^{-}) after smearing " , 200,80,100);

   // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
    TH1D *pT_positive_bin[pT_bins][eta_bins] = {nullptr}; 
    TH1D *pT_negative_bin[pT_bins][eta_bins] = {nullptr} ; 

     // reco distribution for eta bins
    TH1D *pT_positive_bin_reco[pT_bins][eta_bins] = {nullptr} ; 
    TH1D *pT_negative_bin_reco[pT_bins][eta_bins] = {nullptr} ; 

    // smeared pT distributions in each bin
    // deltapT/pT distribution (i) reco wrt gen (ii) sneared wrt gen
    TH1D *delta_pT_positive_bin_reco[pT_bins][eta_bins] = {nullptr} ;
    TH1D *delta_pT_negative_bin_reco[pT_bins][eta_bins] = {nullptr} ;

    TH2D *delta_pT_positive_bin_with_pT;
    TH2D *delta_pT_negative_bin_with_pT;

    void histogram_deltapT_declaration(); 
    // to declare the histograms for pT distributions 
    void histogram_pT_distribution();
    
    // to printout the values of mean of pT distribution before and after smearing
    // to draw and save the all the  distributions
    void saving_histogram();

    // to reconstruct Z mass from the gen Lep instead of taking directly the variable genZ mass ; just done as a sanity check 
    void reconstruct_Zmass() ;
    // both total entries and plot entries functions are to count and display total entries in each pT bins and those events that exceeds pT 200 GeV
    //  used when debug_pT , debug_count and debug_entries are true
    void total_entries(); 
    //  void plot_entries(); 
   float delta_pT_value;
   void delta_pT_positive_bin_reco_first(int i, int j);
   void delta_pT_negative_bin_reco_first(int i, int j);
   void delta_pT_positive_bin_reco_second(int i, int j);
   void delta_pT_negative_bin_reco_second(int i, int j);
 
};

    void DY_2018_macro :: delta_pT_positive_bin_reco_first(int i, int j){
       delta_pT_value = ((pT1- pT1_genFromReco)/pT1_genFromReco);
       delta_pT_positive_bin_reco[i][j]->Fill(delta_pT_value);
       delta_pT_positive_bin_with_pT->Fill(pT1_genFromReco,pterr1/pT1_genFromReco);
    }
 
    void DY_2018_macro :: delta_pT_negative_bin_reco_first(int i, int j){
       delta_pT_value = ((pT2- pT2_genFromReco)/pT2_genFromReco);
       delta_pT_negative_bin_reco[i][j]->Fill(delta_pT_value);
       delta_pT_negative_bin_with_pT->Fill(pT2_genFromReco,pterr2/pT2_genFromReco);
    }
 
    void DY_2018_macro :: delta_pT_positive_bin_reco_second(int i, int j){
       delta_pT_value = ((pT2- pT2_genFromReco)/pT2_genFromReco);
       delta_pT_positive_bin_reco[i][j]->Fill(delta_pT_value);
       delta_pT_positive_bin_with_pT->Fill(pT2_genFromReco,pterr2/pT2_genFromReco);
    }
 
    void DY_2018_macro :: delta_pT_negative_bin_reco_second(int i, int j){
       delta_pT_value = ((pT1- pT1_genFromReco)/pT1_genFromReco);
       delta_pT_negative_bin_reco[i][j]->Fill(delta_pT_value);
       delta_pT_negative_bin_with_pT->Fill(pT1_genFromReco,pterr1/pT1_genFromReco);
    }
  
void DY_2018_macro :: histogram_deltapT_declaration(){

      TString title_p, title_n ;

      title_p = "#frac{#Delta pT}{pT} : #mu^{+}";
      delta_pT_positive_bin_with_pT = new TH2D("delta_pT_positive_bin_with_pT", title_p, 100, 0, 200, 100,0,0.2);
      delta_pT_positive_bin_with_pT->GetXaxis()->SetTitle(" pT(GeV)");
      delta_pT_positive_bin_with_pT->GetYaxis()->SetTitle("#frac{#Delta pT}{pT}");

      title_n = "#frac{#Delta pT}{pT} : #mu^{-}";
      delta_pT_negative_bin_with_pT = new TH2D("delta_pT_negative_bin_with_pT", title_n, 100, 0, 200, 100,0,0.2);
      delta_pT_negative_bin_with_pT->GetXaxis()->SetTitle(" pT(GeV)");
      delta_pT_negative_bin_with_pT->GetYaxis()->SetTitle("#frac{#Delta pT}{pT}");


      TString histogram_name_positive_pt, title_positive;
      TString histogram_name_negative_pt, title_negative;
      for(int i= 0;i<pT_bins; i++){
      for(int j= 0;j<eta_bins; j++){
        histogram_name_positive_pt = TString::Format("delta_pT_positive_bin_reco[%d][%d]",i,j);
        title_positive = "#frac{#Delta pT}{pT} (#mu^{+}) : "+pt_list_symbol[i] + " #leq pT < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        delta_pT_positive_bin_reco[i][j] = new TH1D(histogram_name_positive_pt,title_positive, 400,-0.2,0.2);
      }
     }

     for(int i= 0;i<pT_bins; i++){
      for(int j= 0;j<eta_bins; j++){
        histogram_name_negative_pt = TString::Format("delta_pT_negative_bin_reco[%d][%d]",i,j);
        title_negative = "#frac{#Delta pT}{pT} (#mu^{-}) : "+pt_list_symbol[i] + " #leq pT < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        delta_pT_negative_bin_reco[i][j] = new TH1D(histogram_name_negative_pt,title_negative, 400,-0.2,0.2);
      }
     }

 }      
      
    void DY_2018_macro :: histogram_pT_distribution(){
      TString histogram_name_positive_pT, title_positive;
      TString histogram_name_negative_pT, title_negative;
      TString histogram_name_positive_pT_reco, title_positive_reco;
      TString histogram_name_negative_pT_reco, title_negative_reco;
      float pT_limit[pT_bins]= {10,10,10,10,10,10,10,10,10,20} ;
      int nbins[pT_bins];
       for(int i=0; i<pT_bins;i++){
         nbins[i] = 2*(pt_list[i+1] - pt_list[i]+2*pT_limit[i]) * 5;
       }
     // for positive muons
     for(int i= 0;i<pT_bins; i++){
      for(int j= 0;j<eta_bins; j++){
      histogram_name_positive_pT = TString::Format("pT_positive_bin[%d][%d]",i,j);
      title_positive = "pT -  #mu^{+} :" + pt_list_symbol[i] +" #leq GEN pT < " + pt_list_symbol[i+1] + " GeV , " + eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1]  ; 
      pT_positive_bin[i][j] = new TH1D(histogram_name_positive_pT,title_positive, nbins[i],pt_list[i]-pT_limit[i],pt_list[i+1]+pT_limit[i]);
      pT_positive_bin[i][j]->GetYaxis()->SetTitle("N/100 MeV");
      }
     }
     // for negative muons
     for(int i= 0;i<pT_bins; i++){
      for(int j= 0;j<eta_bins; j++){
       histogram_name_negative_pT = TString::Format("pT_negative_bin[%d][%d]",i, j);
       title_negative = "pT -  #mu^{-} :" + pt_list_symbol[i] +" #leq GEN pT < " + pt_list_symbol[i+1] + " GeV , "+ eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1]  ; 
       pT_negative_bin[i][j] = new TH1D(histogram_name_negative_pT,title_negative, nbins[i],pt_list[i]-pT_limit[i],pt_list[i+1]+pT_limit[i]);
      pT_negative_bin[i][j]->GetYaxis()->SetTitle("N/100 MeV");
     }
    }

    // reco pT distribution in each bin
    for(int i= 0;i<pT_bins; i++){
      for(int j= 0;j<eta_bins; j++){
      histogram_name_positive_pT_reco = TString::Format("pT_positive_bin_reco[%d][%d]",i, j);
      title_positive_reco = "pT(reco) -  #mu^{+} :" + pt_list_symbol[i] +" #leq GEN pT < " + pt_list_symbol[i+1] + " GeV , " + eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1]  ; 
      pT_positive_bin_reco[i][j] = new TH1D(histogram_name_positive_pT_reco,title_positive_reco, nbins[i],pt_list[i]-pT_limit[i],pt_list[i+1]+pT_limit[i]);
      pT_positive_bin_reco[i][j]->GetYaxis()->SetTitle("N/100 MeV");
     }
    }

    // reco pT distribution in each bin for negative muons
    for(int i= 0;i<pT_bins; i++){
      for(int j= 0;j<eta_bins; j++){
      histogram_name_negative_pT_reco = TString::Format("pT_negative_bin_reco[%d][%d]",i, j);
      title_negative_reco = "pT(reco) -  #mu^{-} : " + pt_list_symbol[i] +" #leq GEN pT < " + pt_list_symbol[i+1] + " GeV , " + eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1]; 
      pT_negative_bin_reco[i][j] = new TH1D(histogram_name_negative_pT_reco,title_negative_reco, nbins[i],pt_list[i]-pT_limit[i],pt_list[i+1]+pT_limit[i]);
      pT_negative_bin_reco[i][j]->GetYaxis()->SetTitle("N/100 MeV");
     }
    }
   
  }
 
void DY_2018_macro :: reconstruct_Zmass(){
   genLep1.SetPtEtaPhiM(genLep_pt1, genLep_eta1, genLep_phi1, mass_mu);
   genLep2.SetPtEtaPhiM(genLep_pt2, genLep_eta2, genLep_phi2, mass_mu);

   Lep1.SetPtEtaPhiM(pT1, eta1, phi1, mass_mu);
   Lep2.SetPtEtaPhiM(pT2, eta2, phi2, mass_mu);

   genZ = genLep1 + genLep2;
   recoZ = Lep1+Lep2; 
   massZ_reconstruct = genZ.M();
   massZ_reco = recoZ.M();
}
void DY_2018_macro :: total_entries(){
  total_entries_positive = 0; 
  total_entries_negative = 0; 
  for(int i=0; i<pT_bins;i++){
  total_entries_positive = n_positive_bin[i] + total_entries_positive;
  total_entries_negative = n_negative_bin[i] + total_entries_negative;
    std::cout<<"total number of entries #mu^{+} in a bin  "<<i<<" entries "<<n_positive_bin[i]<<std::endl;
    std::cout<<"total number of entries #mu^{-} in a bin  "<<i<<" entries "<<n_negative_bin[i]<<std::endl;
  }
   if(debug_entries == true ){
    std::cout<<"total number of entries #mu^{+} "<<total_entries_positive<<std::endl;
    std::cout<<"total number of entries #mu^{-} "<<total_entries_negative<<std::endl;
   }
 }

/*void DY_2018_macro :: plot_entries(){
  TGraph *gr_positive_entries = new TGraph(bins,pt_list_mean,n_positive_bin);
  TGraph *gr_negative_entries = new TGraph(bins,pt_list_mean,n_negative_bin);

   gr_positive_entries->GetYaxis()->SetTitle("Entries");
   gr_positive_entries->GetXaxis()->SetTitle("pT (GeV)");
   gr_positive_entries->SetTitle("#mu^{+} before smearing");
 
   gr_negative_entries->GetYaxis()->SetTitle("Entries");
   gr_negative_entries->GetXaxis()->SetTitle("pT (GeV)");
   gr_negative_entries->SetTitle("#mu^{-} before smearing");
  TCanvas *canvas_entries_positive = new TCanvas("canvas_entries_positive","No. of entries in pT bin : #mu^{+}");
  canvas_entries_positive->cd();
  gr_positive_entries->Draw("AP*");

  canvas_entries_positive->SaveAs(saving_path+"entries_positive.pdf");

  TCanvas *canvas_entries_negative = new TCanvas("canvas_entries_negative","No. of entries in pT bin : #mu^{-}");
  canvas_entries_negative->cd();
  gr_negative_entries->Draw("AP*");
  canvas_entries_negative->SaveAs(saving_path+"entries_negative.pdf");
}*/

  void DY_2018_macro:: saving_histogram(){
    TFile *file_output = new TFile("DY_ntuple_without_rochester_correction_resolution_etabins.root","RECREATE"); 
    
     hist_gen_pT1->Write(); 
     hist_gen_pT2->Write(); 
     hist_id1->Write(); 
     hist_id2->Write(); 
     hist_gen_Zmass->Write(); 
    
     hist_gen_positive_mu->Write(); 
     hist_id_positive_mu->Write(); 
     hist_gen_positive_eta->Write(); 
     hist_gen_positive_phi->Write(); 

     hist_gen_negative_mu->Write(); 
     hist_id_negative_mu->Write(); 
     hist_gen_negative_eta->Write(); 
     hist_gen_negative_phi->Write(); 
 
     delta_pT_positive_bin_with_pT->Write();
     delta_pT_negative_bin_with_pT->Write();
     hist_reco_Zmass->Write();
   
    for(int i=0; i<pT_bins; i++){ 
    for(int j=0; j<eta_bins; j++){ 
    pT_positive_bin[i][j]->Write() ; 
    pT_negative_bin[i][j]->Write() ; 
    pT_positive_bin_reco[i][j]->Write() ; 
    pT_negative_bin_reco[i][j]->Write() ; 
    delta_pT_positive_bin_reco[i][j]->Write();
    delta_pT_negative_bin_reco[i][j]->Write();
    }
    }
    file_output->Close();
  }
#endif

#ifdef DY_2018_macro_cxx
DY_2018_macro::DY_2018_macro(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/user/f/ferrico/www/Chenguang/DY_muons_UL_20_s5m0/DY_amcatnlo_2018_skimmed.root");
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/home-c/caruta/HZZ_SkimmedNtuples/DY_madgraph_skimmed_DYJetsToLL_2018.root");
      if (!f || !f->IsOpen()) {

        std::cout<<"**************************************************"<<std::endl;
         f = new TFile("/eos/home-c/caruta/HZZ_SkimmedNtuples/DY_madgraph_skimmed_DYJetsToLL_2018.root");
         //f = new TFile("/eos/user/f/ferrico/www/Chenguang/DY_muons_UL_20_s5m0/DY_amcatnlo_2018_skimmed.root");
      }
      f->GetObject("passedEvents",tree);

   }
   Init(tree);

}

DY_2018_macro::~DY_2018_macro()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DY_2018_macro::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DY_2018_macro::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DY_2018_macro::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   cout<<"started"<<endl;
   fChain->SetBranchAddress("GENmass2l", &GENmass2l, &b_GENmass2l);
   fChain->SetBranchAddress("genLep_pt1", &genLep_pt1, &b_genLep_pt1);
   fChain->SetBranchAddress("genLep_pt2", &genLep_pt2, &b_genLep_pt2);
   fChain->SetBranchAddress("genLep_eta1", &genLep_eta1, &b_genLep_eta1);
   fChain->SetBranchAddress("genLep_eta2", &genLep_eta2, &b_genLep_eta2);
   fChain->SetBranchAddress("genLep_phi1", &genLep_phi1, &b_genLep_phi1);
   fChain->SetBranchAddress("genLep_phi2", &genLep_phi2, &b_genLep_phi2);
   fChain->SetBranchAddress("NInt", &NInt, &b_NInt);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
//   fChain->SetBranchAddress("PvX", &PvX, &b_PvX);
//   fChain->SetBranchAddress("PvY", &PvY, &b_PvY);
//   fChain->SetBranchAddress("PvZ", &PvZ, &b_PvZ);
//   fChain->SetBranchAddress("BsX", &BsX, &b_BsX);
//   fChain->SetBranchAddress("BsY", &BsY, &b_BsY);
//   fChain->SetBranchAddress("BsZ", &BsZ, &b_BsZ);
   fChain->SetBranchAddress("massZ", &massZ, &b_massZ);
   fChain->SetBranchAddress("massZErr", &massZErr, &b_massZErr);
   fChain->SetBranchAddress("massZ_FSR", &massZ_FSR, &b_massZ_FSR);
   fChain->SetBranchAddress("massZErr_FSR", &massZErr_FSR, &b_massZErr_FSR);
 //  fChain->SetBranchAddress("massZ_single_BS", &massZ_single_BS, &b_massZ_single_BS);
 //  fChain->SetBranchAddress("massZ_vtx", &massZ_vtx, &b_massZ_vtx);
 //  fChain->SetBranchAddress("massZ_vtx_FSR", &massZ_vtx_FSR, &b_massZ_vtx_FSR);
 //  fChain->SetBranchAddress("massZErr_vtx", &massZErr_vtx, &b_massZErr_vtx);
 //  fChain->SetBranchAddress("massZErr_vtx_FSR", &massZErr_vtx_FSR, &b_massZErr_vtx_FSR);
 //  fChain->SetBranchAddress("massZ_vtxChi2", &massZ_vtxChi2, &b_massZ_vtxChi2);
   fChain->SetBranchAddress("massZ_vtx_BS", &massZ_vtx_BS, &b_massZ_vtx_BS);
   fChain->SetBranchAddress("massZ_vtx_BS_FSR", &massZ_vtx_BS_FSR, &b_massZ_vtx_BS_FSR);
   fChain->SetBranchAddress("massZErr_vtx_BS", &massZErr_vtx_BS, &b_massZErr_vtx_BS);
   fChain->SetBranchAddress("massZErr_vtx_BS_FSR", &massZErr_vtx_BS_FSR, &b_massZErr_vtx_BS_FSR);
   fChain->SetBranchAddress("massZ_vtxChi2_BS", &massZ_vtxChi2_BS, &b_massZ_vtxChi2_BS);
 //  fChain->SetBranchAddress("PvX1", &PvX1, &b_PvX1);
 //  fChain->SetBranchAddress("PvX2", &PvX2, &b_PvX2);
 //  fChain->SetBranchAddress("PvY1", &PvY1, &b_PvY1);
 //  fChain->SetBranchAddress("PvY2", &PvY2, &b_PvY2);
 //  fChain->SetBranchAddress("PvZ1", &PvZ1, &b_PvZ1);
 //  fChain->SetBranchAddress("PvZ2", &PvZ2, &b_PvZ2);
 //  fChain->SetBranchAddress("BsX1", &BsX1, &b_BsX1);
 //  fChain->SetBranchAddress("BsX2", &BsX2, &b_BsX2);
 //  fChain->SetBranchAddress("BsY1", &BsY1, &b_BsY1);
 //  fChain->SetBranchAddress("BsY2", &BsY2, &b_BsY2);
 //  fChain->SetBranchAddress("BsZ1", &BsZ1, &b_BsZ1);
 //  fChain->SetBranchAddress("BsZ2", &BsZ2, &b_BsZ2);
  // fChain->SetBranchAddress("muonPV_x1", &muonPV_x1, &b_muonPV_x1);
  // fChain->SetBranchAddress("muonPV_y1", &muonPV_y1, &b_muonPV_y1);
  // fChain->SetBranchAddress("muonPV_z1", &muonPV_z1, &b_muonPV_z1);
  // fChain->SetBranchAddress("muonPV_x2", &muonPV_x2, &b_muonPV_x2);
  // fChain->SetBranchAddress("muonPV_y2", &muonPV_y2, &b_muonPV_y2);
  // fChain->SetBranchAddress("muonPV_z2", &muonPV_z2, &b_muonPV_z2);
   fChain->SetBranchAddress("pT1_FromMuonBestTrack", &pT1_FromMuonBestTrack, &b_pT1_FromMuonBestTrack);
   fChain->SetBranchAddress("pT2_FromMuonBestTrack", &pT2_FromMuonBestTrack, &b_pT2_FromMuonBestTrack);
   fChain->SetBranchAddress("eta1_FromMuonBestTrack", &eta1_FromMuonBestTrack, &b_eta1_FromMuonBestTrack);
   fChain->SetBranchAddress("eta2_FromMuonBestTrack", &eta2_FromMuonBestTrack, &b_eta2_FromMuonBestTrack);
   fChain->SetBranchAddress("phi1_FromMuonBestTrack", &phi1_FromMuonBestTrack, &b_phi1_FromMuonBestTrack);
   fChain->SetBranchAddress("phi2_FromMuonBestTrack", &phi2_FromMuonBestTrack, &b_phi2_FromMuonBestTrack);
   fChain->SetBranchAddress("pT1_genFromReco", &pT1_genFromReco, &b_pT1_genFromReco);
   fChain->SetBranchAddress("pT2_genFromReco", &pT2_genFromReco, &b_pT2_genFromReco);
   fChain->SetBranchAddress("pT1", &pT1, &b_pT1);
   fChain->SetBranchAddress("pT2", &pT2, &b_pT2);
   fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("m1", &m1, &b_m1);
   fChain->SetBranchAddress("m2", &m2, &b_m2);
   fChain->SetBranchAddress("pT_FSR1", &pT_FSR1, &b_pT_FSR1);
   fChain->SetBranchAddress("pT_FSR2", &pT_FSR2, &b_pT_FSR2);
   fChain->SetBranchAddress("eta_FSR1", &eta_FSR1, &b_eta_FSR1);
   fChain->SetBranchAddress("eta_FSR2", &eta_FSR2, &b_eta_FSR2);
   fChain->SetBranchAddress("phi_FSR1", &phi_FSR1, &b_phi_FSR1);
   fChain->SetBranchAddress("phi_FSR2", &phi_FSR2, &b_phi_FSR2);
   fChain->SetBranchAddress("m_FSR1", &m_FSR1, &b_m_FSR1);
   fChain->SetBranchAddress("m_FSR2", &m_FSR2, &b_m_FSR2);
   fChain->SetBranchAddress("Tracker1", &Tracker1, &b_Tracker1);
   fChain->SetBranchAddress("Tracker2", &Tracker2, &b_Tracker2);
//   fChain->SetBranchAddress("single_BS_pT1", &single_BS_pT1, &b_single_BS_pT1);
//   fChain->SetBranchAddress("single_BS_pT2", &single_BS_pT2, &b_single_BS_pT2);
//   fChain->SetBranchAddress("single_BS_eta1", &single_BS_eta1, &b_single_BS_eta1);
//   fChain->SetBranchAddress("single_BS_eta2", &single_BS_eta2, &b_single_BS_eta2);
//   fChain->SetBranchAddress("single_BS_phi1", &single_BS_phi1, &b_single_BS_phi1);
//   fChain->SetBranchAddress("single_BS_phi2", &single_BS_phi2, &b_single_BS_phi2);
//   fChain->SetBranchAddress("single_BS_m1", &single_BS_m1, &b_single_BS_m1);
//   fChain->SetBranchAddress("single_BS_m2", &single_BS_m2, &b_single_BS_m2);
//   fChain->SetBranchAddress("vtx_pT1", &vtx_pT1, &b_vtx_pT1);
//   fChain->SetBranchAddress("vtx_pT2", &vtx_pT2, &b_vtx_pT2);
//   fChain->SetBranchAddress("vtx_eta1", &vtx_eta1, &b_vtx_eta1);
//   fChain->SetBranchAddress("vtx_eta2", &vtx_eta2, &b_vtx_eta2);
//   fChain->SetBranchAddress("vtx_phi1", &vtx_phi1, &b_vtx_phi1);
//   fChain->SetBranchAddress("vtx_phi2", &vtx_phi2, &b_vtx_phi2);
//   fChain->SetBranchAddress("vtx_m1", &vtx_m1, &b_vtx_m1);
//   fChain->SetBranchAddress("vtx_m2", &vtx_m2, &b_vtx_m2);
   fChain->SetBranchAddress("vtx_BS_pT1", &vtx_BS_pT1, &b_vtx_BS_pT1);
   fChain->SetBranchAddress("vtx_BS_pT2", &vtx_BS_pT2, &b_vtx_BS_pT2);
   fChain->SetBranchAddress("vtx_BS_eta1", &vtx_BS_eta1, &b_vtx_BS_eta1);
   fChain->SetBranchAddress("vtx_BS_eta2", &vtx_BS_eta2, &b_vtx_BS_eta2);
   fChain->SetBranchAddress("vtx_BS_phi1", &vtx_BS_phi1, &b_vtx_BS_phi1);
   fChain->SetBranchAddress("vtx_BS_phi2", &vtx_BS_phi2, &b_vtx_BS_phi2);
   fChain->SetBranchAddress("vtx_BS_m1", &vtx_BS_m1, &b_vtx_BS_m1);
   fChain->SetBranchAddress("vtx_BS_m2", &vtx_BS_m2, &b_vtx_BS_m2);
//   fChain->SetBranchAddress("vtx_pT_FSR1", &vtx_pT_FSR1, &b_vtx_pT_FSR1);
//   fChain->SetBranchAddress("vtx_pT_FSR2", &vtx_pT_FSR2, &b_vtx_pT_FSR2);
//   fChain->SetBranchAddress("vtx_eta_FSR1", &vtx_eta_FSR1, &b_vtx_eta_FSR1);
//   fChain->SetBranchAddress("vtx_eta_FSR2", &vtx_eta_FSR2, &b_vtx_eta_FSR2);
//   fChain->SetBranchAddress("vtx_phi_FSR1", &vtx_phi_FSR1, &b_vtx_phi_FSR1);
//   fChain->SetBranchAddress("vtx_phi_FSR2", &vtx_phi_FSR2, &b_vtx_phi_FSR2);
//   fChain->SetBranchAddress("vtx_m_FSR1", &vtx_m_FSR1, &b_vtx_m_FSR1);
//   fChain->SetBranchAddress("vtx_m_FSR2", &vtx_m_FSR2, &b_vtx_m_FSR2);
   fChain->SetBranchAddress("vtx_BS_pT_FSR1", &vtx_BS_pT_FSR1, &b_vtx_BS_pT_FSR1);
   fChain->SetBranchAddress("vtx_BS_pT_FSR2", &vtx_BS_pT_FSR2, &b_vtx_BS_pT_FSR2);
   fChain->SetBranchAddress("vtx_BS_eta_FSR1", &vtx_BS_eta_FSR1, &b_vtx_BS_eta_FSR1);
   fChain->SetBranchAddress("vtx_BS_eta_FSR2", &vtx_BS_eta_FSR2, &b_vtx_BS_eta_FSR2);
   fChain->SetBranchAddress("vtx_BS_phi_FSR1", &vtx_BS_phi_FSR1, &b_vtx_BS_phi_FSR1);
   fChain->SetBranchAddress("vtx_BS_phi_FSR2", &vtx_BS_phi_FSR2, &b_vtx_BS_phi_FSR2);
   fChain->SetBranchAddress("vtx_BS_m_FSR1", &vtx_BS_m_FSR1, &b_vtx_BS_m_FSR1);
   fChain->SetBranchAddress("vtx_BS_m_FSR2", &vtx_BS_m_FSR2, &b_vtx_BS_m_FSR2);
//   fChain->SetBranchAddress("d0BS1", &d0BS1, &b_d0BS1);
//   fChain->SetBranchAddress("d0BS2", &d0BS2, &b_d0BS2);
//   fChain->SetBranchAddress("d0BS_vtx_BS1", &d0BS_vtx_BS1, &b_d0BS_vtx_BS1);
//   fChain->SetBranchAddress("d0BS_vtx_BS2", &d0BS_vtx_BS2, &b_d0BS_vtx_BS2);
   fChain->SetBranchAddress("Id1", &Id1, &b_Id1);
   fChain->SetBranchAddress("Id2", &Id2, &b_Id2);
   fChain->SetBranchAddress("pterr1", &pterr1, &b_pterr1);
   fChain->SetBranchAddress("pterr2", &pterr2, &b_pterr2);
//   fChain->SetBranchAddress("pterr1_VX", &pterr1_VX, &b_pterr1_VX);
//   fChain->SetBranchAddress("pterr2_VX", &pterr2_VX, &b_pterr2_VX);
   fChain->SetBranchAddress("pterr1_VX_BS", &pterr1_VX_BS, &b_pterr1_VX_BS);
   fChain->SetBranchAddress("pterr2_VX_BS", &pterr2_VX_BS, &b_pterr2_VX_BS);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("nFSRPhotons", &nFSRPhotons, &b_nFSRPhotons);
   Notify();
}

Bool_t DY_2018_macro::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DY_2018_macro::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DY_2018_macro::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DY_2018_macro_cxx
