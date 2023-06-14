// this file is just to analyze those produced histograms by the DY smaples
//
#include<iostream>
using namespace std;
#include <TROOT.h>  
#include <TChain.h>
#include <TFile.h> 
// Header file for the classes stored in the TTree if any.                                                                                                                          
#include "string" 
#include "vector" 
#include "RooRealVar.h"   
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h" 
#include "RooHistPdf.h" 
#include "RooPolynomial.h"
#include "RooAbsArg.h"   
#include "RooPlot.h"    
//#include "TRatioPlot.h" 
#include "RooAddPdf.h" 
#include "RooFitResult.h"
#include "TAxis.h"      
#include "TH1.h"                                                                                                                                                                    
using namespace RooFit ;
const int pT_bins=13;
const int eta_bins=14;
// base class having all the functions to analyze the histograms
// The functions defined in base class are called by derived class and used in that way
class class_reading{
	public:
   TFile *file = TFile::Open("/afs/cern.ch/user/n/nrawal/work/UL_sample_generation_with_without_rochester/analysis_DYsamples_with_different_variables/2D_plots/DY_ntuple_without_rochester_correction_resolution_etabins.root"); 
   TString saving_path = "/afs/cern.ch/user/n/nrawal/work/UL_sample_generation_with_without_rochester/analysis_DYsamples_with_different_variables/2D_plots/results_without_rochester_correction_resolution/";    
   float pt_list[14] = {0,10,15,20, 25,31,35,40,43,46,50,60,90,200};
   float eta_list[15] = {-2.4, -2.1, -1.9, -1.6, -1.2, -0.8, -0.4, 0,0.4, 0.8, 1.2, 1.6, 1.9, 2.1, 2.4};
   TString pt_list_symbol[14] ={"0","10","15", "20", "25","31","35","40","43","46","50","60","90","200"};  
   TString eta_list_symbol[15] ={"-2.4", "-2.1", "-1.9", "-1.6", "-1.2",  "-0.8", "-0.4", "0","0.4", "0.8", "1.2", "1.6", "1.9", "2.1","2.4"};  
   TString bin_number_pt[13] = {"1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th","11th","12th","13th"};
   TString bin_number_eta[14] = {"1st","2nd","3rd","4th","5th", "6th", "7th", "8th", "9th", "10th","11th", "12th", "13th", "14th"};

//   float pt_list_mean[10] = {12.5,28,33, 38.5,42.5, 44.5,48,55,75,145};
   // float pt_list_mean[10] = {7.5,22.5,35,45,55,65,82.5,100,120,165};
   float mean_positive_delta_pT[pT_bins][eta_bins]={0};
   float mean_positive_error_delta_pT[pT_bins][eta_bins]={0};
   float sigma_positive_delta_pT[pT_bins][eta_bins]={0};
   float sigma_positive_error_delta_pT[pT_bins][eta_bins]={0};

   float mean_negative_delta_pT[pT_bins][eta_bins]={0};
   float mean_negative_error_delta_pT[pT_bins][eta_bins]={0};
   float sigma_negative_delta_pT[pT_bins][eta_bins]={0};
   float sigma_negative_error_delta_pT[pT_bins][eta_bins]={0};
   float chi2_ndf_positive_delta_pT[pT_bins][eta_bins] = {0};
   float chi2_ndf_negative_delta_pT[pT_bins][eta_bins] = {0};


	 // hisotgram declaration to read the histograms from the root file
    TH1F * histogram_gen_Zmass; 
    TH1F * histogram_gen_positive_mu ;
    TH1F * histogram_gen_positive_eta ;
    TH1F * histogram_gen_positive_phi ;
    TH1F * histogram_id_positive_mu; 
    TH1F * histogram_gen_negative_mu; 
    TH1F * histogram_gen_negative_eta ;
    TH1F * histogram_gen_negative_phi ;

    TH1F * histogram_id_negative_mu; 
//		TH1F * histogram_eta_gen_positive_mu; 
//		TH1F * histogram_eta_gen_negative_mu; 
//		TH1F * histogram_phi_gen_positive_mu; 
//		TH1F * histogram_phi_gen_negative_mu; 
    TH1F * histogram_reco_Zmass; 
//    TH1F * histogram_pt_gen_Z; 
//    TH1F * histogram_eta_gen_Z; 
//    TH1F * histogram_phi_gen_Z; 

     // this histograms are defined in a loop since its an array  
  
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   TH2F *histogram_delta_pT_positive_bin_with_pT;  
   TH2F *histogram_delta_pT_negative_bin_with_pT;  
 
    TH1F *histogram_pT_positive_bin[pT_bins][eta_bins] ; 
    TH1F *histogram_pT_negative_bin[pT_bins][eta_bins] ; 
    TH1F *histogram_pT_positive_bin_reco[pT_bins][eta_bins] ; 
    TH1F *histogram_pT_negative_bin_reco[pT_bins][eta_bins] ; 

     TH1F *histogram_delta_pT_positive_bin_reco[pT_bins][eta_bins]; 
     TH1F *histogram_delta_pT_negative_bin_reco[pT_bins][eta_bins]; 

     float pt_error_positive[pT_bins][eta_bins] = {0};
     float pt_error_negative[pT_bins][eta_bins] = {0};
     float pt_measured_positive[pT_bins][eta_bins] = {0}; 
     float pt_measured_negative[pT_bins][eta_bins] = {0}; 
//    TH1F *histogram_pT_negative_corresponding_positive_bin[9] ; 
//    TH1F *histogram_pT_positive_corresponding_negative_bin[9] ; 
//    TH1F *histogram_pT_negative_corresponding_positive_bin_smear[9] ; 
//    TH1F *histogram_pT_positive_corresponding_negative_bin_smear[9] ; 

		// delcared variables to analyze the mean and sigma of the distribution

   // to make a final plot for deltapT mean and pT 

	 // constructor to open the root file
	 void initializing();

  ofstream output_file;
   TFile * output_root_file;
  	// to plot histogram for pT distribution and mass distribution individually 
//   void plotting_hist(TH1F* hist_draw, TString title, TString saving_name, TString title_name_axis);
   void plotting_hist_2D(TH2F* hist_draw, TString title, TString saving_name, TString title_name_axis);
	 	// to count total number of entries in each pT bin and also above 200 GeV and print them in a text file myfile_n_entries                                                          
  	// these 3 functions are for sanity check
    // void plot_entries();   
    // void total_entries();   
    // void evaluating_mean_pT(); 
    // plotting mean and sigma distribution 		

    // to save all the histograms in pdf plots	 
    void saving_histogram_pT();
    void saving_histogram_delta_pT();
    void fitting_delta_pT_distribution();
		// graphs to plot mean and sigma, and difference due to smearing, and  also print them in  text files
    void plotting_mean_delta_pT();
    void temperature_plot_delta_pT();
    void fitting_DSCB_positive(TH1F * , TString saving_name, TString title_name, int, int);
    void fitting_DSCB_negative(TH1F * , TString saving_name, TString title_name, int, int);
};

// declaring other class and its variables - the one which is actually used -=> inherited class from class reading
class derived_class_reading : public class_reading{
	public:

   //TString pt_list_symbol[11] ={"0","25","31","35","40","43","46","50","60","90","200"};  

   //TString eta_list_symbol[5] ={"0","0.9","1.2","2.0","2.4"};  
//  std::pair<float,float>mean_sigma_Zmass_gauss; 
//  std::pair<float,float> mean_sigma_Zmass_gauss_smear; 
  
  void fitting_delta_pT();
  void saving_histograms();
};

void derived_class_reading :: saving_histograms(){
    saving_histogram_pT();
    saving_histogram_delta_pT();
}

void derived_class_reading :: fitting_delta_pT(){
   fitting_delta_pT_distribution();
   plotting_mean_delta_pT();
   temperature_plot_delta_pT();
}

  // declaring base class functions
	//
	// to plot mean Jpsi after smearing
  
// this is to plot mean Jpsi mass before smearing
 /*  void class_reading :: plotting_hist(TH1F * hist_draw, TString title, TString saving_name, TString title_name_axis){
    
    gStyle->SetOptStat(000001111);
    gStyle->SetStatFontSize(0.02); 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1); 

    TCanvas * c = new TCanvas("c",title,900,600);
    c->cd();
    hist_draw->GetXaxis()->SetTitle(title_name_axis+"");
    hist_draw->Draw();
    c->SaveAs(saving_path+saving_name+".pdf");
   }
*/
  void class_reading :: plotting_hist_2D(TH2F * hist_draw, TString title, TString saving_name, TString title_name_axis){
    hist_draw->SetStats(0); 
//    gStyle->SetOptStat(000001111);
//    gStyle->SetStatFontSize(0.02); 
//    gStyle->SetStatW(0.2);
//    gStyle->SetStatH(0.1); 

    TCanvas * c = new TCanvas("c",title,900,600);
    c->cd();
    c->SetRightMargin(0.2);
    //hist_draw->GetXaxis()->SetTitle(title_name_axis+"");
    hist_draw->Draw("colZ");
    hist_draw->GetXaxis()->SetTitle("pT (GeV)");
    hist_draw->GetYaxis()->SetTitle("#eta");
    c->SaveAs(saving_path+saving_name+".pdf");
   }

void class_reading :: saving_histogram_delta_pT(){
   TString string_delta_pT_distribution_negative_reco, string_delta_pT_distribution_positive_reco;  
  // TString string_delta_pT_distribution_negative_smear, string_delta_pT_distribution_positive_smear;  
   // histograms and corresponding other pT distribution
    gStyle->SetOptStat(000001111);
    gStyle->SetStatFontSize(0.02); 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1); 

  TCanvas *canvas_pT_positive_reco;
  TCanvas *canvas_pT_negative_reco;
  for(int i=0; i<pT_bins; i++){
  for(int j=0; j<eta_bins; j++){
  canvas_pT_positive_reco = new TCanvas("canvas_pT_positive_reco", "reco #frac{#Delta pT}{pT} distribution in each bin #mu^{+}", 900,600);
  canvas_pT_positive_reco->cd();
  histogram_delta_pT_positive_bin_reco[i][j]->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
  histogram_delta_pT_positive_bin_reco[i][j]->Draw();
  histogram_delta_pT_positive_bin_reco[i][j]->GetXaxis()->SetTitleOffset(1.2);
  histogram_delta_pT_positive_bin_reco[i][j]->GetXaxis()->CenterTitle(true);
  string_delta_pT_distribution_positive_reco ="delta_pT_positive_" + bin_number_pt[i] + "_eta_"+bin_number_eta[j]+"_bin_reco";
  canvas_pT_positive_reco->SaveAs(saving_path+string_delta_pT_distribution_positive_reco+".pdf");
  canvas_pT_positive_reco->Close();
 
  canvas_pT_negative_reco = new TCanvas("canvas_pT_negative_reco", "reco #frac{#Delta pT}{pT}  distribution in each bin #mu^{+}", 900,600);
  canvas_pT_negative_reco->cd();
  histogram_delta_pT_negative_bin_reco[i][j]->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
  histogram_delta_pT_negative_bin_reco[i][j]->Draw();
  histogram_delta_pT_negative_bin_reco[i][j]->GetXaxis()->SetTitleOffset(1.2);
  histogram_delta_pT_negative_bin_reco[i][j]->GetXaxis()->CenterTitle(true);
  string_delta_pT_distribution_negative_reco ="delta_pT_negative_" + bin_number_pt[i] + "_eta_"+bin_number_eta[j]+"_bin_reco";
  canvas_pT_negative_reco->SaveAs(saving_path+string_delta_pT_distribution_negative_reco+".pdf");
  canvas_pT_negative_reco->Close();
  }
 }

}
 
void class_reading :: saving_histogram_pT(){
   TString string_pT_distribution_negative, string_pT_distribution_positive;  
   TString string_pT_distribution_negative_reco, string_pT_distribution_positive_reco;  
    gStyle->SetOptStat(000001111);
    gStyle->SetStatFontSize(0.02); 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1); 
//  double range_yaxis[10] = {100000, 1000000, 1000000, 1000000, 500000, 80000, 40000, 20000, 15000, 10000}; 
   // histograms and corresponding other pT distribution
TCanvas *canvas_pT_positive;
TCanvas *canvas_pT_negative;

    std::cout<<"entered into saving histogram pT "<<std::endl;
  
  for(int i=0; i<pT_bins; i++){
  for(int j=0; j<eta_bins; j++){
    std::cout<<"entered into saving histogram pT : i  "<<i<<" j "<<j<<std::endl;
  canvas_pT_positive = new TCanvas("canvas_pT_positive", "pT distribution in each bin #mu^{+}", 900,600);
  canvas_pT_positive->cd();

  histogram_pT_positive_bin[i][j]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_positive_bin[i][j]->Draw();
  string_pT_distribution_positive ="pT_positive_" + bin_number_pt[i] + "_eta " + bin_number_eta[j] +"bin";
  canvas_pT_positive->SaveAs(saving_path+string_pT_distribution_positive+".pdf");
  canvas_pT_positive->Close();
  
    std::cout<<"could not finish into saving histogram pT "<<std::endl;
  pt_measured_positive[i][j] =  histogram_pT_positive_bin[i][j]->GetMean();
  std::cout<<"pt measured value : positive "<<pt_measured_positive[i][j]<<std::endl;
  pt_error_positive[i][j] = histogram_pT_positive_bin[i][j]->GetRMS()/histogram_pT_positive_bin[i][j]->Integral();

  canvas_pT_negative = new TCanvas("canvas_pT_negative", "pT distribution in each bin #mu^{-}", 900,600);
  canvas_pT_negative->cd();
  histogram_pT_negative_bin[i][j]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_negative_bin[i][j]->Draw();
  string_pT_distribution_negative ="pT_negative_" + bin_number_pt[i] + "_eta " + bin_number_eta[j] +"bin";
  canvas_pT_negative->SaveAs(saving_path+string_pT_distribution_negative+".pdf");
  canvas_pT_negative->Close();

  pt_measured_negative[i][j] =  histogram_pT_negative_bin[i][j]->GetMean();
  pt_error_negative[i][j] = histogram_pT_negative_bin[i][j]->GetRMS()/histogram_pT_negative_bin[i][j]->Integral();
  std::cout<<" pt values "<<pt_measured_negative[i][j]<<" error "<<pt_error_negative[i][j]<<std::endl;

  }
} 
}
void class_reading :: initializing(){

    output_file.open("mean_sigma.txt");
//initialize all the 9 histograms positive and negative one before and afer smearing
    std::cout<<"started initialization "<<std::endl;
    output_root_file = new TFile("with_rochester_histograms_delta_pT.root","RECREATE");
    histogram_gen_Zmass = (TH1F*) file->Get("hist_gen_Zmass"); 
    histogram_reco_Zmass = (TH1F*) file->Get("hist_reco_Zmass"); 
   histogram_gen_positive_mu= (TH1F*) file->Get("hist_gen_positive_mu") ;
   histogram_id_positive_mu= (TH1F*) file->Get("hist_id_positive_mu"); 
   histogram_gen_negative_mu= (TH1F*) file->Get("hist_gen_negative_mu"); 
   histogram_id_negative_mu= (TH1F*) file->Get("hist_id_negative_mu"); 

   histogram_gen_positive_eta= (TH1F*) file->Get("hist_gen_positive_eta") ;
   histogram_gen_negative_eta= (TH1F*) file->Get("hist_gen_negative_eta") ;
   histogram_gen_positive_phi= (TH1F*) file->Get("hist_gen_positive_phi") ;
   histogram_gen_negative_phi= (TH1F*) file->Get("hist_gen_negative_phi") ;
    
   histogram_delta_pT_positive_bin_with_pT = (TH2F*) file->Get("delta_pT_positive_bin_with_pT") ; 
   histogram_delta_pT_negative_bin_with_pT = (TH2F*) file->Get("delta_pT_negative_bin_with_pT") ; 
//	 histogram_eta_gen_positive_mu= (TH1F*) file->Get("hist_eta_gen_positive_mu"); 
//	 histogram_eta_gen_negative_mu= (TH1F*) file->Get("hist_eta_gen_negative_mu"); 
//	 histogram_phi_gen_positive_mu= (TH1F*) file->Get("hist_phi_gen_positive_mu"); 
//	 histogram_phi_gen_negative_mu= (TH1F*) file->Get("hist_phi_gen_negative_mu"); 
//
    std::cout<<"continuing initialization "<<std::endl;
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution

   for(int i=0; i<pT_bins; i++){
   for(int j=0; j<eta_bins; j++){
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   
		TString title_pT_positive_bin = TString::Format("pT_positive_bin[%d][%d]",i,j);
		TString title_pT_negative_bin = TString::Format("pT_negative_bin[%d][%d]",i,j);
    TString title_pT_positive_bin_reco = TString::Format("pT_positive_bin_reco[%d][%d]",i,j);
		TString title_pT_negative_bin_reco = TString::Format("pT_negative_bin_reco[%d][%d]",i,j);

    std::cout<<"print title "<<title_pT_positive_bin<<std::endl;
    // delta pT distribution
    TString title_delta_pT_positive_bin_reco = TString::Format("delta_pT_positive_bin_reco[%d][%d]",i,j);
		TString title_delta_pT_negative_bin_reco = TString::Format("delta_pT_negative_bin_reco[%d][%d]",i,j);
    histogram_pT_positive_bin[i][j] = (TH1F*) file->Get(title_pT_positive_bin) ; 
    histogram_pT_negative_bin[i][j] = (TH1F*) file->Get(title_pT_negative_bin) ; 
    histogram_pT_positive_bin_reco[i][j] = (TH1F*) file->Get(title_pT_positive_bin_reco) ; 
    histogram_pT_negative_bin_reco[i][j] = (TH1F*) file->Get(title_pT_negative_bin_reco) ; 
    histogram_delta_pT_positive_bin_reco[i][j] = (TH1F*) file->Get(title_delta_pT_positive_bin_reco) ; 
    histogram_delta_pT_negative_bin_reco[i][j] = (TH1F*) file->Get(title_delta_pT_negative_bin_reco) ; 
	}
}

std::cout<<"finished initialization "<<std::endl;
}
  
  void class_reading::fitting_delta_pT_distribution(){
    //fitting both deltapT distributions with a gaussian functions
    gStyle->SetOptFit();
    std::cout<<"declaring canvas 1 fitting deltapTfunction"<<std::endl;
    std::cout<<"inside class reading fitting delta pTfunction"<<std::endl; 
     TDirectory *dir_fit = (TDirectory*) output_root_file->mkdir("fits");
     dir_fit->cd();

    std::pair<float , float> mean_sigma_positive_delta_pT[pT_bins][eta_bins]; 
        output_file<<"positive muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "delta_pT_plots/delta_pT_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "#frac{#Delta pT}{pT} (#mu^{+}) : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        //mean_sigma_positive_delta_pT[i][j] = fitting_DSCB(histogram_delta_pT_positive_bin_reco[i][j], saving_string, title_name);
        fitting_DSCB_positive(histogram_delta_pT_positive_bin_reco[i][j], saving_string, title_name, i, j);
        output_file<<"mean : "<<mean_positive_delta_pT[i][j]<<" sigma "<<sigma_positive_delta_pT[i][j]<<"mean erro : "<<mean_positive_error_delta_pT[i][j]<<" sigma "<<sigma_positive_error_delta_pT[i][j]<<std::endl;
//        std::cout<<"mean sigma pair "<<mean_positive_delta_pT[i][j]<<" second "<<sigma_positive_delta_pT[i][j]<<std::endl;
      }
    }

        output_file<<"negative muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "delta_pT_plots/delta_pT_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "#frac{#Delta pT}{pT} (#mu^{-}) : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        //mean_sigma_positive_delta_pT[i][j] = fitting_DSCB(histogram_delta_pT_positive_bin_reco[i][j], saving_string, title_name);
        fitting_DSCB_negative(histogram_delta_pT_negative_bin_reco[i][j], saving_string, title_name, i, j);
        output_file<<"mean : "<<mean_negative_delta_pT[i][j]<<" sigma "<<sigma_negative_delta_pT[i][j]<<"mean erro : "<<mean_negative_error_delta_pT[i][j]<<" sigma "<<sigma_negative_error_delta_pT[i][j]<<std::endl;
//        std::cout<<"mean sigma pair "<<mean_positive_delta_pT[i][j]<<" second "<<sigma_positive_delta_pT[i][j]<<std::endl;
      }
    }


    std::cout<<"will save 1st deltapT postiive"<<std::endl;

/*
    gStyle->SetOptFit();
    std::cout<<"declaring canvas 1 fitting deltapTfunction"<<std::endl;
    TCanvas *canvas_1;   
    TCanvas *canvas_2;
    int j;
    std::cout<<"inside class reading fitting delta pTfunction"<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        canvas_1 = new TCanvas("canvas_1","positive muon deltapT",1000,1000);
        TF1 * f1 = new TF1("f1","gaus",-0.1,0.1);
        canvas_1->cd(); 
        histogram_delta_pT_positive_bin_reco[i][j]->Fit(f1,"R");
        gStyle->SetOptFit(00111111);
        mean_positive_delta_pT[i][j] = f1->GetParameter(1);
        std::cout<<"mean positive "<<mean_positive_delta_pT[i][j]<<std::endl;
        mean_positive_error_delta_pT[i][j] = f1->GetParError(1);
        sigma_positive_delta_pT[i][j] = f1->GetParameter(2);
        sigma_positive_error_delta_pT[i][j] = f1->GetParError(2);
       canvas_1->SaveAs(saving_path+"delta_pT_plots/delta_pT_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j]+".pdf"); 
       canvas_1->Close();
      }
    }
    std::cout<<"will save 1st deltapT postiive"<<std::endl;
   canvas_1->Delete();

    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
           canvas_2 = new TCanvas("canvas_2","negative muon deltapT",1000,1000);

        TF1 * f2 = new TF1("f2","gaus",-0.1,0.1);
        canvas_2->cd();
        histogram_delta_pT_negative_bin_reco[i][j]->Fit(f2,"R");
        gStyle->SetOptFit(00111111);
        mean_negative_delta_pT[i][j] = f2->GetParameter(1);
        std::cout<<"mean negative "<<mean_negative_delta_pT[i][j]<<std::endl;
        mean_negative_error_delta_pT[i][j] = f2->GetParError(1);
        sigma_negative_delta_pT[i][j] = f2->GetParameter(2);
        sigma_negative_error_delta_pT[i][j] = f2->GetParError(2);
       canvas_2->SaveAs(saving_path+"delta_pT_plots/delta_pT_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j]+".pdf"); 
       canvas_2->Close();


        std::cout<<" mean values "<<mean_negative_delta_pT[i][j]<<" error "<<mean_negative_error_delta_pT[i][j]<<std::endl;
     }
    }
    std::cout<<"will save 2nd deltapT postiive"<<std::endl;
*/
    }

 void class_reading::temperature_plot_delta_pT(){

   TH2F *temp_mean_positive = new TH2F("temp_mean_positive", "Mean (#frac{#Delta pT}{pT}) : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_mean_negative = new TH2F("temp_mean_negative", "Mean (#frac{#Delta pT}{pT}) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_positive = new TH2F("temp_sigma_positive", "#sigma (#frac{#Delta pT}{pT}) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_negative = new TH2F("temp_sigma_negative", "#sigma (#frac{#Delta pT}{pT}) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_positive = new TH2F("chi2_fit_positive", "#chi^2/ndf (#frac{#Delta pT}{pT}) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_negative = new TH2F("chi2_fit_negative", "#chi^2/ndf (#frac{#Delta pT}{pT}) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);

   for(int i=0; i< pT_bins; i++){ 
    for(int j=0; j<eta_bins;j++){
      temp_mean_positive->SetBinContent(i+1,j+1, mean_positive_delta_pT[i][j]);
      temp_sigma_positive->SetBinContent(i+1,j+1, sigma_positive_delta_pT[i][j]);
      temp_mean_negative->SetBinContent(i+1,j+1, mean_negative_delta_pT[i][j]);
      temp_sigma_negative->SetBinContent(i+1,j+1, sigma_negative_delta_pT[i][j]);
      chi2_fit_positive->SetBinContent(i+1, j+1, chi2_ndf_positive_delta_pT[i][j]);
      chi2_fit_negative->SetBinContent(i+1, j+1, chi2_ndf_negative_delta_pT[i][j]);
    }
   }
  // plotting this temperature plots 
  
    temp_mean_positive->GetZaxis()->SetRangeUser(-0.008,0.002);
    temp_sigma_positive->GetZaxis()->SetRangeUser(0,0.06);
    temp_mean_negative->GetZaxis()->SetRangeUser(-0.008,0.002);
    temp_sigma_negative->GetZaxis()->SetRangeUser(0,0.06);

		class_reading::plotting_hist_2D(temp_mean_positive," Mean #frac{#delta p_T}{p_T} : #mu^{+} ","temperature_mean_delta_pT_over_pT_2D_positive_distribution", " ");
		class_reading::plotting_hist_2D(temp_mean_negative," Mean #frac{#delta p_T}{p_T} : #mu^{-} ","temperature_mean_delta_pT_over_pT_2D_negative_distribution", " ");
		class_reading::plotting_hist_2D(temp_sigma_positive," #sigma #frac{#delta p_T}{p_T} : #mu^{+} ","temperature_sigma_delta_pT_over_pT_2D_positive_distribution", " ");
		class_reading::plotting_hist_2D(temp_sigma_negative," #sigma #frac{#delta p_T}{p_T} : #mu^{-} ","temperature_sigma_delta_pT_over_pT_2D_negative_distribution", " ");
		class_reading::plotting_hist_2D(chi2_fit_positive," #chi^2/ndf (#frac{#delta p_T}{p_T}) : #mu^{+} ","chi2_fit_positive", " ");
		class_reading::plotting_hist_2D(chi2_fit_negative," #chi^2/ndf (#frac{#delta p_T}{p_T}) : #mu^{-} ","chi2_fit_negative", " ");

    // done with temperature plots
 }
 
  // plotting mean and sigma from deltapT/pT in a final plot with mean and mean error
  void class_reading::plotting_mean_delta_pT(){

    TDirectory *dir_graph = (TDirectory*) output_root_file->mkdir("mean_resolution");
    dir_graph->cd();
    TGraphErrors *gr_mean1; 
    TGraphErrors *gr_mean2; 
    TGraphErrors *gr1; 
    TGraphErrors *gr2; 

    TCanvas *c1_mean_graph;
    TCanvas *c2_mean_graph;
    TCanvas *c1_graph;
    TCanvas *c2_graph;
 
    TString string_title;
    std::cout<<"strating function plotting mean delta pT"<<std::endl; 
    // mean and sigma plots for different eta bins
    ofstream new_text;
    new_text.open("final_file.txt");
   
    TMultiGraph *mg_positive =  new TMultiGraph();  
    TMultiGraph *mg_negative =  new TMultiGraph();  
    TMultiGraph *mg_sigma_positive =  new TMultiGraph();  
    TMultiGraph *mg_sigma_negative =  new TMultiGraph();  

     auto  * legend1 = new TLegend(0.15,0.15,0.30,0.35);
     auto  * legend2 = new TLegend(0.15,0.15,0.30,0.35);
     auto  * legend3 = new TLegend(0.11,0.50,0.23,0.70);
     auto  * legend4 = new TLegend(0.11,0.60,0.23,0.80);
     TString title_legend; 
    
     int colour[14] = {2,3, 4,6,7,8,9,11,12,13,14,15,16,17}; 
     for(int j=0 ; j< eta_bins; j++){
      float new_pt_measured_positive[pT_bins] = {0}; 
      float new_pt_error_positive[pT_bins] = {0}; 
      float new_pt_measured_negative[pT_bins] = {0}; 
      float new_pt_error_negative[pT_bins] = {0}; 

      float new_mean_positive_delta_pT[pT_bins] = {0}; 
      float new_mean_negative_delta_pT[pT_bins] = {0}; 
      float new_sigma_positive_delta_pT[pT_bins] = {0}; 
      float new_sigma_negative_delta_pT[pT_bins] = {0}; 

      float new_mean_positive_error_delta_pT[pT_bins] = {0}; 
      float new_mean_negative_error_delta_pT[pT_bins] = {0}; 
      float new_sigma_positive_error_delta_pT[pT_bins] = {0}; 
      float new_sigma_negative_error_delta_pT[pT_bins] = {0}; 

      for(int i=0; i<pT_bins; i++){
        std::cout<<"first bin of eta : pT mean value  "<<pt_measured_positive[i][j]<<std::endl; 
        new_pt_measured_positive[i]= pt_measured_positive[i][j];
        new_pt_error_positive[i]= pt_error_positive[i][j];
        new_pt_measured_negative[i]= pt_measured_negative[i][j];
        new_pt_error_negative[i]= pt_error_negative[i][j];

        std::cout<<"first bin of eta : deltapT mean "<<mean_positive_delta_pT[i][j]<<std::endl; 
        new_mean_positive_delta_pT[i]= mean_positive_delta_pT[i][j];
        new_mean_negative_delta_pT[i]= mean_negative_delta_pT[i][j];
        new_mean_positive_error_delta_pT[i]= mean_positive_error_delta_pT[i][j];
        new_mean_negative_error_delta_pT[i]= mean_negative_error_delta_pT[i][j];

        new_sigma_positive_delta_pT[i]= sigma_positive_delta_pT[i][j];
        new_sigma_negative_delta_pT[i]= sigma_negative_delta_pT[i][j];
        new_sigma_positive_error_delta_pT[i]= sigma_positive_error_delta_pT[i][j];
        new_sigma_negative_error_delta_pT[i]= sigma_negative_error_delta_pT[i][j];

        new_text<<"positive muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_positive_delta_pT[i]<<" pt "<<new_pt_measured_positive[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_positive_delta_pT[i]<<" pt error "<<new_pt_error_positive[i]<<std::endl;

        new_text<<"negative muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_negative_delta_pT[i]<<" pt "<<new_pt_measured_negative[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_negative_delta_pT[i]<<" pt error "<<new_pt_error_negative[i]<<std::endl;
        

      }
    
       std::cout<<"started graph plotting "<<std::endl; 
       gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_mean_positive_delta_pT,new_pt_error_positive, new_mean_positive_error_delta_pT);

    string_title = TString::Format("Mean #frac{#Delta pT}{pT} (#mu^{+}) : %.1f  #leq #eta < %.1f", eta_list[j], eta_list[j+1]);

    c1_mean_graph = new TCanvas("c1_mean_graph","Mean #frac{#Delta pT}{pT} mu^{+}",1000,1000);
    c1_mean_graph->cd();
    c1_mean_graph->SetLeftMargin(0.15);
    c1_mean_graph->SetTopMargin(0.2);
    gr_mean1->SetTitle(string_title);
    gr_mean1->GetXaxis()->SetTitleOffset(1.2);
    gr_mean1->GetYaxis()->SetTitleOffset(1.2);
    gr_mean1->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean1->GetYaxis()->SetRangeUser(-0.008,0.002);
    gr_mean1->Draw("AP*");
    gr_mean1->SetName("mean_delta_pT_positive_"+bin_number_eta[j]);
    gr_mean1->GetXaxis()->SetRangeUser(0,200);
    c1_mean_graph->Update();
    c1_mean_graph->SaveAs(saving_path+"graph_resolution_pT_mean_positive_eta_"+bin_number_eta[j]+".pdf");
    c1_mean_graph->Close();

    std::cout<<"ended 1st  graph plotting "<<std::endl; 
    c2_mean_graph = new TCanvas("c2_mean_graph","Mean #frac{#Delta pT}{pT} mu^{-}",1000,1000);
    c2_mean_graph->cd();
    c2_mean_graph->SetLeftMargin(0.15);
    c2_mean_graph->SetTopMargin(0.2);
    gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_mean_negative_delta_pT, new_pt_error_negative, new_mean_negative_error_delta_pT);

    string_title = TString::Format("Mean #frac{#Delta pT}{pT} (#mu^{-}) : %.1f  #leq #eta < %.1f ", eta_list[j], eta_list[j+1]);
    gr_mean2->SetTitle(string_title);
    gr_mean2->GetXaxis()->SetTitleOffset(1.2);
    gr_mean2->GetYaxis()->SetTitleOffset(1.2);
    gr_mean2->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean2->GetYaxis()->SetRangeUser(-0.008,0.002);
    gr_mean2->GetXaxis()->SetRangeUser(0,200);
    gr_mean2->SetName("mean_delta_pT_negative_"+bin_number_eta[j]);
    gr_mean2->Draw("AP*");
    c2_mean_graph->Update();
    c2_mean_graph->SaveAs(saving_path+"graph_resolution_pT_mean_negative_eta_"+bin_number_eta[j]+".pdf");
    c2_mean_graph->Close();

    gr_mean1->SetMarkerColor(colour[j]);
    gr_mean2->SetMarkerColor(colour[j]);
    mg_positive->Add(gr_mean1);
    mg_negative->Add(gr_mean2);

    //title_legend = TString::Format("%0.0f #leq #eta < %0.0f", eta_list_symbol[j],eta_list_symbol[j+1]);
    title_legend =   eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
    legend1->AddEntry(gr_mean1,title_legend,"P");
    legend2->AddEntry(gr_mean2,title_legend,"P");
    legend3->AddEntry(gr1,title_legend,"P");
    legend4->AddEntry(gr2,title_legend,"P");

    c1_graph = new TCanvas("c1_graph","#sigma #frac{#Delta pT}{pT} mu^{+}",1000,1000);
    c1_graph->cd();
    c1_graph->SetLeftMargin(0.15);
    c1_graph->SetTopMargin(0.2);
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_sigma_positive_delta_pT,new_pt_error_positive, new_sigma_positive_error_delta_pT);

    string_title = TString::Format("#sigma #frac{#Delta pT}{pT} (#mu^{+}) : %.1f  #leq #eta < %.1f", eta_list[j], eta_list[j+1]);
    gr1->SetTitle(string_title);
    gr1->GetXaxis()->SetTitleOffset(1.2);
    gr1->GetYaxis()->SetTitleOffset(1.2);
    gr1->GetXaxis()->SetTitle("pT(GEN)");
    gr1->GetYaxis()->SetRangeUser(0,0.06);
    gr1->GetXaxis()->SetRangeUser(0,200);
    gr1->SetName("sigma_delta_pT_positive_"+bin_number_eta[j]);
    gr1->Draw("AP*");
    c1_graph->Update();
    c1_graph->SaveAs(saving_path+"graph_resolution_pT_sigma_positive_eta_"+bin_number_eta[j]+".pdf");
    c1_graph->Close();

    c2_graph = new TCanvas("c2_graph","#sigma #frac{#Delta pT}{pT} mu^{-}",1000,1000);
    c2_graph->cd();
    c2_graph->SetLeftMargin(0.15);
    c2_graph->SetTopMargin(0.2);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_sigma_negative_delta_pT, new_pt_error_negative,  new_sigma_negative_error_delta_pT);
    string_title = TString::Format("#sigma #frac{#Delta pT}{pT} (#mu^{-}) : %.1f  #leq #eta < %.1f", eta_list[j], eta_list[j+1]);
    gr2->SetTitle(string_title);
    gr2->GetXaxis()->SetTitleOffset(1.2);
    gr2->GetYaxis()->SetTitleOffset(1.2);
    gr2->GetXaxis()->SetTitle("pT(GEN)");
    gr2->GetYaxis()->SetRangeUser(0,0.06);
    gr2->GetXaxis()->SetRangeUser(0,200);
    gr2->Draw("AP*");
    gr2->SetName("sigma_delta_pT_negative_"+bin_number_eta[j]);
    c2_graph->Update();
    c2_graph->SaveAs(saving_path+"graph_resolution_pT_sigma_negative_eta_"+bin_number_eta[j]+".pdf");
    c2_graph->Close(); 

    gr1->SetMarkerColor(colour[j]);
    gr2->SetMarkerColor(colour[j]);
   
    gr1->Write();
    gr2->Write();
    gr_mean1->Write();
    gr_mean2->Write();
    
    mg_sigma_positive->Add(gr1);
    mg_sigma_negative->Add(gr2);

	} // end of eta bins and plots
 
   TCanvas *mg_first_positive = new TCanvas("mg_first_positive", "Mean #frac{#Delta pT}{pT} mu^{+}",1000,1000); 
   TCanvas *mg_first_negative = new TCanvas("mg_first_negative", "Mean #frac{#Delta pT}{pT} mu^{+}",1000,1000); 
   TCanvas *mg_second_positive = new TCanvas("mg_second_positive", "#sigma #frac{#Delta pT}{pT} mu^{+}",1000,1000); 
   TCanvas *mg_second_negative = new TCanvas("mg_second_negative", "#sigma #frac{#Delta pT}{pT} mu^{-}",1000,1000); 

   mg_first_positive->cd();
   mg_positive->Draw("AP");
   mg_positive->SetTitle("Mean (#frac{#Delta pT}{pT}) : #mu^{+} ");
   mg_positive->GetXaxis()->SetTitle("pT (GeV) ");
   mg_positive->GetXaxis()->SetRangeUser(0,200);
   mg_positive->GetYaxis()->SetRangeUser(-0.008,0.002);
   legend1->Draw();
   mg_first_positive->SaveAs(saving_path+"mean_delta_pT_positive.pdf");

   mg_first_negative->cd();
   mg_negative->Draw("AP");
   mg_negative->SetTitle("Mean (#frac{#Delta pT}{pT}) : #mu^{-} ");
   mg_negative->GetXaxis()->SetTitle("pT (GeV) ");
   mg_negative->GetXaxis()->SetRangeUser(0,200);
   mg_negative->GetYaxis()->SetRangeUser(-0.008,0.002);
   legend2->Draw();
   mg_first_negative->SaveAs(saving_path+"mean_delta_pT_negative.pdf");

   mg_second_positive->cd();
   mg_sigma_positive->Draw("AP");
   mg_sigma_positive->SetTitle("#sigma (#frac{#Delta pT}{pT}) : #mu^{+} ");
   mg_sigma_positive->GetXaxis()->SetTitle("pT (GeV) ");
   mg_sigma_positive->GetXaxis()->SetRangeUser(0,200);
   mg_sigma_positive->GetYaxis()->SetRangeUser(0,0.06);
   //legend3->Draw();
   mg_second_positive->SaveAs(saving_path+"sigma_delta_pT_positive.pdf");

   mg_second_negative->cd();
   mg_sigma_negative->Draw("AP");
   mg_sigma_negative->SetTitle("#sigma (#frac{#Delta pT}{pT}) : #mu^{-} ");
   mg_sigma_negative->GetXaxis()->SetTitle("pT (GeV) ");
   mg_sigma_negative->GetXaxis()->SetRangeUser(0,200);
   mg_sigma_negative->GetYaxis()->SetRangeUser(0.0,0.06);
   //legend4->Draw();
   mg_second_negative->SaveAs(saving_path+"sigma_delta_pT_negative.pdf");

  }

    // implementing DSCB fit
    void class_reading :: fitting_DSCB_positive(TH1F * hist_fit, TString saving_name, TString title_name, int pt_bin_number, int eta_bin_number){
      TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
      TF1* f1 = new TF1("f1", "gaus", -0.02,0.02);
      hist_fit_clone->Fit(f1,"R+");

      float mean_gaus =0; float sigma_gaus = 0; 
      mean_gaus = f1->GetParameter(1);
      sigma_gaus = f1->GetParameter(2);
 
      RooRealVar delta_pT("delta_pT", "deltapT", -0.2,0.2, "");
      RooDataHist histo("histo","dataset with var",delta_pT,hist_fit);
      RooRealVar Mean("Mean", "Mean",mean_gaus, -0.01, 0.01);
      RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.005,0.1);//sigma[decay]);
      RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
      RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 200);//expL[decay]);
      RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.1, 50);//alphaR[decay]);
      RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 200);//expR[decay]);
      RooMyPDF_DSCB DSCB("DSCB", "DSCB", delta_pT, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);


/*     RooRealVar delta_pT("delta_pT", "deltapT", -0.2,0.2, "");
     RooDataHist histo("histo","dataset with var",delta_pT,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, -0.02, 0.02);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.0001,0.1);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.0001, 500);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 500);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.0001, 500);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 500);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", delta_pT, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
*/ 
      output_file<<"fitting positive muon : mean "<<mean_gaus<<" sigma : "<<sigma_gaus<<std::endl;
     // Fit the mass into DSCB

     TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
     c_MC->SetFrameFillColor(0);
 
     RooPlot* xframe = delta_pT.frame(RooFit::Title(title_name));
     histo.plotOn(xframe);
 
     Int_t color = kRed+2;
     Double_t size_text = 0.020;
     DSCB.fitTo(histo, Range(-0.15,0.15));
     DSCB.plotOn(xframe, RooFit::LineColor(color),Name("#frac{#delta p_T}{p_T}"));
     DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
     c_MC->cd();
     int pt_number = pt_bin_number+1;
     int eta_number = eta_bin_number+1;
     TString title_s =TString::Format("delta_pT_fit_positive_mu_pt_%d_eta_%d",pt_number, eta_number); 
     c_MC->SetName(title_s);
     c_MC->SetTopMargin(1.5);
     xframe->getAttText()->SetTextSize(size_text);
     xframe->getAttText()->SetTextColor(color);
//     xframe->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
//     xframe->GetXaxis()->CenterTitle(true);
     xframe->GetYaxis()->SetTitle("");
     xframe->GetXaxis()->SetTitle("");
//     xframe->GetXaxis()->SetTitleOffset(1.4);
     xframe->chiSquare();
     //xframe->SetStats();
     xframe->Draw();
    
     
     //const TH1* histogram = histo.createHistogram("delta_pT");

    TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
      leg2->SetFillColor(kWhite);
      leg2->SetLineColor(kBlack);
      leg2->AddEntry("histo","#frac{#Delta pT}{pT}", "EP");
      leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
      leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
//      leg2->AddEntry("histogram->GetBinContent(0)",Form("Underflow= %.0f",histogram->GetBinContent(0)),"");
//      leg2->AddEntry("histogram->GetBinContent(histogram->GetNbinsX() + 1)",Form("Overflow= %.0f",histogram->GetBinContent(histogram->GetNbinsX() + 1)),"");
      leg2->Draw("same");
 
//    std::cout<<"Underflow= "<<histogram->GetBinContent(0)<<std::endl;
//    std::cout<<"Overflow= "<<histogram->GetBinContent(histogram->GetNbinsX()+1)<<std::endl;
     gStyle->SetOptStat();
     c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
     c_MC->Write(); 
     c_MC->Close();
     std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;

     
     mean_positive_delta_pT[pt_bin_number][eta_bin_number] = Mean.getVal();
     mean_positive_error_delta_pT[pt_bin_number][eta_bin_number] = Mean.getError();

     sigma_positive_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getVal();
     sigma_positive_error_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getError();

     chi2_ndf_positive_delta_pT[pt_bin_number][eta_bin_number] =  xframe->chiSquare();
 }
      // for negative muons
    void class_reading :: fitting_DSCB_negative(TH1F * hist_fit, TString saving_name, TString title_name, int pt_bin_number, int eta_bin_number){
      TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
      TF1* f1 = new TF1("f1", "gaus", -0.02,0.02);
      hist_fit_clone->Fit(f1,"R+");

      float mean_gaus =0; float sigma_gaus = 0; 
      mean_gaus = f1->GetParameter(1);
      sigma_gaus = f1->GetParameter(2);

      output_file<<"fitting negative muon : mean "<<mean_gaus<<" sigma : "<<sigma_gaus<<std::endl;
      RooRealVar delta_pT("delta_pT", "deltapT", -0.2,0.2, "");
      RooDataHist histo("histo","dataset with var",delta_pT,hist_fit);
      RooRealVar Mean("Mean", "Mean",mean_gaus, -0.01, 0.01);
      RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.005,0.1);//sigma[decay]);
      RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
      RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 200);//expL[decay]);
      RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.1, 50);//alphaR[decay]);
      RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 200);//expR[decay]);
      RooMyPDF_DSCB DSCB("DSCB", "DSCB", delta_pT, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);



      // Fit the mass into DSCB
/*     RooRealVar delta_pT("delta_pT", "deltapT", -0.2,0.2, "");
     RooDataHist histo("histo","dataset with var",delta_pT,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, -0.02, 0.02);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.0001,0.1);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.0001, 500);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 500);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.0001, 500);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 500);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", delta_pT, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
*/ 
     TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
     c_MC->SetFrameFillColor(0);
 
     RooPlot* xframe = delta_pT.frame(RooFit::Title(title_name));
     histo.plotOn(xframe);
 
     Int_t color = kRed+2;
     Double_t size_text = 0.020;
     DSCB.fitTo(histo, Range(-0.15,0.15));
     DSCB.plotOn(xframe, RooFit::LineColor(color),Name("#frac{#delta p_T}{p_T}"));
     DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
     c_MC->cd();
     c_MC->SetBottomMargin(1.5);
     int pt_number = pt_bin_number+1;
     int eta_number = eta_bin_number+1;
     TString title_s =TString::Format("delta_pT_fit_negative_mu_pt_%d_eta_%d",pt_number, eta_number); 
     c_MC->SetName(title_s);
     xframe->getAttText()->SetTextSize(size_text);
     xframe->getAttText()->SetTextColor(color);
 //    xframe->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
 //    xframe->GetXaxis()->CenterTitle(true);
     xframe->GetYaxis()->SetTitle("");
     xframe->GetXaxis()->SetTitle("");
  //   xframe->GetXaxis()->SetTitleOffset(1.4);
     xframe->chiSquare();
     //xframe->SetStats();
     xframe->Draw();
    
     
     //const TH1* histogram = histo.createHistogram("delta_pT");

    TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
      leg2->SetFillColor(kWhite);
      leg2->SetLineColor(kBlack);
      leg2->AddEntry("histo","#frac{#Delta pT}{pT}", "EP");
      leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
      leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
//      leg2->AddEntry("histogram->GetBinContent(0)",Form("Underflow= %.0f",histogram->GetBinContent(0)),"");
//      leg2->AddEntry("histogram->GetBinContent(histogram->GetNbinsX() + 1)",Form("Overflow= %.0f",histogram->GetBinContent(histogram->GetNbinsX() + 1)),"");
      leg2->Draw("same");
 
//    std::cout<<"Underflow= "<<histogram->GetBinContent(0)<<std::endl;
//    std::cout<<"Overflow= "<<histogram->GetBinContent(histogram->GetNbinsX()+1)<<std::endl;
     gStyle->SetOptStat();
     c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
     c_MC->Write(); 
     c_MC->Close();
     std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;

     mean_negative_delta_pT[pt_bin_number][eta_bin_number] = Mean.getVal();
     mean_negative_error_delta_pT[pt_bin_number][eta_bin_number] = Mean.getError();

     sigma_negative_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getVal();
     sigma_negative_error_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getError();
     
     chi2_ndf_negative_delta_pT[pt_bin_number][eta_bin_number] =  xframe->chiSquare();
 }


 void delta_pT_analysis(){ 
 
   // this code is to analyze all the histograms produced for the Jpsi sample : histograms ar     e Z mass histogram in all the pT categories, for both mu^+ and mu^-   
   // derived class reading is the class where I am performing operation to fit the histogram     s , and dervied class reading is using functions from base class reading to perform this oper     ations   
   gROOT->SetBatch(kTRUE);    
   std::cout<<" initiated a class obj"<<std::endl;   
   derived_class_reading obj; 
 
   // initializing is to open the root file and to read all the histograms 
   obj.initializing();  
 
   // plotting every histogram on canvas with the plot_histograms function   
//   obj.plot_histograms();  
  
   // plotting inclusive and binned distribution on same canvas  
   //obj.plotting_inclusive_bin(); 
 
   // to do BW and DSCB fits and find mean and sigma values  and this mean and sigma values a     re further stored and used to make final plots in the functio nmean_sigma_calculation 
//   obj.fitting_histograms(); 
 
   // saving text file also has the calculation for mean and sigma which will be used further     e 
//   obj.saving_text_file(); 
   
   // graph to plot mean z mass with pT 
//   obj.mean_sigma_calculation(); 
    
   // to save the pT distributions for individual histograms 
   // calling a derived class functions available with these functions 
 
   std::cout<<" class obj declared successfully"<<std::endl;   
   obj.saving_histograms(); 
   std::cout<<" saving histogram delta pT"<<std::endl;   
// 
   obj.fitting_delta_pT(); 
  } 




