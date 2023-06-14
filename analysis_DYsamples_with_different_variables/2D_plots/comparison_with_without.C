void comparison_with_without(){
TString input_path = "/afs/cern.ch/user/n/nrawal/work/UL_sample_generation_with_without_rochester/analysis_DYsamples_with_different_variables/2D_plots/";
TFile *file1 = TFile::Open(input_path+"without_rochester_histograms_delta_pT.root");
TFile *file2 = TFile::Open(input_path+"with_rochester_histograms_delta_pT.root");

TDirectoryFile *dir1 = (TDirectoryFile*) file1->Get("mean_resolution");
TDirectoryFile *dir2 = (TDirectoryFile*) file2->Get("mean_resolution");
TString list[14] = {"1st","2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th", "11th", "12th", "13th", "14th"};
TString eta_list[15] ={"-2.4", "-2.1", "-1.9", "-1.6", "-1.2",  "-0.8", "-0.4", "0","0.4", "0.8", "1.2", "1.6", "1.9", "2.1","2.4"};
//TString eta_list[5] = {"0.0", "0.9", "1.2", "2.0","2.4"};
TString title_mean_positive;
TString title_mean_negative;
TString title_sigma_positive;
TString title_sigma_negative;
TGraphErrors *gr1_mean_positive;
TGraphErrors *gr1_sigma_positive;
TGraphErrors *gr2_mean_positive;
TGraphErrors *gr2_sigma_positive;
TGraphErrors *gr1_mean_negative;
TGraphErrors *gr1_sigma_negative;
TGraphErrors *gr2_mean_negative;
TGraphErrors *gr2_sigma_negative;

TMultiGraph *mg_mean_positive;
TMultiGraph *mg_sigma_positive;
TMultiGraph *mg_mean_negative;
TMultiGraph *mg_sigma_negative;


TCanvas *c_positive ;
TCanvas *c_negative ;
TCanvas *c_positive_2 ;
TCanvas *c_negative_2 ;
TString title_s;
for(int i =0 ; i<14; i++){
 title_mean_positive = "mean_delta_pT_positive_"+list[i];
 title_sigma_positive = "sigma_delta_pT_positive_"+list[i];
 gr1_mean_positive = (TGraphErrors*) dir1->Get(title_mean_positive);
 gr1_sigma_positive = (TGraphErrors*) dir1->Get(title_sigma_positive);
 gr1_mean_positive->SetMarkerColor(2);
 gr1_sigma_positive->SetMarkerColor(2);
/*
  c_negative_2 = new TCanvas("c_negative","canvas",1000,1000);
  c_negative_2->cd();
  gr1_mean_positive->Draw("AP");
  c_negative_2->SaveAs("tmp_delta_pT_negative_bin_"+list[i]+".pdf");
  c_negative_2->Close();
*/


  title_mean_positive = "mean_delta_pT_positive_"+list[i];
  title_sigma_positive = "sigma_delta_pT_positive_"+list[i];
  TGraphErrors *gr2_mean_positive = (TGraphErrors*) dir2->Get(title_mean_positive);
  TGraphErrors *gr2_sigma_positive = (TGraphErrors*) dir2->Get(title_sigma_positive);
 gr2_mean_positive->SetMarkerColor(4);
 gr2_sigma_positive->SetMarkerColor(4);

 title_mean_negative = "mean_delta_pT_negative_"+list[i];
 title_sigma_negative = "sigma_delta_pT_negative_"+list[i];
 TGraphErrors *gr1_mean_negative = (TGraphErrors*) dir1->Get(title_mean_negative);
 TGraphErrors *gr1_sigma_negative = (TGraphErrors*) dir1->Get(title_sigma_negative);
 gr1_mean_negative->SetMarkerColor(2);
 gr1_sigma_negative->SetMarkerColor(2);

  title_mean_negative = "mean_delta_pT_negative_"+list[i];
  title_sigma_negative = "sigma_delta_pT_negative_"+list[i];
  gr2_mean_negative = (TGraphErrors*) dir2->Get(title_mean_negative);
  gr2_sigma_negative = (TGraphErrors*) dir2->Get(title_sigma_negative);
 gr2_mean_negative->SetMarkerColor(4);
 gr2_sigma_negative->SetMarkerColor(4);

  mg_mean_positive = new TMultiGraph();
  mg_sigma_positive = new TMultiGraph();
  mg_mean_negative = new TMultiGraph();
  mg_sigma_negative = new TMultiGraph();

  auto legend1 = new TLegend(0.70,0.70,0.89,0.89);
  auto legend2 = new TLegend(0.70,0.70,0.89,0.89);
  auto legend3 = new TLegend(0.70,0.70,0.89,0.89);
  auto legend4 = new TLegend(0.70,0.70,0.89,0.89);

  legend1->AddEntry(gr1_mean_positive," without rochester correction ","P");
  legend1->AddEntry(gr2_mean_positive," with rochester correction ","P");

  legend2->AddEntry(gr1_mean_negative," without rochester correction ","P");
  legend2->AddEntry(gr2_mean_negative," with rochester correction ","P");

  legend3->AddEntry(gr1_sigma_positive," without rochester correction ","P");
  legend3->AddEntry(gr2_sigma_positive," with rochester correction ","P");

  legend4->AddEntry(gr1_sigma_negative," without rochester correction ","P");
  legend4->AddEntry(gr2_sigma_negative," with rochester correction ","P");

  mg_mean_positive->Add(gr1_mean_positive);
  mg_mean_positive->Add(gr2_mean_positive);

  mg_sigma_positive->Add(gr1_sigma_positive);
  mg_sigma_positive->Add(gr2_sigma_positive);

  mg_mean_negative->Add(gr1_mean_negative);
  mg_mean_negative->Add(gr2_mean_negative);

  mg_sigma_negative->Add(gr1_sigma_negative);
  mg_sigma_negative->Add(gr2_sigma_negative);
 
  title_s = "Mean #frac{#Delta pT}{pT} (#mu^{+}) :"+eta_list[i]+" #leq #eta < "+eta_list[i+1];
  mg_mean_positive->SetTitle(title_s);
  title_s = "Mean #frac{#Delta pT}{pT} (#mu^{-}) :"+eta_list[i]+" #leq #eta < "+eta_list[i+1];
  mg_mean_negative->SetTitle(title_s);
  title_s = "#sigma #frac{#Delta pT}{pT} (#mu^{+}) :"+eta_list[i]+" #leq #eta < "+eta_list[i+1];
  mg_sigma_positive->SetTitle(title_s);
  title_s = "#sigma #frac{#Delta pT}{pT} (#mu^{-}) :"+eta_list[i]+" #leq #eta < "+eta_list[i+1];
  mg_sigma_negative->SetTitle(title_s);

  // plot them on same canvas
  c_positive = new TCanvas("c_positive","canvas",1000,1000);
  c_positive->cd();
  mg_mean_positive->Draw("AP");
  mg_mean_positive->GetYaxis()->SetRangeUser(-0.007,0.002);
  legend1->Draw();
  c_positive->SaveAs("mean_delta_pT_positive_bin_"+list[i]+".pdf");
  c_positive->Close();

  c_positive_2 = new TCanvas("c_positive","canvas",1000,1000);
  c_positive_2->cd();
  mg_sigma_positive->Draw("AP");
  mg_sigma_positive->GetYaxis()->SetRangeUser(0,0.06);
  legend3->Draw();
  c_positive_2->SaveAs("sigma_delta_pT_positive_bin_"+list[i]+".pdf");
  c_positive_2->Close();

  c_negative = new TCanvas("c_negative","canvas",1000,1000);
  c_negative->cd();
  mg_mean_negative->Draw("AP");
  mg_mean_negative->GetYaxis()->SetRangeUser(-0.007,0.002);
  legend2->Draw();
  c_negative->SaveAs("mean_delta_pT_negative_bin_"+list[i]+".pdf");
  c_negative->Close();

  c_negative_2 = new TCanvas("c_negative","canvas",1000,1000);
  c_negative_2->cd();
  mg_sigma_negative->Draw("AP");
  mg_sigma_negative->GetYaxis()->SetRangeUser(0,0.06);
  legend4->Draw();
  c_negative_2->SaveAs("sigma_delta_pT_negative_bin_"+list[i]+".pdf");
  c_negative_2->Close();
  
}


}
