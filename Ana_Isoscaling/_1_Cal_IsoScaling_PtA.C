#include "TMath.h"
using namespace TMath;

//this script will draw all the different particle together.
void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle, int MarkerStyle, int MarkerColor);
void Set_1D_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle);

TH1D* Select_Pt_A(string SystemName, string ParticleName, TH2D* h2_Pt_A_Rapidity,double* Range);
TF1* Extend_Pt_A(string SystemName, string ParticleName, double* RapidityRange, TH1D* h1_Pt_A, double* FitRange, double* ExtendRange);

TH1D* Get_Ratio_T_3He(string NameTem,TH1D* h1_PtA_Lab_T, TH1D* h1_PtA_Lab_3He);
TH1D* Get_PseudoNeutron(string NameTem,TH1D* h1_Ratio_T_3He, TH1D* h1_PtA_Lab_P);
TH1D* Get_CIp(string NameTem, int ParticleNum,TH1D** h1_PtA_Lab);
TH1D* Get_CIn(string NameTem, int ParticleNum,TH1D** h1_PtA_Lab,TH1D* h1_PseudoNeutron);
TH1D* GetTemperature(string NameTem,TH1D** h1_PtA_Lab);

void RemovePoint_LargeRelError(TH1D* h1,double RelErr);
void Cal_Alpha_Beta(int Index, TCanvas* c1_N, TCanvas* c1_Z, int ParticleNum, double* Z, double* N, double* R21, double* R21_Err,double* Alpha, double* Alpha_Err, double* Beta, double* Beta_Err, double* C, double* C_Err);

void _1_Cal_IsoScaling_PtA()
{
  gStyle->SetErrorX(0);

  string SystemTag_1 = "Ca40Ni58E140"; 
  string SystemTag_2 = "Ca48Ni64E140"; 
  
  double RelErr_Cut = 0.05;

  bool IsExtend = 0;
  int RebinNum = 20;
  bool Is_R21Err_Applied=1;
  
  int iEff = 1;
  int ib = 0;
  double RapidityRange[2] = {0.4,0.6};
  
  const int EffOpt = 3;
  string EffCor_Tag[] = {"_noEff", "_GeoEff",  "_GeoReactionEff"};
  
  const int ParticleNum = 6;
  string ParticleName[]={"p","d","t","3He","4He","6He","6Li","7Li"};
  int LineColor[] = {1,1,1,2,2,2};
//  int MarkerStyle[] = {1,1,1,1,1,1};
  int MarkerStyle[] = {21,25,36,22,26,33};
  
  //const int ImpactParaNum = 3;
  //string ImpactParaTag[] = {"b_0.0_3.0","b_3.0_5.0","b_5.0_7.0"};
  
  const int ImpactParaNum = 3;
  string ImpactParaTag[] = {"bHat_0.00_0.40","bHat_0.40_0.70","bHat_0.70_1.10"};
  string ImpactPara_Select_Tag = "bHat";//"ImpactPars"
  
  char NameTem[200];
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_1.c_str(),ImpactParaTag[ib].c_str());
  TFile* f1_MergedESpec_1 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_2.c_str(),ImpactParaTag[ib].c_str());
  TFile* f1_MergedESpec_2 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;

  TH2D* h2_y_PtA_Lab[2][ParticleNum];
  TH1D* h1_PtA_Lab[2][ParticleNum];
  TF1* f1_PtA_D[2];
  TF1* f1_PtA_T[2];
  TF1* f1_PtA_3He[2];
  
  TH1D* h1_R21[ParticleNum];
  TH1D* h1_Temperature[2];
  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    cout<<"Read : "<<NameTem<<endl;
    f1_MergedESpec_1->cd();
    h2_y_PtA_Lab[0][iP] = (TH2D*) f1_MergedESpec_1->Get(NameTem);
    sprintf(NameTem,"h2_%s_%s_%s_y_PtA_Lab%s_bhat",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
    h2_y_PtA_Lab[0][iP]->SetName(NameTem);
    
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    f1_MergedESpec_2->cd();
    h2_y_PtA_Lab[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem);
    sprintf(NameTem,"h2_%s_%s_%s_y_PtA_Lab%s_bhat",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
    h2_y_PtA_Lab[1][iP]->SetName(NameTem);
    
    //the below is for projecting Pt_A.
    f1_MergedESpec_1->cd();
    h1_PtA_Lab[0][iP] = Select_Pt_A(SystemTag_1, ParticleName[iP], h2_y_PtA_Lab[0][iP],RapidityRange);
    f1_MergedESpec_2->cd();
    h1_PtA_Lab[1][iP] = Select_Pt_A(SystemTag_2, ParticleName[iP], h2_y_PtA_Lab[1][iP],RapidityRange);
    
    h1_PtA_Lab[0][iP]->SetLineColor(LineColor[iP]);
    h1_PtA_Lab[1][iP]->SetLineColor(LineColor[iP]);
    h1_PtA_Lab[0][iP]->SetMarkerColor(LineColor[iP]);
    h1_PtA_Lab[1][iP]->SetMarkerColor(LineColor[iP]);
    h1_PtA_Lab[0][iP]->SetMarkerStyle(MarkerStyle[iP]);
    h1_PtA_Lab[1][iP]->SetMarkerStyle(MarkerStyle[iP]);
    
    //the below is for extending the spectrum.
    if(iP==1 && IsExtend==1) //Deuteron
    {
      //double FitRange[2] = {210,380}; double ExtendRange[2] = {370,600};
      double FitRange[2] = {210,380}; double ExtendRange[2] = {370,600};
      f1_MergedESpec_1->cd();
      f1_PtA_D[0] = Extend_Pt_A(SystemTag_1, ParticleName[iP], RapidityRange, h1_PtA_Lab[0][iP], FitRange, ExtendRange);
      f1_MergedESpec_2->cd();
      f1_PtA_D[1] = Extend_Pt_A(SystemTag_2, ParticleName[iP], RapidityRange, h1_PtA_Lab[1][iP], FitRange, ExtendRange);
    }
    else if(iP==2&& IsExtend==1) //Triton
    {
      //double FitRange[2] = {210,320}; double ExtendRange[2] = {310,600};
      double FitRange[2] = {210,320}; double ExtendRange[2] = {310,600};
      f1_MergedESpec_1->cd();
      f1_PtA_T[0] = Extend_Pt_A(SystemTag_1, ParticleName[iP], RapidityRange, h1_PtA_Lab[0][iP], FitRange, ExtendRange);
      f1_MergedESpec_2->cd();
      f1_PtA_T[1] = Extend_Pt_A(SystemTag_2, ParticleName[iP], RapidityRange, h1_PtA_Lab[1][iP], FitRange, ExtendRange);
    }
    else if(iP==3&& IsExtend==1) //3He
    {
      //double FitRange[2] = {210,500}; double ExtendRange[2] = {300,600};
      double FitRange[2] = {210,500}; double ExtendRange[2] = {300,600};
      f1_MergedESpec_1->cd();
      f1_PtA_3He[0] = Extend_Pt_A(SystemTag_1, ParticleName[iP], RapidityRange, h1_PtA_Lab[0][iP], FitRange, ExtendRange);
      f1_MergedESpec_2->cd();
      f1_PtA_3He[1] = Extend_Pt_A(SystemTag_2, ParticleName[iP], RapidityRange, h1_PtA_Lab[1][iP], FitRange, ExtendRange);
    }
    
    //the below is for rebin the data.
    h1_PtA_Lab[0][iP]->Rebin(RebinNum);
    h1_PtA_Lab[0][iP]->Scale(1.0/RebinNum);
    h1_PtA_Lab[1][iP]->Rebin(RebinNum);
    h1_PtA_Lab[1][iP]->Scale(1.0/RebinNum);
  }
  
  //Draw the Pt/A~Rapidity
  sprintf(NameTem,"c1_Pt_A_Rapidity_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Rapidity_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Sys1->Divide(3,2,0,0);
  sprintf(NameTem,"c1_Pt_A_Rapidity_%s_%s_%s_y_%.1f_%.1f",SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Rapidity_Sys2 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Sys2->Divide(3,2,0,0);
  for(int i=0;i<ParticleNum;i++)
  {
    c1_Pt_A_Rapidity_Sys1->cd(i+1);//->SetGrid(1,1);
    h2_y_PtA_Lab[0][i]->Draw("col");
    c1_Pt_A_Rapidity_Sys2->cd(i+1);//->SetGrid(1,1);
    h2_y_PtA_Lab[1][i]->Draw("col");
  }
  
  
  
  //Draw the Pt/A(MeV/c/u)
  sprintf(NameTem,"c1_Pt_A_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Sys1 = new TCanvas(NameTem,NameTem,1);
  sprintf(NameTem,"c1_Pt_A_%s_%s_%s_y_%.1f_%.1f",SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Sys2 = new TCanvas(NameTem,NameTem,1);
  TLegend* legend_PID = new TLegend(0.6,0.6,0.8,0.8);
  c1_Pt_A_Sys1->Divide(2,1);
  for(int i=0;i<ParticleNum-1;i++)
  {
    c1_Pt_A_Sys1->cd(1)->SetLogy(1);
    if(i==0) 
    { 
      h1_PtA_Lab[0][i]->Draw();
      h1_PtA_Lab[0][i]->GetXaxis()->SetRangeUser(130,400);
      h1_PtA_Lab[0][i]->GetYaxis()->SetRangeUser(2E-5,5E-2);
    }
    else { h1_PtA_Lab[0][i]->Draw("same"); }
    h1_PtA_Lab[0][i]->GetXaxis()->SetTitle("P_{t}/A(MeV/c)");
    h1_PtA_Lab[0][i]->GetYaxis()->SetTitle("d^{2}M/dP_{t}d_{y}(MeV^{-1}c^{-1})");
    c1_Pt_A_Sys2->cd()->SetGrid(1,1);
    if(i==0) 
    { 
      h1_PtA_Lab[1][i]->Draw(); 
      h1_PtA_Lab[1][i]->GetXaxis()->SetRangeUser(100,400);
    }
    else { h1_PtA_Lab[1][i]->Draw("same"); }
    
    legend_PID->AddEntry(h1_PtA_Lab[1][i],ParticleName[i].c_str(),"lp");
  }
  c1_Pt_A_Sys1->cd(1);  legend_PID->Draw();
  c1_Pt_A_Sys2->cd();  legend_PID->Draw();
  
  //Draw R21
  sprintf(NameTem,"c1_R21_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_R21 = new TCanvas(NameTem,NameTem,1);

  
  for(int i=0;i<ParticleNum;i++)
  {
    //c1_R21->cd()->SetGrid(1,1);
    c1_Pt_A_Sys1->cd(2);
    sprintf(NameTem,"h1_R21_%s_%s_%s_%s_%s_y_%.1f_%.1f",ParticleName[i].c_str(),SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
    h1_R21[i] = (TH1D*) h1_PtA_Lab[1][i]->Clone(NameTem);
    h1_R21[i]->Divide(h1_PtA_Lab[0][i]);
    h1_R21[i]->GetYaxis()->SetTitle("R21");
    h1_R21[i]->GetXaxis()->SetTitle("P_{t}/A(MeV/c)");
    
    RemovePoint_LargeRelError(h1_R21[i],RelErr_Cut);
    if(i==0)
    {
      h1_R21[i]->Draw();
      h1_R21[i]->GetXaxis()->SetRangeUser(100,400);
      h1_R21[i]->GetYaxis()->SetTitle("R21");
      h1_R21[i]->GetYaxis()->SetRangeUser(0.5,2);
    }
    else if(i!=5) { h1_R21[i]->Draw("same"); }
  }
  //legend_PID->Draw();
  
  //Draw temperature
  sprintf(NameTem,"c1_Temperature_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Temperature = new TCanvas(NameTem,NameTem,1);
  sprintf(NameTem,"h1_Temperature_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  h1_Temperature[0] = GetTemperature(NameTem,h1_PtA_Lab[0]);
  h1_Temperature[0]->Draw();
  h1_Temperature[0]->SetMarkerColor(2); h1_Temperature[0]->SetMarkerStyle(20);
  RemovePoint_LargeRelError(h1_Temperature[0],RelErr_Cut);
  
  sprintf(NameTem,"h1_Temperature_%s_%s_%s_y_%.1f_%.1f",SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  h1_Temperature[1] = GetTemperature(NameTem,h1_PtA_Lab[1]);
  h1_Temperature[1]->Draw("same");
  h1_Temperature[1]->SetMarkerColor(8); h1_Temperature[1]->SetMarkerStyle(25);
  RemovePoint_LargeRelError(h1_Temperature[1],RelErr_Cut);
  
  TLegend* legend_Temperature = new TLegend(0.6,0.6,0.8,0.8);
  legend_Temperature->AddEntry(h1_Temperature[0],SystemTag_1.c_str(),"lp");
  legend_Temperature->AddEntry(h1_Temperature[1],SystemTag_2.c_str(),"lp");
  legend_Temperature->Draw();
  
  //show different R21
  sprintf(NameTem,"c1_R21_NZ_NotEqul_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_NZ_NotEqul = new TCanvas(NameTem,NameTem,1);
  TLegend* legend_R21_NZ_NotEqul_Cluster = new TLegend(0.6,0.6,0.8,0.8);
  
  TH1D* h1_R21_T = (TH1D*) h1_PtA_Lab[1][2]->Clone("h1_R21_T"); h1_R21_T->Divide(h1_PtA_Lab[0][2]);
  h1_R21_T->GetXaxis()->SetRangeUser(100,400);
  h1_R21_T->GetYaxis()->SetRangeUser(0.5,2.0);
  h1_R21_T->GetYaxis()->SetTitle("R21");
  h1_R21_T->SetMarkerColor(6); h1_R21_T->SetLineColor(6); h1_R21_T->SetMarkerStyle(23);
  h1_R21_T->Draw();
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_T,"T","lp");
  RemovePoint_LargeRelError(h1_R21_T,RelErr_Cut);
  
  TH1D* h1_R21_Alpha_Divide_P = (TH1D*) h1_PtA_Lab[1][4]->Clone("h1_R21_Alpha_Divide_P");
  h1_R21_Alpha_Divide_P->Divide(h1_PtA_Lab[0][4]);
  h1_R21_Alpha_Divide_P->Multiply(h1_PtA_Lab[0][0]);
  h1_R21_Alpha_Divide_P->Divide(h1_PtA_Lab[1][0]);
  h1_R21_Alpha_Divide_P->SetMarkerColor(6); h1_R21_Alpha_Divide_P->SetLineColor(6); h1_R21_Alpha_Divide_P->SetMarkerStyle(26);
  h1_R21_Alpha_Divide_P->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_Alpha_Divide_P,"#alpha/P","lp");
  RemovePoint_LargeRelError(h1_R21_Alpha_Divide_P,RelErr_Cut);
  
  ////////////////////////////////////////////////////////////////////////////////
  
  TH1D* h1_R21_Alpha_Divide_3He = (TH1D*) h1_PtA_Lab[1][4]->Clone("h1_R21_Alpha_Divide_3He");
  h1_R21_Alpha_Divide_3He->Divide(h1_PtA_Lab[0][4]);
  h1_R21_Alpha_Divide_3He->Multiply(h1_PtA_Lab[0][3]);
  h1_R21_Alpha_Divide_3He->Divide(h1_PtA_Lab[1][3]);
  h1_R21_Alpha_Divide_3He->SetMarkerColor(1); h1_R21_Alpha_Divide_3He->SetLineColor(1); h1_R21_Alpha_Divide_3He->SetMarkerStyle(24);
  h1_R21_Alpha_Divide_3He->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_Alpha_Divide_3He,"#alpha/3He","lp");
  RemovePoint_LargeRelError(h1_R21_Alpha_Divide_3He,RelErr_Cut);
  
  TH1D* h1_R21_D_Divide_P = (TH1D*) h1_PtA_Lab[1][1]->Clone("h1_R21_D_Divide_P");
  h1_R21_D_Divide_P->Divide(h1_PtA_Lab[0][1]);
  h1_R21_D_Divide_P->Multiply(h1_PtA_Lab[0][0]);
  h1_R21_D_Divide_P->Divide(h1_PtA_Lab[1][0]);
  h1_R21_D_Divide_P->SetMarkerColor(1); h1_R21_D_Divide_P->SetLineColor(1); h1_R21_D_Divide_P->SetMarkerStyle(32);
  h1_R21_D_Divide_P->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_D_Divide_P,"D/P","lp");
  RemovePoint_LargeRelError(h1_R21_D_Divide_P,RelErr_Cut);
  
  TH1D* h1_R21_T_Divide_D = (TH1D*) h1_PtA_Lab[1][2]->Clone("h1_R21_T_Divide_D");
  h1_R21_T_Divide_D->Divide(h1_PtA_Lab[0][2]);
  h1_R21_T_Divide_D->Multiply(h1_PtA_Lab[0][1]);
  h1_R21_T_Divide_D->Divide(h1_PtA_Lab[1][1]);
  h1_R21_T_Divide_D->SetMarkerColor(1); h1_R21_T_Divide_D->SetLineColor(1); h1_R21_T_Divide_D->SetMarkerStyle(22);
  h1_R21_T_Divide_D->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_T_Divide_D,"T/D","lp");
  RemovePoint_LargeRelError(h1_R21_T_Divide_D,RelErr_Cut);
  
  //////////////////////////////////////////////////////////////////////////
  
  TH1D* h1_R21_D_Multiply_P = (TH1D*) h1_PtA_Lab[1][1]->Clone("h1_R21_D_Multiply_P");
  h1_R21_D_Multiply_P->Divide(h1_PtA_Lab[0][1]);
  h1_R21_D_Multiply_P->Multiply(h1_PtA_Lab[1][0]);
  h1_R21_D_Multiply_P->Divide(h1_PtA_Lab[0][0]);
  h1_R21_D_Multiply_P->SetMarkerColor(4); h1_R21_D_Multiply_P->SetLineColor(4); h1_R21_D_Multiply_P->SetMarkerStyle(30);
  h1_R21_D_Multiply_P->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_D_Multiply_P,"D#bulletP","lp");
  RemovePoint_LargeRelError(h1_R21_D_Multiply_P,RelErr_Cut);
  
  TH1D* h1_R21_3He = (TH1D*) h1_PtA_Lab[1][3]->Clone("h1_R21_3He"); h1_R21_3He->Divide(h1_PtA_Lab[0][3]);
  h1_R21_3He->SetMarkerColor(4); h1_R21_3He->SetLineColor(4); h1_R21_3He->SetMarkerStyle(29);
  h1_R21_3He->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_3He,"3He","lp");
  RemovePoint_LargeRelError(h1_R21_3He,RelErr_Cut);
  
  /////////////////////////////////////////////////////////////////////////
  
  TH1D* h1_R21_P = (TH1D*) h1_PtA_Lab[1][0]->Clone("h1_R21_P"); h1_R21_P->Divide(h1_PtA_Lab[0][0]);
  h1_R21_P->SetMarkerColor(2); h1_R21_P->SetLineColor(2);  h1_R21_P->SetMarkerStyle(21);
  h1_R21_P->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_P,"P","lp");
  RemovePoint_LargeRelError(h1_R21_P,RelErr_Cut);
  
  TH1D* h1_R21_Alpha_Divide_T = (TH1D*) h1_PtA_Lab[1][4]->Clone("h1_R21_Alpha_Divide_T");
  h1_R21_Alpha_Divide_T->Divide(h1_PtA_Lab[0][4]);
  h1_R21_Alpha_Divide_T->Multiply(h1_PtA_Lab[0][2]);
  h1_R21_Alpha_Divide_T->Divide(h1_PtA_Lab[1][2]);
  h1_R21_Alpha_Divide_T->SetMarkerColor(2); h1_R21_Alpha_Divide_T->SetLineColor(2); h1_R21_Alpha_Divide_T->SetMarkerStyle(25);
  h1_R21_Alpha_Divide_T->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_Alpha_Divide_T,"#alpha/T","lp");
  RemovePoint_LargeRelError(h1_R21_Alpha_Divide_T,RelErr_Cut);
  
  
  TH1D* h1_R21_3He_Divide_D = (TH1D*) h1_PtA_Lab[1][3]->Clone("h1_R21_3He_Divide_D");
  h1_R21_3He_Divide_D->Divide(h1_PtA_Lab[0][3]);
  h1_R21_3He_Divide_D->Multiply(h1_PtA_Lab[0][1]);
  h1_R21_3He_Divide_D->Divide(h1_PtA_Lab[1][1]);
  h1_R21_3He_Divide_D->SetMarkerColor(2); h1_R21_3He_Divide_D->SetLineColor(2); h1_R21_3He_Divide_D->SetMarkerStyle(24);
  h1_R21_3He_Divide_D->Draw("same");
  legend_R21_NZ_NotEqul_Cluster->AddEntry(h1_R21_3He_Divide_D,"3He/D","lp");
  RemovePoint_LargeRelError(h1_R21_3He_Divide_D,RelErr_Cut);
  
  ///////////////////////////////////////////////////////////////////////////////
  
  legend_R21_NZ_NotEqul_Cluster->Draw();
  
  //for the NZ equal
  ////////////////////////////////////////////////////////////////////
  sprintf(NameTem,"c1_R21_NZ_Equl_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_NZ_Equl = new TCanvas(NameTem,NameTem,1);
  TLegend* legend_R21_NZ_Equl_Cluster = new TLegend(0.6,0.6,0.8,0.8);
  
  TH1D* h1_R21_3He_Multiply_T = (TH1D*) h1_PtA_Lab[1][3]->Clone("h1_R21_3He_Multiply_T");
  h1_R21_3He_Multiply_T->Divide(h1_PtA_Lab[0][3]);
  h1_R21_3He_Multiply_T->Multiply(h1_PtA_Lab[1][2]);
  h1_R21_3He_Multiply_T->Divide(h1_PtA_Lab[0][2]);
  h1_R21_3He_Multiply_T->SetMarkerColor(1); h1_R21_3He_Multiply_T->SetLineColor(1); h1_R21_3He_Multiply_T->SetMarkerStyle(34);
  h1_R21_3He_Multiply_T->Draw();
  h1_R21_3He_Multiply_T->GetXaxis()->SetRangeUser(100,400);
  h1_R21_3He_Multiply_T->GetYaxis()->SetRangeUser(0.5,2.0);
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_3He_Multiply_T,"3He#bulletT","lp");
  RemovePoint_LargeRelError(h1_R21_3He_Multiply_T,RelErr_Cut);
  
  
  TH1D* h1_R21_4He_Multiply_D = (TH1D*) h1_PtA_Lab[1][4]->Clone("h1_R21_4He_Multiply_D");
  h1_R21_4He_Multiply_D->Divide(h1_PtA_Lab[0][4]);
  h1_R21_4He_Multiply_D->Multiply(h1_PtA_Lab[1][1]);
  h1_R21_4He_Multiply_D->Divide(h1_PtA_Lab[0][1]);
  h1_R21_4He_Multiply_D->SetMarkerColor(1); h1_R21_4He_Multiply_D->SetLineColor(1); h1_R21_4He_Multiply_D->SetMarkerStyle(28);
  h1_R21_4He_Multiply_D->Draw("same");
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_4He_Multiply_D,"#alpha#bulletD","lp");
  RemovePoint_LargeRelError(h1_R21_4He_Multiply_D,RelErr_Cut);
  
  ////////////////////////////////////////////////////////////////////
  
  TH1D* h1_R21_T_Multiply_P = (TH1D*) h1_PtA_Lab[1][2]->Clone("h1_R21_T_Multiply_P");
  h1_R21_T_Multiply_P->Divide(h1_PtA_Lab[0][2]);
  h1_R21_T_Multiply_P->Multiply(h1_PtA_Lab[1][0]);
  h1_R21_T_Multiply_P->Divide(h1_PtA_Lab[0][0]);
  h1_R21_T_Multiply_P->SetMarkerColor(6); h1_R21_T_Multiply_P->SetLineColor(6); h1_R21_T_Multiply_P->SetMarkerStyle(3);
  h1_R21_T_Multiply_P->Draw("same");
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_T_Multiply_P,"T#bulletP","lp");
  RemovePoint_LargeRelError(h1_R21_T_Multiply_P,RelErr_Cut);
  
  
  TH1D* h1_R21_4He = (TH1D*) h1_PtA_Lab[1][4]->Clone("h1_R21_4He"); h1_R21_4He->Divide(h1_PtA_Lab[0][4]);
  h1_R21_4He->SetMarkerColor(6); h1_R21_4He->SetLineColor(6); h1_R21_4He->SetMarkerStyle(33);
  h1_R21_4He->Draw("same");
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_4He,"#alpha","lp");
  RemovePoint_LargeRelError(h1_R21_4He,RelErr_Cut);
  
  TH1D* h1_R21_D_Multiply_D = (TH1D*) h1_PtA_Lab[1][1]->Clone("h1_R21_D_Multiply_D");
  h1_R21_D_Multiply_D->Divide(h1_PtA_Lab[0][1]);
  h1_R21_D_Multiply_D->Multiply(h1_PtA_Lab[1][1]);
  h1_R21_D_Multiply_D->Divide(h1_PtA_Lab[0][1]);
  h1_R21_D_Multiply_D->SetMarkerColor(6); h1_R21_D_Multiply_D->SetLineColor(6); h1_R21_D_Multiply_D->SetMarkerStyle(27);
  h1_R21_D_Multiply_D->Draw("same");
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_D_Multiply_D,"D#bulletD","lp");
  RemovePoint_LargeRelError(h1_R21_D_Multiply_D,RelErr_Cut);
  
  ////////////////////////////////////////////////////////////////////
  
  TH1D* h1_R21_D = (TH1D*) h1_PtA_Lab[1][1]->Clone("h1_R21_D"); h1_R21_D->Divide(h1_PtA_Lab[0][1]);
  h1_R21_D->GetXaxis()->SetRangeUser(100,400);
  h1_R21_D->GetYaxis()->SetRangeUser(0.5,2.0);
  h1_R21_D->GetYaxis()->SetTitle("R21");
  h1_R21_D->SetMarkerColor(4); h1_R21_D->SetLineColor(4); h1_R21_D->SetMarkerStyle(22);
  h1_R21_D->Draw("same");
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_D,"D","lp");
  RemovePoint_LargeRelError(h1_R21_D,RelErr_Cut);
  
  
  TH1D* h1_R21_3He_Divide_P = (TH1D*) h1_PtA_Lab[1][3]->Clone("h1_R21_3He_Divide_P");
  h1_R21_3He_Divide_P->Divide(h1_PtA_Lab[0][3]);
  h1_R21_3He_Divide_P->Multiply(h1_PtA_Lab[0][0]);
  h1_R21_3He_Divide_P->Divide(h1_PtA_Lab[1][0]);
  h1_R21_3He_Divide_P->SetMarkerColor(4); h1_R21_3He_Divide_P->SetLineColor(4); h1_R21_3He_Divide_P->SetMarkerStyle(26);
  h1_R21_3He_Divide_P->Draw("same");
  legend_R21_NZ_Equl_Cluster->AddEntry(h1_R21_3He_Divide_P,"3He/P","lp");
  RemovePoint_LargeRelError(h1_R21_3He_Divide_P,RelErr_Cut);
  
  legend_R21_NZ_Equl_Cluster->Draw();
  
  //the below is for calculating the single ratio of n-like/p-like and t-like/3He-like.
  TCanvas* c1_SR_T_3He = new TCanvas("c1_SR_T_3He","c1_DR_T_3He",1);
  c1_SR_T_3He->Divide(2,1);
  TLegend* legend_SR_T_3He = new TLegend(0.6,0.6,0.8,0.8);
  
  c1_SR_T_3He->cd(1);
  TH1D* h1_SR_T_3He_Sys1 = (TH1D*) h1_PtA_Lab[0][2]->Clone("h1_SR_T_3He_Sys1");
  h1_SR_T_3He_Sys1->Divide(h1_PtA_Lab[0][3]);
  h1_SR_T_3He_Sys1->Draw();
  h1_SR_T_3He_Sys1->GetXaxis()->SetRangeUser(150,400);
  h1_SR_T_3He_Sys1->GetYaxis()->SetRangeUser(0,3);
  TH1D* h1_SR_T_3He_Sys2 = (TH1D*) h1_PtA_Lab[1][2]->Clone("h1_SR_T_3He_Sys2");
  h1_SR_T_3He_Sys2->Divide(h1_PtA_Lab[1][3]);
  h1_SR_T_3He_Sys2->Draw("same");
  Set_1D_HistoStyle(h1_SR_T_3He_Sys1, "Pt/A(MeV/u/c)", "SR(T/3He)", 20, 2);
  Set_1D_HistoStyle(h1_SR_T_3He_Sys2, "Pt/A(MeV/u/c)", "SR(T/3He)", 20, 1);
  
  
  c1_SR_T_3He->cd(2);
  TH1D* h1_SR_D_P2_Sys1 = (TH1D*) h1_PtA_Lab[0][1]->Clone("h1_SR_D_P2_Sys1");
  h1_SR_D_P2_Sys1->Divide(h1_PtA_Lab[0][0]);
  h1_SR_D_P2_Sys1->Divide(h1_PtA_Lab[0][0]);
  h1_SR_D_P2_Sys1->Draw();
  h1_SR_D_P2_Sys1->GetXaxis()->SetRangeUser(150,450);
  h1_SR_D_P2_Sys1->GetYaxis()->SetRangeUser(0,50);
  TH1D* h1_SR_D_P2_Sys2 = (TH1D*) h1_PtA_Lab[1][1]->Clone("h1_SR_D_P2_Sys2");
  h1_SR_D_P2_Sys2->Divide(h1_PtA_Lab[1][0]);
  h1_SR_D_P2_Sys2->Divide(h1_PtA_Lab[1][0]);
  h1_SR_D_P2_Sys2->Draw("same");
  Set_1D_HistoStyle(h1_SR_D_P2_Sys1, "Pt/A(MeV/u/c)", "SR(D/P^{2})", 20, 2);
  Set_1D_HistoStyle(h1_SR_D_P2_Sys2, "Pt/A(MeV/u/c)", "SR(D/P^{2})", 20, 1);

  //legend_SR_T_3He->AddEntry(h1_SR_T_3He_Sys1,"T/3He","p");
  
  
  //the below is for calculating the double ratio of n-like/p-like and t-like/3He-like.
  TCanvas* c1_DR_T_3He = new TCanvas("c1_DR_T_3He","c1_DR_T_3He",1);
  
  TLegend* legend_DR_np_T_3He = new TLegend(0.6,0.6,0.8,0.8);
  
  TH1D* h1_DR_T_3He = (TH1D*) h1_R21_T->Clone("h1_DR_T_3He");
  h1_DR_T_3He->Divide(h1_R21_3He);
  h1_DR_T_3He->Draw();
  legend_DR_np_T_3He->AddEntry(h1_DR_T_3He,"T/3He","p");
  
  TH1D* h1_DR_Alpha_P_D_P = (TH1D*) h1_R21_Alpha_Divide_P->Clone("h1_DR_Alpha_P_D_P");
  h1_DR_Alpha_P_D_P->Divide(h1_R21_D_Multiply_P);
  h1_DR_Alpha_P_D_P->Draw("same");
  legend_DR_np_T_3He->AddEntry(h1_DR_Alpha_P_D_P,"#frac{#alpha/P}{D#bulletP}","p");
  
  TH1D* h1_DR_T_D_Alpha_T = (TH1D*) h1_R21_T_Divide_D->Clone("h1_DR_T_D_Alpha_T");
  h1_DR_T_D_Alpha_T->Divide(h1_R21_Alpha_Divide_T);
  h1_DR_T_D_Alpha_T->Draw("same");
  legend_DR_np_T_3He->AddEntry(h1_DR_T_D_Alpha_T,"#frac{T/D}{#alpha/T}","p");
  
  TH1D* h1_DR_D_P_P = (TH1D*) h1_R21_D_Divide_P->Clone("h1_DR_D_P_P");
  h1_DR_D_P_P->Divide(h1_R21[0]);
  h1_DR_D_P_P->Draw("same");
  legend_DR_np_T_3He->AddEntry(h1_DR_D_P_P,"#frac{D/P}{P}","p");
  legend_DR_np_T_3He->Draw();
  
  //The below is for calculating the Alpha and Beta.
  double Z[] = {1,1,1,2,2,2}; double N[] = {0,1,2,1,2,4};
  int PtA_Num = 0;
  double PtA[100] = {0}; double PtA_Err[100] = {0}; 
  double Alpha[100] = {0}; double Alpha_Err[100] = {0}; double Beta[100] = {0}; double Beta_Err[100] = {0};
  double C[100] = {0}; double C_Err[100] = {0};
  double Alpha_Plus_Beta[100] = {0}; double Alpha_Plus_Beta_Err[100] = {0};
  double Alpha_Subtract_Beta[100] = {0}; double Alpha_Subtract_Beta_Err[100] = {0};
  
  sprintf(NameTem,"c1_N_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_N = new TCanvas(NameTem,NameTem,1);
  c1_N->Divide(3,4);
  sprintf(NameTem,"c1_Z_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Z = new TCanvas(NameTem,NameTem,1);
  c1_Z->Divide(3,4);
  
  //int FirstBinIndex = h1_R21[0]->FindFirstBinAbove(0);
  int XBinNum = h1_R21[0]->GetNbinsX();
  for(int i=1;i<=XBinNum;i++) 
  {
    int IsValue = 1; //if each particle have R21, then the alpha and beta will be calculated.
    
    double PtA_R21[6] = {0}; double PtA_R21_Err[6] = {0};
    for(int iP=0;iP<ParticleNum-1;iP++)
    {
      PtA_R21[iP] = h1_R21[iP]->GetBinContent(i);
      if(Is_R21Err_Applied==1) { PtA_R21_Err[iP] = h1_R21[iP]->GetBinError(i); }
      else { PtA_R21_Err[iP] = 0.02*PtA_R21[iP]; }
      if(PtA_R21[iP]==0) { IsValue=0; }
    }
    
    if(IsValue==1) 
    {
      PtA[PtA_Num] = h1_R21[0]->GetBinCenter(i);
      cout<<endl<<endl<<" --> Pt/A(MeV/c) : "<<PtA[PtA_Num]<<endl;
      cout.width(7);
      cout<<"      "<<"R21"<<"  "<<"R21_Err"<<endl;
      for(int iP=0;iP<ParticleNum-1;iP++)
      {
        cout.width(7);
        cout<<ParticleName[iP]<<"  "<<PtA_R21[iP]<<"  "<<PtA_R21_Err[iP]<<endl;
      }
      
      Cal_Alpha_Beta(PtA_Num, c1_N, c1_Z, ParticleNum-1, Z, N, PtA_R21, PtA_R21_Err, Alpha, Alpha_Err, Beta, Beta_Err, C, C_Err);
      Alpha_Plus_Beta[PtA_Num] = Alpha[PtA_Num]+Beta[PtA_Num];
      Alpha_Plus_Beta_Err[PtA_Num] = Sqrt(Power(Alpha_Err[PtA_Num],2)+Power(Beta_Err[PtA_Num],2));
      Alpha_Subtract_Beta[PtA_Num] = Alpha[PtA_Num]-Beta[PtA_Num];
      Alpha_Subtract_Beta_Err[PtA_Num] = Sqrt(Power(Alpha_Err[PtA_Num],2)+Power(Beta_Err[PtA_Num],2));
      PtA_Num++; 
    }
    
    
  }
  
  cout<<"PtA Num : "<<PtA_Num<<endl;
  
  TGraphErrors* gr1_Alpha = new TGraphErrors(PtA_Num,PtA,Alpha,PtA_Err,Alpha_Err);
  TGraphErrors* gr1_Beta = new TGraphErrors(PtA_Num,PtA,Beta,PtA_Err,Beta_Err);
  TGraphErrors* gr1_Alpha_Plus_Beta = new TGraphErrors(PtA_Num,PtA,Alpha_Plus_Beta,PtA_Err,Alpha_Plus_Beta_Err);
  TGraphErrors* gr1_Alpha_Subtract_Beta = new TGraphErrors(PtA_Num,PtA,Alpha_Subtract_Beta,PtA_Err,Alpha_Subtract_Beta_Err);
  
  sprintf(NameTem,"c1_Alpha_Beta_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Alpha_Beta = new TCanvas(NameTem,NameTem,1);
  
  gr1_Alpha->Draw("A*"); gr1_Alpha->SetMarkerStyle(20); gr1_Alpha->SetMarkerColor(2);
  gr1_Beta->Draw("same*"); gr1_Beta->SetMarkerStyle(21); gr1_Beta->SetMarkerColor(4);
  gr1_Alpha_Plus_Beta->Draw("same*"); gr1_Alpha_Plus_Beta->SetMarkerStyle(22); gr1_Alpha_Plus_Beta->SetMarkerColor(6);
  //gr1_Alpha_Subtract_Beta->Draw("same*"); gr1_Alpha_Subtract_Beta->SetMarkerStyle(23); gr1_Alpha_Subtract_Beta->SetMarkerColor(6);
  gr1_Alpha->GetYaxis()->SetRangeUser(-0.4,0.4);
  
  Set_1D_GraphStyle(gr1_Alpha,"Pt/A(MeV/c)", "#Delta #mu/T");
  Set_1D_GraphStyle(gr1_Beta,"Pt/A(MeV/c)", "#Delta #mu/T");
  Set_1D_GraphStyle(gr1_Alpha_Plus_Beta,"Pt/A(MeV/c)", "#Delta #mu/T");
  gr1_Alpha->SetMarkerStyle(20); gr1_Alpha->SetMarkerColor(1); 
  gr1_Beta->SetMarkerStyle(24); gr1_Beta->SetMarkerColor(1);
  gr1_Alpha_Plus_Beta->SetMarkerStyle(34); gr1_Alpha_Plus_Beta->SetMarkerColor(1);
  
  TLegend* legend_Alpha_Beta = new TLegend(0.6,0.6,0.8,0.8);
  legend_Alpha_Beta->SetNColumns(3);
  legend_Alpha_Beta->AddEntry(gr1_Alpha,"#Delta#mu_{n}/T","p");
  legend_Alpha_Beta->AddEntry(gr1_Beta,"#Delta#mu_{p}/T","p");
  legend_Alpha_Beta->AddEntry(gr1_Alpha_Plus_Beta,"(#Delta#mu_{n}+#Delta#mu_{p})/T","p");
  //legend_Alpha_Beta->AddEntry(gr1_Alpha_Plus_Beta,"(#Delta#mu_{p}+#Delta#mu_{n})/T","p");
  //legend_Alpha_Beta->AddEntry(gr1_Alpha_Subtract_Beta,"(#Delta#mu_{p}-#Delta#mu_{n})/T","p");
  legend_Alpha_Beta->Draw();
  
  //the below is for calculating the R21 from the alpha and beta.
  //c1_R21->cd();
  c1_Pt_A_Sys1->cd(2);
  TGraph* gr1_R21_Alpha_Beta[ParticleNum] = {0};
  for(int iP=0;iP<ParticleNum-1;iP++)
  {
    double R21_Particle[100] = {0};
    for(int i=0;i<PtA_Num;i++)
    {
      R21_Particle[i] = Exp((Alpha[i]*N[iP]+Beta[i]*Z[iP]+C[i]));
    }
    gr1_R21_Alpha_Beta[iP] = new TGraph(PtA_Num,PtA,R21_Particle);
    gr1_R21_Alpha_Beta[iP]->Draw("sameL");
    gr1_R21_Alpha_Beta[iP]->SetLineColor(LineColor[iP]);
    gr1_R21_Alpha_Beta[iP]->SetLineWidth(2);
    gr1_R21_Alpha_Beta[iP]->SetLineStyle(2);
  }
  
}

TH1D* Select_Pt_A(string SystemName, string ParticleName, TH2D* h2_Pt_A_Rapidity,double* Range)
{
  int XBin1 = h2_Pt_A_Rapidity->GetXaxis()->FindBin(Range[0]);
  int XBin2 = h2_Pt_A_Rapidity->GetXaxis()->FindBin(Range[1]);
  
  string h2_Name = h2_Pt_A_Rapidity->GetName();
  h2_Name = h2_Name.substr(2);
  char NameTem[200];
  sprintf(NameTem,"h1_%s_%s_Pt_A_y_%.1f_%.1f",SystemName.c_str(), ParticleName.c_str(), Range[0],Range[1]);
  
  TH1D* h1_Pt_A = (TH1D*) h2_Pt_A_Rapidity->ProjectionY(NameTem,XBin1,XBin2);
  h1_Pt_A->Scale(1.0/(Range[1]-Range[0])); //the rapidity normalization.
  Set_1D_HistoStyle(h1_Pt_A,"Pt/A(MeV/u/c)","d^{2}M/d_{Pt}d_{y}",20,1);
  cout<<" Produce " << h1_Pt_A->GetName()<<endl;
  
return h1_Pt_A;
}



TF1* Extend_Pt_A(string SystemName, string ParticleName, double* RapidityRange, TH1D* h1_Pt_A, double* FitRange, double* ExtendRange)
{
  char NameTem[200];
  sprintf(NameTem,"f1_Pt_A_%s_%s_y_%.1f_%.1f",SystemName.c_str(),ParticleName.c_str(),RapidityRange[0],RapidityRange[1]);
  TF1* f1_Pt_A = new TF1(NameTem,"TMath::Exp(-[0]*TMath::Power(x-[1],2)+[2])",200,600);
  h1_Pt_A->Fit(NameTem,"","",FitRange[0],FitRange[1]);
  
  int Bin1 = h1_Pt_A->FindBin(ExtendRange[0]);
  int Bin2 = h1_Pt_A->FindBin(ExtendRange[1]);
  
  for(int i=Bin1;i<=Bin2;i++)
  {
    double BinCenter = h1_Pt_A->GetBinCenter(i);
    double Value = f1_Pt_A->Eval(BinCenter);
    h1_Pt_A->SetBinContent(i,Value);
    h1_Pt_A->SetBinError(i,Value*0.1); //now, this 0.1 error is just an a
  }
  
  f1_Pt_A->SetLineColor(4);
  
return f1_Pt_A;
}

TH1D* Get_CIp(string NameTem, int ParticleNum,TH1D** h1_PtA_Lab)
{
  TH1D* h1_CIp_PtA_Lab = (TH1D*) h1_PtA_Lab[0]->Clone(NameTem.c_str());
  h1_CIp_PtA_Lab->Add(h1_PtA_Lab[1],1);
  h1_CIp_PtA_Lab->Add(h1_PtA_Lab[2],1);
  h1_CIp_PtA_Lab->Add(h1_PtA_Lab[3],2);
  h1_CIp_PtA_Lab->Add(h1_PtA_Lab[4],2);
  h1_CIp_PtA_Lab->Add(h1_PtA_Lab[5],2);
return h1_CIp_PtA_Lab;
}

TH1D* Get_Ratio_T_3He(string NameTem,TH1D* h1_PtA_Lab_T, TH1D* h1_PtA_Lab_3He)
{
  TH1D* h1_Ratio_T_3He = (TH1D*) h1_PtA_Lab_T->Clone(NameTem.c_str());
  h1_Ratio_T_3He->Divide(h1_PtA_Lab_3He);
return h1_Ratio_T_3He;
}

TH1D* Get_PseudoNeutron(string NameTem,TH1D* h1_Ratio_T_3He, TH1D* h1_PtA_Lab_P)
{
  TH1D* h1_PseudoNeutron_PtA_Lab = (TH1D*) h1_PtA_Lab_P->Clone(NameTem.c_str());
  h1_PseudoNeutron_PtA_Lab->Multiply(h1_Ratio_T_3He);
  
  return h1_PseudoNeutron_PtA_Lab;
}

TH1D* Get_CIn(string NameTem, int ParticleNum,TH1D** h1_PtA_Lab,TH1D* h1_PseudoNeutron)
{
  TH1D* h1_CIn_PtA_Lab = (TH1D*) h1_PtA_Lab[1]->Clone(NameTem.c_str());
  h1_CIn_PtA_Lab->Add(h1_PtA_Lab[2],2);
  h1_CIn_PtA_Lab->Add(h1_PtA_Lab[3],1);
  h1_CIn_PtA_Lab->Add(h1_PtA_Lab[4],2);
  h1_CIn_PtA_Lab->Add(h1_PtA_Lab[5],4);
  h1_CIn_PtA_Lab->Add(h1_PseudoNeutron,1);
  
return h1_CIn_PtA_Lab;
}

TH1D* GetTemperature(string NameTem,TH1D** h1_PtA_Lab)
{
  TH1D* h1_Temperature = (TH1D*) h1_PtA_Lab[1]->Clone(NameTem.c_str()); //D
  h1_Temperature->Multiply(h1_PtA_Lab[4]); //4He
  h1_Temperature->Divide(h1_PtA_Lab[2]); //T
  h1_Temperature->Divide(h1_PtA_Lab[3]); //3He
  int BinNum = h1_Temperature->GetNbinsX();
  for(int i=1;i<=BinNum;i++)
  {
    double BinContent = h1_Temperature->GetBinContent(i);
    if(BinContent<=0) { continue; }
    double T = 14.3/Log(1.59*BinContent);
    h1_Temperature->SetBinContent(i,T);
  }
return h1_Temperature;
}

void RemovePoint_LargeRelError(TH1D* h1,double RelErr)
{
  int BinNum = h1->GetNbinsX();
  for(int i=1;i<=BinNum;i++)
  {
    double BinContent = h1->GetBinContent(i);
    double BinContentErr = h1->GetBinError(i);
    if((BinContent!=0 && Abs(BinContentErr/BinContent)>RelErr) || (BinContent==0 && Abs(BinContentErr)>RelErr))
    {
      h1->SetBinContent(i,0);
      h1->SetBinError(i,0);
    }
  }
}

void Cal_Alpha_Beta(int Index, TCanvas* c1_N, TCanvas* c1_Z, int ParticleNum, double* Z, double* N, double* R21, double* R21_Err,double* Alpha, double* Alpha_Err, double* Beta, double* Beta_Err, double* C, double* C_Err)
{
  double Z_Err[20] = {0}; double N_Err[20] = {0};
  double Ln_R21[20] = {0}; double Ln_R21_Err[20] = {0};
  for(int i=0;i<ParticleNum;i++)
  { 
    Z_Err[i] = 0; N_Err[i] = 0;
    Ln_R21[i] = Log(R21[i]); 
    Ln_R21_Err[i] = R21_Err[i]/R21[i];
  }
  
  char NameTem[200];
  sprintf(NameTem,"gr2_%d",Index);
  TGraph2DErrors* gr2 = new TGraph2DErrors(5);
  gr2->SetName(NameTem);
  for(int i=0;i<ParticleNum;i++)
  {
    gr2->SetPoint(i,N[i],Z[i],Ln_R21[i]);
    gr2->SetPointError(i,N_Err[i],Z_Err[i],Ln_R21_Err[i]); //attention: if the Ln_R21_Err[i] is set to 0, then it will show a error/bug.
  }
  
  sprintf(NameTem,"f2_%d",Index);
  TF2* f2 = new TF2(NameTem,"[0]*x+[1]*y+[2]",-1,3,1,3); //3 parameters
  f2->SetParameters(0.3,-0.3,0.1);
  gr2->Fit(f2);
  TF2 *fit2 = (TF2*)gr2->FindObject(NameTem);
  double Results[3]; double ResultsErr[3];
  f2->GetParameters(Results);
  ResultsErr[0] = f2->GetParError(0); ResultsErr[1] = f2->GetParError(1); ResultsErr[2] = f2->GetParError(2);
  Alpha[Index] = Results[0]; Beta[Index] = Results[1]; C[Index] = Results[2];
  Alpha_Err[Index] = ResultsErr[0]; Beta_Err[Index] = ResultsErr[1]; C_Err[Index] = ResultsErr[2];
  
  //the below is for drawing the subset of function.
  
  TGraphErrors* gr1_R21_N = new TGraphErrors(ParticleNum,N,R21,N_Err,R21_Err);
  TGraphErrors* gr1_R21_Z = new TGraphErrors(ParticleNum,Z,R21,Z_Err,R21_Err);
  
  c1_N->cd(Index+1)->SetLogy(1);
  gr1_R21_N->Draw("A*");
  gr1_R21_N->GetXaxis()->SetLimits(0,6);
  gr1_R21_N->GetYaxis()->SetRangeUser(0.5,2);
  
  sprintf(NameTem,"TMath::Exp((x*%.3f+%.3f*%d+%.3f))",Results[0],Results[1],1,Results[2]);
  TF1* f1_R21_Z_H = new TF1("f1_R21_Z_H",NameTem,-0.2,2.2);
  f1_R21_Z_H->Draw("same");
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d+%.3f)",Results[0],Results[1],2,Results[2]);
  TF1* f1_R21_Z_He = new TF1("f1_R21_Z_He",NameTem,0.8,4.2);
  f1_R21_Z_He->Draw("same");
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d+%.3f)",Results[0],Results[1],3,Results[2]);
  TF1* f1_R21_Z_Li = new TF1("f1_R21_Z_Li",NameTem,2.8,5.2);
  f1_R21_Z_Li->Draw("same"); f1_R21_Z_Li->SetLineStyle(7);
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d+%.3f)",Results[0],Results[1],4,Results[2]);
  TF1* f1_R21_Z_Be = new TF1("f1_R21_Z_Be",NameTem,2.8,6.2);
  f1_R21_Z_Be->Draw("same"); f1_R21_Z_Be->SetLineStyle(7);
  
  c1_Z->cd(Index+1)->SetLogy(1);
  gr1_R21_Z->Draw("A*");
  gr1_R21_Z->GetXaxis()->SetLimits(0,6);
  gr1_R21_Z->GetYaxis()->SetRangeUser(0.5,2);
  
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x+%.3f)",1,Results[0],Results[1],Results[2]);
  TF1* f1_R21_N1_Z1 = new TF1("f1_R21_N1_Z1",NameTem,0.8,2.2);
  f1_R21_N1_Z1->Draw("same");
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x+%.3f)",2,Results[0],Results[1],Results[2]);
  TF1* f1_R21_N2_Z1 = new TF1("f1_R21_N2_Z1",NameTem,0.8,2.2);
  f1_R21_N2_Z1->Draw("same");
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x+%.3f)",4,Results[0],Results[1],Results[2]);
  TF1* f1_R21_N4_Z2 = new TF1("f1_R21_N4_Z2",NameTem,1.8,3.2);
  f1_R21_N4_Z2->Draw("same"); f1_R21_N4_Z2->SetLineStyle(7);
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x+%.3f)",3,Results[0],Results[1],Results[2]);
  TF1* f1_R21_N3_Z3 = new TF1("f1_R21_N3_Z3",NameTem,2.8,4.2);
  f1_R21_N3_Z3->Draw("same"); f1_R21_N3_Z3->SetLineStyle(7);
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x+%.3f)",5,Results[0],Results[1],Results[2]);
  TF1* f1_R21_N5_Z3 = new TF1("f1_R21_N5_Z3",NameTem,2.8,4.2);
  f1_R21_N5_Z3->Draw("same"); f1_R21_N5_Z3->SetLineStyle(7);
}

void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle, int MarkerStyle, int MarkerColor)
{
  h1->GetXaxis()->SetTitle(XTitle.c_str());
  h1->GetYaxis()->SetTitle(YTitle.c_str());
  
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetLabelSize(0.05);
  
  h1->SetLineWidth(2);
  h1->SetMarkerSize(1.5);
  h1->SetMarkerColor(MarkerColor);
  h1->SetMarkerStyle(MarkerStyle);
}

void Set_1D_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle)
{
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();
  gr1->GetXaxis()->SetTitle(XTitle.c_str());
  gr1->GetYaxis()->SetTitle(YTitle.c_str());
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetXaxis()->SetTitleOffset(0.8);
  gr1->GetYaxis()->SetTitleOffset(0.8);
  gr1->GetXaxis()->SetLabelSize(0.05);
  gr1->GetYaxis()->SetLabelSize(0.05);
  gr1->SetLineWidth(2);
  gr1->SetMarkerSize(1.5);
}
