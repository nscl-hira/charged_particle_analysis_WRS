#include "TMath.h"
using namespace TMath;

//this script will draw all the different particle together.
void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle);
void Set_2D_HistoStyle(TH2D* h2, string XTitle, string YTitle);
void Set_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle);
TH2D* Get_Rapidity_PtA(string h2_Rapidity_PtA_Name, TH2D* h2_Ekin_Theta_Lab, int ParticleA, double ParticleMass, double Rapidity_Beam_inLabFrame,double* ThetaCut, double* EkinCut);
TH1D* Select_Pt_A(string SystemName, string ParticleName, TH2D* h2_Pt_A_Rapidity,double* Range,bool Is_PtA_Theta_RawCount);
void RemovePoint_LargeRelError(TH1D* h1,double RelErr);
void Cal_Avg(TH1D* h1_Ratio,double* Avg_Range,double* Avg_Sys);
void Scale_Uncertainty(TH1D* h1, double* ScaleFactor);
void PrintPointValue(TH1D* h1);
TH1D* Rebin_Histo(TH1D* h1); // this will change according to temporary decision.

void Cal_StdDev_Histo(TH1D* h1_Avg,int HistoNum,TH1D** h1_DR_Array,int PtA_PointNum, double* PtA_Array,double* StdDev);

void ProduceHistoErr(TH1D* h1_PtA_Lab_Err);

void _4_Cal_DR_CIn_CIp_Rebin()
{
  gStyle->SetErrorX(0);
  
  string SystemTag_1 = "Ca40Ni58E140"; 
  string SystemTag_2 = "Ca48Ni64E140";
  
  bool Is_d_p2 = 1;
  
  Exp_RunInfo* RunInfo = new Exp_RunInfo();
  double Rapidity_Beam_inLabFrame[2] = { RunInfo->Get_BeamRapidity_Lab(SystemTag_1), RunInfo->Get_BeamRapidity_Lab(SystemTag_2)};
  Hira_ReactionLost* Hira_ReactionLost_Corrector = new Hira_ReactionLost();
  bool Is_Corr_ReactionLost = 0; //Attention: this correction is based on the (Ekin, Theta_lab)
  
  string ConditionFolder = "Noise_Theta_30_75";
  double RelErr_Cut = 0.25;

  bool IsExtend = 0;
  int RebinNum = 30;
  //2 approaches can be used to produce the Pt/A: from Ekin_Theta_Lab or from the merged Pt/A.
  //if from mergered Pt/A, the reaction lost correction cannot be applied.
  bool Is_Start_From_Ekin_ThetaLab = 1;
  bool Is_PtA_Theta_RawCount = 0; //this rawcount means no GeoEff and no normalized. This option is only used to study the statistic error.
  
  int iEff = 1;
  int ib = 0;
  double RapidityRange[2] = {0.4,0.6};
  
  const int EffOpt = 3;
  string EffCor_Tag[] = {"_noEff", "_GeoEff",  "_GeoReactionEff"};
  
  double Mass_1u = 931.49410242;
  const int ParticleNum = 5;
  string ParticleName[ParticleNum] = {"p","d","t","3He","4He"};
  int ParticleZ[ParticleNum] = {1,1,1,2,2};
  int ParticleN[ParticleNum] = {0,1,2,1,2};
  double ParticleMass[ParticleNum] = {938.272,1876.179,2809.514, //H
                      2809.496,3728.510}; //unit: MeV/c2
  
  double ThetaCut[ParticleNum][2] = {{30,75},{30,75},{30,75},
                                     {30,75},{30,75}};
                                        
  //double EkinCut[ParticleNum][2] = {{20,198.0/1.0},{15,263.0/2.0},{12,312.0/3.0},{20,200},{18,200}};
  //Large Attention: here 3He is using the Ekin_cut as Triton for a reasonable T/3He
  //double EkinCut[ParticleNum][2] = {{20,198.0/1.0},{15,263.0/2.0},{12,312.0/3.0},{20,312.0/3.0},{18,200}}; 
  //double EkinCut[ParticleNum][2] = {{20,312.0/3.0},{15,312.0/3.0},{12,312.0/3.0},{20,312.0/3.0},{18,312.0/3.0}}; 
  double EkinCut[ParticleNum][2] = {{20,263.0/2.0},{15,263.0/2.0},{12,312.0/3.0},{20,312.0/3.0},{18,312.0/3.0}}; 
  
  int LineColor[] = {1,1,1,2,2,2};
  int MarkerStyle[] = {20,24,21,20,24,21};
  
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
  
  TH2D* h2_Ekin_Theta_Lab_Normalized_GeoEff[2][ParticleNum];
  TH2D* h2_y_PtA_Lab[2][ParticleNum];
  TH1D* h1_PtA_Lab[2][ParticleNum];
  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    //for system 1
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    cout<<"Read : "<<NameTem<<endl;
    f1_MergedESpec_1->cd();
    h2_Ekin_Theta_Lab_Normalized_GeoEff[0][iP] = (TH2D*) f1_MergedESpec_1->Get(NameTem);
    sprintf(NameTem,"h2_%s_%s_%s_Ekin_Theta_Lab%s_bhat",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
    h2_Ekin_Theta_Lab_Normalized_GeoEff[0][iP]->SetName(NameTem);
    
    if(Is_Corr_ReactionLost==1) { Hira_ReactionLost_Corrector->CorrESpec_perA(h2_Ekin_Theta_Lab_Normalized_GeoEff[0][iP], ParticleZ[iP], ParticleZ[iP]+ParticleN[iP]); }
    
    if(Is_Start_From_Ekin_ThetaLab==1)
    {
 sprintf(NameTem,"h2_%s_%s_%s_y_PtA_Lab%s_bhat",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
      h2_y_PtA_Lab[0][iP] = Get_Rapidity_PtA(NameTem, h2_Ekin_Theta_Lab_Normalized_GeoEff[0][iP],(ParticleZ[iP]+ParticleN[iP]),ParticleMass[iP],Rapidity_Beam_inLabFrame[0],ThetaCut[iP], EkinCut[iP]);;
    }
    else
    {
      if(Is_PtA_Theta_RawCount==1)
      { sprintf(NameTem,"h2_Rapidity_PtA_Lab_RawCount_%s",ParticleName[iP].c_str()); }
      else
      { sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str()); }
      h2_y_PtA_Lab[0][iP] = (TH2D*) f1_MergedESpec_1->Get(NameTem);
      sprintf(NameTem,"h2_%s_%s_%s_y_PtA_Lab%s_bhat",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
      h2_y_PtA_Lab[0][iP]->SetName(NameTem);
    }
    
    //for system 2
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    f1_MergedESpec_2->cd();
    h2_Ekin_Theta_Lab_Normalized_GeoEff[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem);
    sprintf(NameTem,"h2_%s_%s_%s_Ekin_Theta_Lab%s_bhat",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
    h2_Ekin_Theta_Lab_Normalized_GeoEff[1][iP]->SetName(NameTem);
    if(Is_Corr_ReactionLost==1) { Hira_ReactionLost_Corrector->CorrESpec_perA(h2_Ekin_Theta_Lab_Normalized_GeoEff[1][iP], ParticleZ[iP], ParticleZ[iP]+ParticleN[iP]); }
    
    if(Is_Start_From_Ekin_ThetaLab==1)
    {
sprintf(NameTem,"h2_%s_%s_%s_y_PtA_Lab%s_bhat",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
      h2_y_PtA_Lab[1][iP] = Get_Rapidity_PtA(NameTem, h2_Ekin_Theta_Lab_Normalized_GeoEff[1][iP],(ParticleZ[iP]+ParticleN[iP]),ParticleMass[iP],Rapidity_Beam_inLabFrame[1],ThetaCut[iP], EkinCut[iP]);
    }
    else
    {
      if(Is_PtA_Theta_RawCount==1)
      { sprintf(NameTem,"h2_Rapidity_PtA_Lab_RawCount_%s",ParticleName[iP].c_str());}
      else
      { sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str()); }
      
      h2_y_PtA_Lab[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem); 
sprintf(NameTem,"h2_%s_%s_%s_y_PtA_Lab%s_bhat",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
      h2_y_PtA_Lab[1][iP]->SetName(NameTem);
    }
    
    
    //the below is for projecting Pt_A.
    f1_MergedESpec_1->cd();
    h1_PtA_Lab[0][iP] = Select_Pt_A(SystemTag_1, ParticleName[iP], h2_y_PtA_Lab[0][iP],RapidityRange,Is_PtA_Theta_RawCount);
    f1_MergedESpec_2->cd();
    h1_PtA_Lab[1][iP] = Select_Pt_A(SystemTag_2, ParticleName[iP], h2_y_PtA_Lab[1][iP],RapidityRange,Is_PtA_Theta_RawCount);
    
    h1_PtA_Lab[0][iP]->SetLineColor(LineColor[iP]);
    h1_PtA_Lab[1][iP]->SetLineColor(LineColor[iP]);
    h1_PtA_Lab[0][iP]->SetMarkerColor(LineColor[iP]);
    h1_PtA_Lab[1][iP]->SetMarkerColor(LineColor[iP]);
    h1_PtA_Lab[0][iP]->SetMarkerStyle(MarkerStyle[iP]);
    h1_PtA_Lab[1][iP]->SetMarkerStyle(MarkerStyle[iP]);
    
    //the below is for rebin the data.
    
    h1_PtA_Lab[0][iP]->Rebin(RebinNum); //the error will be calculated by the sqrt(x1^2 + x2^2).
    if(Is_PtA_Theta_RawCount!=1) {h1_PtA_Lab[0][iP]->Scale(1.0/RebinNum); } //if need or any confusing rise, please check the error bar calculation again.
    h1_PtA_Lab[1][iP]->Rebin(RebinNum);
    if(Is_PtA_Theta_RawCount!=1) { h1_PtA_Lab[1][iP]->Scale(1.0/RebinNum); }
    
    //Modify the last 3 points
    h1_PtA_Lab[0][iP] = Rebin_Histo(h1_PtA_Lab[0][iP]);
    h1_PtA_Lab[1][iP] = Rebin_Histo(h1_PtA_Lab[1][iP]);
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
  //for system1
  sprintf(NameTem,"c1_Pt_A_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Sys1 = new TCanvas(NameTem,NameTem,1);
  sprintf(NameTem,"c1_Pt_A_Err_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Err_Sys1 = new TCanvas(NameTem,NameTem,1);
  //for system2
  sprintf(NameTem,"c1_Pt_A_%s_%s_%s_y_%.1f_%.1f",SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Sys2 = new TCanvas(NameTem,NameTem,1);
  sprintf(NameTem,"c1_Pt_A_Err_%s_%s_%s_y_%.1f_%.1f",SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),RapidityRange[0],RapidityRange[1]);
  TCanvas* c1_Pt_A_Err_Sys2 = new TCanvas(NameTem,NameTem,1);
  
  //Draw the Pt/A and error of Pt/A on the canvas
  TLegend* legend_PID = new TLegend(0.6,0.6,0.8,0.8);
  for(int i=0;i<ParticleNum;i++)
  {
    TH1D* h1_PtA_Lab_Err_Sys1 = (TH1D*) h1_PtA_Lab[0][i]->Clone();
    ProduceHistoErr(h1_PtA_Lab_Err_Sys1);
    TH1D* h1_PtA_Lab_Err_Sys2 = (TH1D*) h1_PtA_Lab[1][i]->Clone();
    ProduceHistoErr(h1_PtA_Lab_Err_Sys2);
    
    if(i==0) 
    { 
      c1_Pt_A_Sys1->cd()->SetGrid(1,1);
      h1_PtA_Lab[0][i]->Draw();
      h1_PtA_Lab[0][i]->GetXaxis()->SetRangeUser(100,400);
      
      c1_Pt_A_Err_Sys1->cd()->SetGrid(1,1);
      h1_PtA_Lab_Err_Sys1->Draw();
      h1_PtA_Lab_Err_Sys1->GetXaxis()->SetRangeUser(100,400);
    }
    else 
    { 
      c1_Pt_A_Sys1->cd();
      h1_PtA_Lab[0][i]->Draw("same"); 
      
      c1_Pt_A_Err_Sys1->cd();
      h1_PtA_Lab_Err_Sys1->Draw("same");
    }
    
    
    if(i==0) 
    { 
      c1_Pt_A_Sys2->cd()->SetGrid(1,1);
      h1_PtA_Lab[1][i]->Draw(); 
      h1_PtA_Lab[1][i]->GetXaxis()->SetRangeUser(100,400);
      
      c1_Pt_A_Err_Sys2->cd()->SetGrid(1,1);
      h1_PtA_Lab_Err_Sys2->Draw();
      h1_PtA_Lab_Err_Sys2->GetXaxis()->SetRangeUser(100,400);
    }
    else 
    { 
      c1_Pt_A_Sys2->cd();
      h1_PtA_Lab[1][i]->Draw("same");
      
      c1_Pt_A_Err_Sys2->cd();
      h1_PtA_Lab_Err_Sys2->Draw("same");
    }
    cout<<"*********"<<SystemTag_1<<"  "<<ParticleName[i]<<endl;
    PrintPointValue(h1_PtA_Lab[0][i]);
    cout<<endl;
    cout<<"*********"<<SystemTag_2<<"  "<<ParticleName[i]<<endl;
    PrintPointValue(h1_PtA_Lab[1][i]);
    cout<<endl;
    
    legend_PID->AddEntry(h1_PtA_Lab[1][i],ParticleName[i].c_str(),"lp");
  }
  c1_Pt_A_Sys1->cd();  legend_PID->Draw();
  c1_Pt_A_Sys2->cd();  legend_PID->Draw();
  c1_Pt_A_Err_Sys1->cd();  legend_PID->Draw();
  c1_Pt_A_Err_Sys2->cd();  legend_PID->Draw();

  
  TH1D* h1_SR_T_3He[2];
  TH1D* h1_SR_D_P2[2];
  sprintf(NameTem,"h1_%s_SR_T_3He",SystemTag_1.c_str());
  h1_SR_T_3He[0] = (TH1D*) h1_PtA_Lab[0][2]->Clone(NameTem);
  h1_SR_T_3He[0]->Divide(h1_PtA_Lab[0][3]);
  sprintf(NameTem,"h1_%s_SR_T_3He",SystemTag_2.c_str());
  h1_SR_T_3He[1] = (TH1D*) h1_PtA_Lab[1][2]->Clone(NameTem);
  h1_SR_T_3He[1]->Divide(h1_PtA_Lab[1][3]);
  
  sprintf(NameTem,"h1_%s_SR_D_P2",SystemTag_1.c_str());
  h1_SR_D_P2[0] = (TH1D*) h1_PtA_Lab[0][1]->Clone(NameTem);
  h1_SR_D_P2[0]->Divide(h1_PtA_Lab[0][0]);
  h1_SR_D_P2[0]->Divide(h1_PtA_Lab[0][0]);
  sprintf(NameTem,"h1_%s_SR_D_P2",SystemTag_2.c_str());
  h1_SR_D_P2[1] = (TH1D*) h1_PtA_Lab[1][1]->Clone(NameTem);
  h1_SR_D_P2[1]->Divide(h1_PtA_Lab[1][0]);
  h1_SR_D_P2[1]->Divide(h1_PtA_Lab[1][0]);
  
  sprintf(NameTem,"c1_%s_%s_SR_T_3He",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_SR_T_3He = new TCanvas(NameTem,NameTem,1);
  c1_SR_T_3He->Divide(2,2);
  c1_SR_T_3He->cd(1); 
  h1_SR_T_3He[0]->Draw(); h1_SR_T_3He[1]->Draw("same"); 
  Set_1D_HistoStyle(h1_SR_T_3He[0], "Pt/A(MeV/c)", "SR(T/3He)");
  Set_1D_HistoStyle(h1_SR_T_3He[1], "Pt/A(MeV/c)", "SR(T/3He)");
  h1_SR_T_3He[0]->SetMarkerColor(2);
  h1_SR_T_3He[1]->SetMarkerColor(1);
  
  c1_SR_T_3He->cd(2);
  h1_SR_D_P2[0]->Draw(); h1_SR_D_P2[1]->Draw("same"); 
  Set_1D_HistoStyle(h1_SR_D_P2[0], "Pt/A(MeV/c)", "SR(D/P^{2})");
  Set_1D_HistoStyle(h1_SR_D_P2[1], "Pt/A(MeV/c)", "SR(D/P^{2})");
  h1_SR_D_P2[0]->SetMarkerColor(2);
  h1_SR_D_P2[1]->SetMarkerColor(1);
  
  c1_SR_T_3He->cd(3);
  TH1D* h1_Ratio_T_3He_D_P_P[2];
  h1_Ratio_T_3He_D_P_P[0] = (TH1D*) h1_SR_T_3He[0]->Clone("h1_Ratio_T_3He_D_P_P_Sys1");
  h1_Ratio_T_3He_D_P_P[0]->Divide(h1_SR_D_P2[0]);
  h1_Ratio_T_3He_D_P_P[0]->Draw();
  double Avg_Sys1[2] = {0};
  double Avg_Range[2] = {180,350};
  Cal_Avg(h1_Ratio_T_3He_D_P_P[0],Avg_Range,Avg_Sys1);
  cout<<SystemTag_1<<" : (t/3He)/(d/p2) = "<<Avg_Sys1[0]<<" StdDev = "<<Avg_Sys1[1]<<endl;
  
  h1_Ratio_T_3He_D_P_P[1] = (TH1D*) h1_SR_T_3He[1]->Clone("h1_Ratio_T_3He_D_P_P_Sys2");
  h1_Ratio_T_3He_D_P_P[1]->Divide(h1_SR_D_P2[1]);
  h1_Ratio_T_3He_D_P_P[1]->Draw("same");
  double Avg_Sys2[2] = {0};
  Cal_Avg(h1_Ratio_T_3He_D_P_P[1],Avg_Range,Avg_Sys2);
  cout<<SystemTag_2<<" : (t/3He)/(d/p2) = "<<Avg_Sys2[0]<<" StdDev = "<<Avg_Sys2[1]<<endl;
  
  TH1D* h1_Pseudo_n[2];
  sprintf(NameTem,"h1_%s_pseudo_n",SystemTag_1.c_str());
  if(Is_d_p2!=1) 
  {
    h1_Pseudo_n[0] = (TH1D*) h1_PtA_Lab[0][0]->Clone(NameTem);
    h1_Pseudo_n[0]->Multiply(h1_SR_T_3He[0]); 
  }
  else
  {
    h1_Pseudo_n[0] = (TH1D*) h1_PtA_Lab[0][1]->Clone(NameTem);
    h1_Pseudo_n[0]->Divide(h1_PtA_Lab[0][0]);
    //h1_Pseudo_n[0]->Scale(Avg_Sys1[0]);
    Scale_Uncertainty(h1_Pseudo_n[0],Avg_Sys1);
  }
  
  sprintf(NameTem,"h1_%s_pseudo_n",SystemTag_2.c_str());
  if(Is_d_p2!=1) 
  {
    h1_Pseudo_n[1] = (TH1D*) h1_PtA_Lab[1][0]->Clone(NameTem);
    h1_Pseudo_n[1]->Multiply(h1_SR_T_3He[1]);
  }
  else
  {
    h1_Pseudo_n[1] = (TH1D*) h1_PtA_Lab[1][1]->Clone(NameTem);
    h1_Pseudo_n[1]->Divide(h1_PtA_Lab[1][0]);
    //h1_Pseudo_n[1]->Scale(Avg_Sys2[0]);
    Scale_Uncertainty(h1_Pseudo_n[1],Avg_Sys2);
  }
  
  sprintf(NameTem,"c1_%s_%s_pseudo_n",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_pseudo_n = new TCanvas(NameTem,NameTem,1);
  c1_pseudo_n->cd();
  h1_Pseudo_n[0]->Draw(); h1_Pseudo_n[1]->Draw("same");
  Set_1D_HistoStyle(h1_Pseudo_n[0], "Pt/A(MeV/c)", "d^{2}M/dydPt");
  Set_1D_HistoStyle(h1_Pseudo_n[1], "Pt/A(MeV/c)", "d^{2}M/dydPt");
  h1_Pseudo_n[0]->SetMarkerColor(2);
  h1_Pseudo_n[1]->SetMarkerColor(1);
  
  TH1D* h1_CI_n[2]; TH1D* h1_CI_p[2];
  sprintf(NameTem,"h1_%s_CI_n",SystemTag_1.c_str());
  h1_CI_n[0] = (TH1D*) h1_Pseudo_n[0]->Clone(NameTem);
  h1_CI_n[0]->Add(h1_PtA_Lab[0][1],1); h1_CI_n[0]->Add(h1_PtA_Lab[0][2],2); 
  h1_CI_n[0]->Add(h1_PtA_Lab[0][3],1); h1_CI_n[0]->Add(h1_PtA_Lab[0][4],2);
  
  sprintf(NameTem,"h1_%s_CI_n",SystemTag_2.c_str());
  h1_CI_n[1] = (TH1D*) h1_Pseudo_n[1]->Clone(NameTem);
  h1_CI_n[1]->Add(h1_PtA_Lab[1][1],1); h1_CI_n[1]->Add(h1_PtA_Lab[1][2],2); 
  h1_CI_n[1]->Add(h1_PtA_Lab[1][3],1); h1_CI_n[1]->Add(h1_PtA_Lab[1][4],2);
  
  sprintf(NameTem,"h1_%s_CI_p",SystemTag_1.c_str());
  h1_CI_p[0] = (TH1D*) h1_PtA_Lab[0][0]->Clone(NameTem);
  h1_CI_p[0]->Add(h1_PtA_Lab[0][1],1); h1_CI_p[0]->Add(h1_PtA_Lab[0][2],1); 
  h1_CI_p[0]->Add(h1_PtA_Lab[0][3],2); h1_CI_p[0]->Add(h1_PtA_Lab[0][4],2);
  
  sprintf(NameTem,"h1_%s_CI_p",SystemTag_2.c_str());
  h1_CI_p[1] = (TH1D*) h1_PtA_Lab[1][0]->Clone(NameTem);
  h1_CI_p[1]->Add(h1_PtA_Lab[1][1],1); h1_CI_p[1]->Add(h1_PtA_Lab[1][2],1); 
  h1_CI_p[1]->Add(h1_PtA_Lab[1][3],2); h1_CI_p[1]->Add(h1_PtA_Lab[1][4],2);
  
  sprintf(NameTem,"c1_%s_%s_CI_n",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_CI_n = new TCanvas(NameTem,NameTem,1);
  h1_CI_n[0]->Draw(); h1_CI_n[1]->Draw("same");
  Set_1D_HistoStyle(h1_CI_n[0], "Pt/A(MeV/c)", "d^{2}M/dydPt");
  Set_1D_HistoStyle(h1_CI_n[1], "Pt/A(MeV/c)", "d^{2}M/dydPt");
  h1_CI_n[0]->SetMarkerColor(2);
  h1_CI_n[1]->SetMarkerColor(1);
  
  sprintf(NameTem,"c1_%s_%s_CI_p",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_CI_p = new TCanvas(NameTem,NameTem,1);
  h1_CI_p[0]->Draw(); h1_CI_p[1]->Draw("same");
  Set_1D_HistoStyle(h1_CI_p[0], "Pt/A(MeV/c)", "d^{2}M/dydPt");
  Set_1D_HistoStyle(h1_CI_p[1], "Pt/A(MeV/c)", "d^{2}M/dydPt");
  h1_CI_p[0]->SetMarkerColor(2);
  h1_CI_p[1]->SetMarkerColor(1);
  
  sprintf(NameTem,"h1_%s_%s_DR_CIn_CIp",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_CIn_CIp = (TH1D*) h1_CI_n[1]->Clone(NameTem);
  h1_DR_CIn_CIp->Divide(h1_CI_n[0]);
  h1_DR_CIn_CIp->Multiply(h1_CI_p[0]);
  h1_DR_CIn_CIp->Divide(h1_CI_p[1]);
  
  sprintf(NameTem,"c1_%s_%s_DR_CIn_CIp",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_DR_CIn_CIp = new TCanvas(NameTem,NameTem,1);
  h1_DR_CIn_CIp->Draw();
  Set_1D_HistoStyle(h1_DR_CIn_CIp, "Pt/A(MeV/c)", "DR(CI-n/CI-p)");
  RemovePoint_LargeRelError(h1_DR_CIn_CIp, RelErr_Cut);
  PrintPointValue(h1_DR_CIn_CIp);
  
  //the below is for calculate the DR(n-like/p-like)
  sprintf(NameTem,"c1_%s_%s_DR_nlike_plike",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_DR_nlike_plike = new TCanvas(NameTem,NameTem,1);
  
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_Avg",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_Avg = (TH1D*) h1_PtA_Lab[0][0]->Clone(NameTem);
  h1_DR_nlike_plike_Avg->Reset();
  TLegend* legend_DR_nlike_plike = new TLegend(0.6,0.6,0.8,0.8);
  legend_DR_nlike_plike->SetNColumns(3);
  
  

  //T/3He
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_T_3He",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_T_3He = (TH1D*) h1_PtA_Lab[1][2]->Clone(NameTem);
  h1_DR_nlike_plike_T_3He->Divide(h1_PtA_Lab[0][2]);
  h1_DR_nlike_plike_T_3He->Divide(h1_PtA_Lab[1][3]);
  h1_DR_nlike_plike_T_3He->Multiply(h1_PtA_Lab[0][3]);
  h1_DR_nlike_plike_T_3He->Draw();
  h1_DR_nlike_plike_T_3He->SetMarkerStyle(20);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_T_3He,"#frac{T}{^{3}He}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_T_3He, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  h1_DR_nlike_plike_T_3He->SetMarkerColor(2);
  h1_DR_nlike_plike_T_3He->GetYaxis()->SetRangeUser(1,2);
  h1_DR_nlike_plike_T_3He->GetXaxis()->SetRangeUser(120,375);
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_T_3He);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_T_3He, RelErr_Cut);
  
  
  //(D/P)/P
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_D_P_P",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_D_P_P = (TH1D*) h1_PtA_Lab[1][1]->Clone(NameTem);
  h1_DR_nlike_plike_D_P_P->Divide(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_D_P_P->Divide(h1_PtA_Lab[1][0]); h1_DR_nlike_plike_D_P_P->Divide(h1_PtA_Lab[1][0]);
  h1_DR_nlike_plike_D_P_P->Multiply(h1_PtA_Lab[0][0]); h1_DR_nlike_plike_D_P_P->Multiply(h1_PtA_Lab[0][0]);
  h1_DR_nlike_plike_D_P_P->Draw("same");
  h1_DR_nlike_plike_D_P_P->SetMarkerStyle(24);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_D_P_P,"#frac{D/P}{P}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_D_P_P, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_D_P_P);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_D_P_P, RelErr_Cut);
  
  //(T/D)/P
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_T_D_P",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_T_D_P = (TH1D*) h1_PtA_Lab[1][2]->Clone(NameTem);
  h1_DR_nlike_plike_T_D_P->Divide(h1_PtA_Lab[0][2]);
  h1_DR_nlike_plike_T_D_P->Divide(h1_PtA_Lab[1][0]); h1_DR_nlike_plike_T_D_P->Divide(h1_PtA_Lab[1][1]);
  h1_DR_nlike_plike_T_D_P->Multiply(h1_PtA_Lab[0][0]); h1_DR_nlike_plike_T_D_P->Multiply(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_T_D_P->Draw("same");
  h1_DR_nlike_plike_T_D_P->SetMarkerStyle(25);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_T_D_P,"#frac{T/D}{P}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_T_D_P, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_T_D_P);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_T_D_P, RelErr_Cut);
  
  //(4He/3He)/(3He/D)
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_4He_3He_3He_D",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_4He_3He_3He_D = (TH1D*) h1_PtA_Lab[1][4]->Clone(NameTem);
  h1_DR_nlike_plike_4He_3He_3He_D->Divide(h1_PtA_Lab[0][4]);
  h1_DR_nlike_plike_4He_3He_3He_D->Multiply(h1_PtA_Lab[1][1]);
  h1_DR_nlike_plike_4He_3He_3He_D->Divide(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_4He_3He_3He_D->Multiply(h1_PtA_Lab[0][3]); h1_DR_nlike_plike_4He_3He_3He_D->Multiply(h1_PtA_Lab[0][3]);
  h1_DR_nlike_plike_4He_3He_3He_D->Divide(h1_PtA_Lab[1][3]); h1_DR_nlike_plike_4He_3He_3He_D->Divide(h1_PtA_Lab[1][3]);
  h1_DR_nlike_plike_4He_3He_3He_D->Draw("same");
  h1_DR_nlike_plike_4He_3He_3He_D->SetMarkerStyle(26);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_4He_3He_3He_D,"#frac{^{4}He/^{3}He}{^{3}He/D}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_4He_3He_3He_D, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_4He_3He_3He_D);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_4He_3He_3He_D, RelErr_Cut);
  
  //(4He/3He)/P
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_4He_3He_P",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_4He_3He_P = (TH1D*) h1_PtA_Lab[1][4]->Clone(NameTem);
  h1_DR_nlike_plike_4He_3He_P->Divide(h1_PtA_Lab[1][3]);
  h1_DR_nlike_plike_4He_3He_P->Divide(h1_PtA_Lab[1][0]);
  h1_DR_nlike_plike_4He_3He_P->Multiply(h1_PtA_Lab[0][3]);
  h1_DR_nlike_plike_4He_3He_P->Multiply(h1_PtA_Lab[0][0]);
  h1_DR_nlike_plike_4He_3He_P->Divide(h1_PtA_Lab[0][4]);
  h1_DR_nlike_plike_4He_3He_P->Draw("same");
  h1_DR_nlike_plike_4He_3He_P->SetMarkerStyle(27);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_4He_3He_P,"#frac{^{4}He/^{3}He}{P}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_4He_3He_P, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_4He_3He_P);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_4He_3He_P, RelErr_Cut);
  
  //(D/P)/(4He/T)
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_D_P_4He_T",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_D_P_4He_T = (TH1D*) h1_PtA_Lab[1][1]->Clone(NameTem);
  h1_DR_nlike_plike_D_P_4He_T->Multiply(h1_PtA_Lab[1][2]);
  h1_DR_nlike_plike_D_P_4He_T->Divide(h1_PtA_Lab[1][0]);
  h1_DR_nlike_plike_D_P_4He_T->Divide(h1_PtA_Lab[1][4]);
  h1_DR_nlike_plike_D_P_4He_T->Multiply(h1_PtA_Lab[0][0]);
  h1_DR_nlike_plike_D_P_4He_T->Multiply(h1_PtA_Lab[0][4]);
  h1_DR_nlike_plike_D_P_4He_T->Divide(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_D_P_4He_T->Divide(h1_PtA_Lab[0][2]);
  h1_DR_nlike_plike_D_P_4He_T->Draw("same");
  h1_DR_nlike_plike_D_P_4He_T->SetMarkerStyle(28);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_D_P_4He_T,"#frac{D/P}{^{4}He/T}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_D_P_4He_T, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_D_P_4He_T);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_D_P_4He_T, RelErr_Cut);
  
  
  //(D/P)/(3He/D)
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_D_P_3He_D",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_D_P_3He_D = (TH1D*) h1_PtA_Lab[1][1]->Clone(NameTem);
  h1_DR_nlike_plike_D_P_3He_D->Multiply(h1_PtA_Lab[1][1]);
  h1_DR_nlike_plike_D_P_3He_D->Divide(h1_PtA_Lab[1][0]);
  h1_DR_nlike_plike_D_P_3He_D->Divide(h1_PtA_Lab[1][3]);
  h1_DR_nlike_plike_D_P_3He_D->Multiply(h1_PtA_Lab[0][0]);
  h1_DR_nlike_plike_D_P_3He_D->Multiply(h1_PtA_Lab[0][3]);
  h1_DR_nlike_plike_D_P_3He_D->Divide(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_D_P_3He_D->Divide(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_D_P_3He_D->Draw("same");
  h1_DR_nlike_plike_D_P_3He_D->SetMarkerStyle(30);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_D_P_3He_D,"#frac{D/P}{^{3}He/D}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_D_P_3He_D, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_D_P_3He_D);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_D_P_3He_D, RelErr_Cut);
  
  
  //(T/D)/(4He/T)
  sprintf(NameTem,"h1_%s_%s_DR_nlike_plike_T_D_4He_T",SystemTag_1.c_str(),SystemTag_2.c_str());
  TH1D* h1_DR_nlike_plike_T_D_4He_T = (TH1D*) h1_PtA_Lab[1][2]->Clone(NameTem);
  h1_DR_nlike_plike_T_D_4He_T->Multiply(h1_PtA_Lab[1][2]);
  h1_DR_nlike_plike_T_D_4He_T->Divide(h1_PtA_Lab[1][1]);
  h1_DR_nlike_plike_T_D_4He_T->Divide(h1_PtA_Lab[1][4]);
  h1_DR_nlike_plike_T_D_4He_T->Multiply(h1_PtA_Lab[0][1]);
  h1_DR_nlike_plike_T_D_4He_T->Multiply(h1_PtA_Lab[0][4]);
  h1_DR_nlike_plike_T_D_4He_T->Divide(h1_PtA_Lab[0][2]);
  h1_DR_nlike_plike_T_D_4He_T->Divide(h1_PtA_Lab[0][2]);
  h1_DR_nlike_plike_T_D_4He_T->Draw("same");
  h1_DR_nlike_plike_T_D_4He_T->SetMarkerStyle(32);
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_T_D_4He_T,"#frac{T/D}{^{4}He/T}","p");
  Set_1D_HistoStyle(h1_DR_nlike_plike_T_D_4He_T, "Pt/A(MeV/c)", "DR(n-like/p-like)");
  
  h1_DR_nlike_plike_Avg->Add(h1_DR_nlike_plike_T_D_4He_T);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_T_D_4He_T, RelErr_Cut);
  
  h1_DR_nlike_plike_Avg->Scale(1.0/8.0);
  RemovePoint_LargeRelError(h1_DR_nlike_plike_Avg, RelErr_Cut);
  h1_DR_nlike_plike_Avg->Draw("sameE5");
  legend_DR_nlike_plike->AddEntry(h1_DR_nlike_plike_Avg,"Average","p");
  h1_DR_nlike_plike_Avg->SetLineWidth(3);
  h1_DR_nlike_plike_Avg->SetLineColor(6);
  h1_DR_nlike_plike_Avg->SetMarkerStyle(21);
  h1_DR_nlike_plike_Avg->SetMarkerColor(6);
//  h1_DR_nlike_plike_Avg->SetFillColor(6);
  h1_DR_nlike_plike_Avg->SetFillColorAlpha(6, 0.35);
//  h1_DR_nlike_plike_Avg->SetFillStyle(4050);
  //the below is for calculating the StdDev
  int HistoNum = 8;
  TH1D* h1_DR_Array[8];
  h1_DR_Array[0] = h1_DR_nlike_plike_T_3He;
  h1_DR_Array[1] = h1_DR_nlike_plike_D_P_P;
  h1_DR_Array[2] = h1_DR_nlike_plike_T_D_P;
  h1_DR_Array[3] = h1_DR_nlike_plike_4He_3He_3He_D;
  h1_DR_Array[4] = h1_DR_nlike_plike_4He_3He_P;
  h1_DR_Array[5] = h1_DR_nlike_plike_D_P_4He_T;
  h1_DR_Array[6] = h1_DR_nlike_plike_D_P_3He_D;
  h1_DR_Array[7] = h1_DR_nlike_plike_T_D_4He_T;
  
  double PtA_Start = 130; //MeV/c
  int PtA_PointNum = 13; 
  double PtA_Array[100] = {0};
  for(int i=0;i<PtA_PointNum;i++) { PtA_Array[i] = PtA_Start+i*RebinNum; }
  double StdDevValue[100] = {0};
  Cal_StdDev_Histo(h1_DR_nlike_plike_Avg,HistoNum,h1_DR_Array,PtA_PointNum,PtA_Array,StdDevValue);
//  Set_Err_Histo(h1_DR_nlike_plike_Avg,PtA_Array,StdDev);
  
  legend_DR_nlike_plike->Draw();
  
  c1_SR_T_3He->cd(4);
  h1_DR_nlike_plike_T_3He->Draw();
  h1_DR_nlike_plike_D_P_P->Draw("same");
  
}

TH2D* Get_Rapidity_PtA(string h2_Rapidity_PtA_Name, TH2D* h2_Ekin_Theta_Lab, int ParticleA, double ParticleMass, double Rapidity_Beam_inLabFrame,double* ThetaCut, double* EkinCut)
{
  int BinXNum = h2_Ekin_Theta_Lab->GetNbinsX();
  int BinYNum = h2_Ekin_Theta_Lab->GetNbinsY();
  
  TAxis* X_Axis = h2_Ekin_Theta_Lab->GetXaxis();
  TAxis* Y_Axis = h2_Ekin_Theta_Lab->GetYaxis();
  
  double XBinSize = X_Axis->GetBinCenter(2)-X_Axis->GetBinCenter(1);
  double YBinSize = Y_Axis->GetBinCenter(2)-Y_Axis->GetBinCenter(1);
  
  TH2D* h2_Rapidity_PtA_tem = new TH2D(h2_Rapidity_PtA_Name.c_str(),";y/y_{Beam-Lab};Pt/A(MeV/c)",100,0,1,600,0,600);
  TAxis* X_Axis_Rapidity = h2_Rapidity_PtA_tem->GetXaxis();
  TAxis* Y_Axis_PtA = h2_Rapidity_PtA_tem->GetYaxis();
  
  for(int iY=1;iY<BinYNum;iY++)
  {
    for(int iX=1;iX<BinYNum;iX++)
    {
      double Ekin = X_Axis->GetBinCenter(iX)+XBinSize*gRandom->Uniform(-0.5,0.5); //this random does not affect the yield and ratio.
      double Theta = Y_Axis->GetBinCenter(iY)+YBinSize*gRandom->Uniform(-0.5,0.5); //this random does not affect the yield and ratio.
      
      if(Theta<ThetaCut[0] || Theta>ThetaCut[1]) { continue; } //here is make a cut on Ekin/A
      if(Ekin<EkinCut[0] || Ekin>EkinCut[1]) { continue; }
      
      Ekin = Ekin*ParticleA;
      
      double BinContent = h2_Ekin_Theta_Lab->GetBinContent(iX,iY);
      double BinErr = h2_Ekin_Theta_Lab->GetBinError(iX,iY);
      
      double P_Mag_Lab = Sqrt(Ekin*Ekin+2*Ekin*ParticleMass);
      TVector3 P_3D_Lab(0,0,P_Mag_Lab);
      P_3D_Lab.SetTheta(Theta*DegToRad());
      P_3D_Lab.SetPhi(0.0*DegToRad());
      
      double Pt_Lab = P_3D_Lab.Pt();//get the transverse momentum.
      double PtA_Lab = Pt_Lab/ParticleA;
      
      TLorentzVector P_4D_CM(P_3D_Lab,ParticleMass+Ekin);
      double Rapidity_Lab = P_4D_CM.Rapidity()/Rapidity_Beam_inLabFrame;
      
      int Rapidity_Bin = X_Axis_Rapidity->FindBin(Rapidity_Lab);
      int PtA_Bin = Y_Axis_PtA->FindBin(PtA_Lab);
      double BinError_Rapidity_PtA = h2_Rapidity_PtA_tem->GetBinError(Rapidity_Bin,PtA_Bin);
      
      h2_Rapidity_PtA_tem->Fill(Rapidity_Lab,PtA_Lab,BinContent);
      h2_Rapidity_PtA_tem->SetBinError(Rapidity_Bin,PtA_Bin,Sqrt(BinErr*BinErr+BinError_Rapidity_PtA*BinError_Rapidity_PtA));
    }
  }
return h2_Rapidity_PtA_tem;
}

TH1D* Select_Pt_A(string SystemName, string ParticleName, TH2D* h2_Pt_A_Rapidity,double* Range, bool Is_PtA_Theta_RawCount)
{
  int XBin1 = h2_Pt_A_Rapidity->GetXaxis()->FindBin(Range[0]);
  int XBin2 = h2_Pt_A_Rapidity->GetXaxis()->FindBin(Range[1]);
  
  string h2_Name = h2_Pt_A_Rapidity->GetName();
  h2_Name = h2_Name.substr(2);
  char NameTem[200];
  sprintf(NameTem,"h1_%s_%s_Pt_A_y_%.1f_%.1f",SystemName.c_str(), ParticleName.c_str(), Range[0],Range[1]);
  
  TH1D* h1_Pt_A = (TH1D*) h2_Pt_A_Rapidity->ProjectionY(NameTem,XBin1,XBin2);
  if(Is_PtA_Theta_RawCount!=1) { h1_Pt_A->Scale(1.0/(Range[1]-Range[0])); } //the rapidity normalization.
  Set_1D_HistoStyle(h1_Pt_A,"Pt/A(MeV/u/c)","d^{2}M/d_{Pt}d_{y}");
  cout<<" Produce " << h1_Pt_A->GetName()<<endl;
  
return h1_Pt_A;
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

void ProduceHistoErr(TH1D* h1)
{
  int BinNum = h1->GetNbinsX();
  for(int i=1;i<=BinNum;i++)
  {
    double BinContent = h1->GetBinContent(i);
    double BinContentErr = h1->GetBinError(i);
    if(BinContent!=0)
    {
      h1->SetBinContent(i,BinContentErr/BinContent);
      h1->SetBinError(i,0);
    }
  }
}

void Cal_StdDev_Histo(TH1D* h1_Avg,int HistoNum,TH1D** h1_DR_Array,int PtA_PointNum,double* PtA_Array,double* StdDevValue)
{
  double Value[100][100];
  
  for(int i=0;i<HistoNum;i++)
  {
    for(int iPoint=0;iPoint<PtA_PointNum;iPoint++)
    {
      int Bin1 = h1_DR_Array[i]->FindBin(PtA_Array[iPoint]);
      double Content = h1_DR_Array[i]->GetBinContent(Bin1);
      Value[iPoint][i] = Content;
    }
  }
  
  for(int iPoint=0;iPoint<PtA_PointNum;iPoint++)
  {
    StdDevValue[iPoint] = StdDev(HistoNum,Value[iPoint]);
    cout<<StdDevValue[iPoint]<<"  ";//确认StdDev的计算方法。
  }
  cout<<endl;
  
  cout<<"Set the StdDev as the error of the averaged histogram"<<endl;
  for(int iPoint=0;iPoint<PtA_PointNum;iPoint++)
  {
    int Bin1 = h1_Avg->FindBin(PtA_Array[iPoint]);
    h1_Avg->SetBinError(Bin1,StdDevValue[iPoint]);
  }
}

void Cal_Avg(TH1D* h1_Ratio,double* Avg_Range,double* Avg_Sys)
{
  int Bin1 = h1_Ratio->FindBin(Avg_Range[0]);
  int Bin2 = h1_Ratio->FindBin(Avg_Range[1]);
  cout<<"Bin1 : "<<Bin1<<" Bin2 : "<<Bin2<<endl;
  
  
  int Num = 0;
  double Avg = 0;
  double Value[1000] = {0};
  for(int i=Bin1;i<=Bin2;i++)
  {
    Value[Num] = h1_Ratio->GetBinContent(i);
    Avg += Value[Num];
    Num++;
  }
  
  Avg = Avg/Num;
  double StdDev = RMS(Num,Value);
  
  Avg_Sys[0] = Avg;
  Avg_Sys[1] = StdDev;
}

void Scale_Uncertainty(TH1D* h1, double* ScaleFactor)
{
  int Nbin = h1->GetNbinsX();
  for(int i=1;i<=Nbin;i++)
  {
    double Value = h1->GetBinContent(i);
    double Err = h1->GetBinError(i);
    
    double ScaleValue = Value*ScaleFactor[0];
    double ScaleValueErr = 0;
    if(Value!=0) { ScaleValueErr = ScaleValue*Sqrt( Power(Err/Value,2) + Power(ScaleFactor[1]/ScaleFactor[0],2) ); }
    
    h1->SetBinContent(i,ScaleValue);
    h1->SetBinError(i,ScaleValueErr);
  }
}

void PrintPointValue(TH1D* h1)
{
  int Nbin = h1->GetNbinsX();
  for(int i=1;i<=Nbin;i++)
  {
    double BinCenter = h1->GetBinCenter(i);
    double Value = h1->GetBinContent(i);
    double Err = h1->GetBinError(i);
    
    cout<<i<<"  "<<BinCenter<<"  "<<Value<<"  "<<Err<<endl;
  }
}

TH1D* Rebin_Histo(TH1D* h1)
{
  int Nbin = h1->GetNbinsX();
  int BinIndex = h1->FindBin(405)-1;
  string HistoName = h1->GetName();
  HistoName = HistoName+"_RebinLastPoints";
  
  double XBins_New[10000] = {0};
  for(int i=1;i<=BinIndex;i++)
  {
    XBins_New[i-1] = h1->GetBinLowEdge(i);
  }
  XBins_New[BinIndex] = h1->GetBinLowEdge(BinIndex)+3*h1->GetBinWidth(Nbin);
  
  for(int i=BinIndex+1;i<=Nbin;i++)
  {
    XBins_New[i] = XBins_New[i-1]+h1->GetBinWidth(Nbin);
  }
  
  TH1D* h1_tem = (TH1D*) h1->Rebin(Nbin-2,HistoName.c_str(),XBins_New);
  double Value_Tem = h1_tem->GetBinContent(BinIndex); 
  double Err_Tem = h1_tem->GetBinError(BinIndex);
  h1_tem->SetBinContent(BinIndex,Value_Tem/3.0);
  h1_tem->SetBinError(BinIndex,Err_Tem/3.0);

return h1_tem;
}

void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle)
{
  h1->GetXaxis()->SetTitle(XTitle.c_str());
  h1->GetYaxis()->SetTitle(YTitle.c_str());
  
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetLabelSize(0.05);
  
  h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetYaxis()->SetTitleOffset(1.0);
  
  h1->SetLineWidth(2);
  h1->SetLineColor(1);
  h1->SetMarkerSize(1.5);
  h1->SetMarkerColor(1);
//  h1->SetMarkerStyle(5);
}

void Set_2D_HistoStyle(TH2D* h2, string XTitle, string YTitle)
{
  if(XTitle!="") { h2->GetXaxis()->SetTitle(XTitle.c_str()); }
  if(YTitle!="") { h2->GetYaxis()->SetTitle(YTitle.c_str()); }
  
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
}

void Set_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle)
{
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();
  gr1->GetXaxis()->SetTitle(XTitle.c_str());
  gr1->GetYaxis()->SetTitle(YTitle.c_str());
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetXaxis()->SetTitleOffset(0.8);
  gr1->GetYaxis()->SetTitleOffset(1.3);
  gr1->GetXaxis()->SetLabelSize(0.05);
  gr1->GetYaxis()->SetLabelSize(0.05);
  //gr1->SetLineWidth(2);
  gr1->SetMarkerSize(1.5);
  
  gr1->GetXaxis()->SetLimits(-0.5,0.5);
  gr1->GetYaxis()->SetRangeUser(-0.5,0.5);
}
