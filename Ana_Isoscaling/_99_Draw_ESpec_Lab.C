#include "TMath.h"
using namespace TMath;

void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle);
TH1D* Select(TH2D* h2_Ekin_Theta_Lab, double Theta, string HistoName, double DeltaTheta);

void _99_Draw_ESpec_Lab()
{
  string SystemTag_1 = "Ca40Ni58E140"; 
  string SystemTag_2 = "Ca48Ni64E140"; 
  
  string ThetaTag = "Noise_Theta_30_75";
  string EffCor_Tag = "GeoEff";
  string ImpactParaTag = "bHat_0.00_0.40";
  
  const int ParticleNum = 6;
  string ParticleName[]={"p","d","t","3He","4He","6He","6Li","7Li"};
  int LineColor[] = {1,1,1,2,2,2};
  int MarkerStyle[] = {20,24,21,20,24,21};
  
  TH2D* h2_Ekin_Theta_Lab[2][ParticleNum];
  
  char NameTem[200];
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_1.c_str(),ImpactParaTag.c_str());
  TFile* f1_MergedESpec_1 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_2.c_str(),ImpactParaTag.c_str());
  TFile* f1_MergedESpec_2 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;
  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Sum_GeoEff_%s",ParticleName[iP].c_str());
    cout<<"Read : "<<NameTem<<endl;
    f1_MergedESpec_1->cd();
    h2_Ekin_Theta_Lab[0][iP] = (TH2D*) f1_MergedESpec_1->Get(NameTem);
    sprintf(NameTem,"h2_Ekin_Theta_Lab_%s_%s_%s_%s",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_Ekin_Theta_Lab[0][iP]->SetName(NameTem);
    
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Sum_GeoEff_%s",ParticleName[iP].c_str());
    f1_MergedESpec_2->cd();
    h2_Ekin_Theta_Lab[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem);
    sprintf(NameTem,"h2_Ekin_Theta_Lab_%s_%s_%s_%s",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_Ekin_Theta_Lab[1][iP]->SetName(NameTem);
  }
  
  sprintf(NameTem,"c1_Ekin_Theta_%s_%s_%s",SystemTag_1.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Ekin_Theta_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_Ekin_Theta_Sys1->Divide(3,2);
  
  sprintf(NameTem,"c1_Ekin_Theta_%s_%s_%s",SystemTag_2.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Ekin_Theta_Sys2 = new TCanvas(NameTem,NameTem,1);
  c1_Ekin_Theta_Sys2->Divide(3,2);
  
  for(int i=0;i<ParticleNum;i++)
  {
    c1_Ekin_Theta_Sys1->cd(i+1)->SetGrid(1,1);
    h2_Ekin_Theta_Lab[0][i]->Draw("colz");
    
    c1_Ekin_Theta_Sys2->cd(i+1)->SetGrid(1,1);
    h2_Ekin_Theta_Lab[1][i]->Draw("colz");
  }
  
  //Select energy spectrum in the lab
  const int AngleNum = 5;
  double AngleValue[AngleNum] = {30,40,50,60,70};
  int AngleLineColor[] = {1,2,4,6,7};
  
  TH1D* h1_Proton_ESpec_Lab[2][AngleNum];
  TH1D* h1_Triton_ESpec_Lab[2][AngleNum];
  
  double DeltaTheta = 2; //Deg.
  for(int i=0;i<AngleNum;i++)
  {
    sprintf(NameTem,"h1_Proton_ESpec_Lab_%s_Theta%.1f",SystemTag_1.c_str(),AngleValue[i]);
    h1_Proton_ESpec_Lab[0][i] = Select(h2_Ekin_Theta_Lab[0][0],AngleValue[i],NameTem,DeltaTheta);
    Set_1D_HistoStyle(h1_Proton_ESpec_Lab[0][i], "Ekin(MeV)", "d^{2}M/d#OmegadE(Sr^{-1}MeV^{-1})");
    
    sprintf(NameTem,"h1_Triton_ESpec_Lab_%s_Theta%.1f",SystemTag_1.c_str(),AngleValue[i]);
    h1_Triton_ESpec_Lab[0][i] = Select(h2_Ekin_Theta_Lab[0][2],AngleValue[i],NameTem,DeltaTheta);
    Set_1D_HistoStyle(h1_Triton_ESpec_Lab[0][i], "Ekin/A(MeV)", "d^{2}M/d#OmegadE(Sr^{-1}MeV^{-1})");
    
    sprintf(NameTem,"h1_Proton_ESpec_Lab_%s_Theta%.1f",SystemTag_2.c_str(),AngleValue[i]);
    h1_Proton_ESpec_Lab[1][i] = Select(h2_Ekin_Theta_Lab[1][0],AngleValue[i],NameTem,DeltaTheta);
    Set_1D_HistoStyle(h1_Proton_ESpec_Lab[1][i], "Ekin(MeV)", "d^{2}M/d#OmegadE(Sr^{-1}MeV^{-1})");
    
    sprintf(NameTem,"h1_Triton_ESpec_Lab_%s_Theta%.1f",SystemTag_2.c_str(),AngleValue[i]);
    h1_Triton_ESpec_Lab[1][i] = Select(h2_Ekin_Theta_Lab[1][2],AngleValue[i],NameTem,DeltaTheta);
    Set_1D_HistoStyle(h1_Triton_ESpec_Lab[1][i], "Ekin/A(MeV)", "d^{2}M/d#OmegadE(Sr^{-1}MeV^{-1})");
  
    h1_Proton_ESpec_Lab[0][i]->Rebin(5);
    h1_Proton_ESpec_Lab[0][i]->Scale(1.0/5.0);
    h1_Proton_ESpec_Lab[1][i]->Rebin(5);
    h1_Proton_ESpec_Lab[1][i]->Scale(1.0/5.0);
    
    h1_Triton_ESpec_Lab[0][i]->Rebin(5);
    h1_Triton_ESpec_Lab[0][i]->Scale(1.0/5.0);
    h1_Triton_ESpec_Lab[1][i]->Rebin(5);
    h1_Triton_ESpec_Lab[1][i]->Scale(1.0/5.0);
  
  }
  sprintf(NameTem,"c1_ESpecLab_Proton_%s",SystemTag_1.c_str());
  TCanvas* c1_ESpecLab_Proton_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_ESpecLab_Proton_Sys1->cd();
  for(int i=0;i<AngleNum;i++)
  {
    h1_Proton_ESpec_Lab[0][i]->SetLineColor(AngleLineColor[i]);
    
    if(i==0) 
    { 
      h1_Proton_ESpec_Lab[0][i]->Draw("hist");
      h1_Proton_ESpec_Lab[0][i]->GetXaxis()->SetRangeUser(30,200);
    }
    else {h1_Proton_ESpec_Lab[0][i]->Draw("histsame");}
  }
  
  int AngleIndex_ForCompare = 3;
  sprintf(NameTem,"c1_ESpecLab_Proton_%s_%s",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_ESpecLab_Proton_Sys1_Sys2 = new TCanvas(NameTem,NameTem,1);
  //h1_Proton_ESpec_Lab[0][1]->Divide(h1_Proton_ESpec_Lab[1][1]);
  h1_Proton_ESpec_Lab[0][AngleIndex_ForCompare ]->Draw("hist");
  h1_Proton_ESpec_Lab[1][AngleIndex_ForCompare ]->Draw("histsame");
  
  sprintf(NameTem,"c1_ESpecLab_Triton_%s_%s",SystemTag_1.c_str(),SystemTag_2.c_str());
  TCanvas* c1_ESpecLab_Triton_Sys1_Sys2 = new TCanvas(NameTem,NameTem,1);
  h1_Triton_ESpec_Lab[0][AngleIndex_ForCompare ]->Draw("hist");
  h1_Triton_ESpec_Lab[1][AngleIndex_ForCompare ]->Draw("histsame");
  
}

TH1D* Select(TH2D* h2_Ekin_Theta_Lab, double Theta, string HistoName, double DeltaTheta)
{
  TAxis* Y_Axis = h2_Ekin_Theta_Lab->GetYaxis();
  double Theta_1 = Theta-0.5*DeltaTheta;
  double Theta_2 = Theta+0.5*DeltaTheta;
  int Bin1 = Y_Axis->FindBin(Theta_1);
  int Bin2 = Y_Axis->FindBin(Theta_2);
  
  TH1D* h1_ESpec = (TH1D*) h2_Ekin_Theta_Lab->ProjectionX(HistoName.c_str(),Bin1,Bin2);
  double SolidAngle = 2*Pi()*DeltaTheta*DegToRad()*Sin(Theta*DegToRad());
  h1_ESpec->Scale(1.0/SolidAngle);
return h1_ESpec;
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
  
  h1->SetLineWidth(3);
  h1->SetMarkerSize(1.5);
}
