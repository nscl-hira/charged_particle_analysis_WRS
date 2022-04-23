#include "TMath.h"
using namespace TMath;

//this script will draw all the different particle together.
void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle);
void Set_2D_HistoStyle(TH2D* h2, string XTitle, string YTitle);
void Set_1D_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle);
TH1D* GetTemperature(string NameTem,TH1D** h1_PtA_Lab);
TH2D* Rebin_ESpec(TH2D* h2_y_PtA_Lab, string NameTem, double* RapidityRange,double Rapidity_BinSize,double* PtARange,double PtA_BinSize);
void RemovePoint_LargeRelError(TH1D* h1,double RelErr);
void Cal_Alpha_Beta(int ParticleNum, double* Z, double* N, double* R21, double* R21_Err, double* FitResults);
double GetYield(TH2D* h2_y_PtA_Lab,double* RapidityRange, double* PtARange,double* Err);

void _3_Cal_IsoScaling_1Region()
{
  gStyle->SetErrorX(0);

  string SystemTag_1 = "Ca40Ni58E140"; 
  string SystemTag_2 = "Ca48Ni64E140";
  
  string ThetaTag = "Noise_Theta_30_75";
  string EffCor_Tag = "GeoEff";
  string ImpactParaTag = "bHat_0.00_0.40";
  
  double RelErr_Cut = 0.05;
  
  double RapidityRange[2] = {0.4,0.6};
  double PtARange[2] = {180,200.0};
  
  const int ParticleNum = 5;
  string ParticleName[]={"p","d","t","3He","4He","6He","6Li","7Li"};
  int LineColor[] = {1,1,1,2,2,2,4,4};
  int MarkerStyle[] = {20,24,21,20,24,21,20,24};
  
  TH2D* h2_y_PtA_Lab[2][ParticleNum];
  double Yield_1Region[2][ParticleNum];
  double Yield_Err_1Region[2][ParticleNum];
  double R21_1Region[ParticleNum];
  double R21_Err_1Region[ParticleNum];
  
  
  char NameTem[200];
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_1.c_str(),ImpactParaTag.c_str());
  TFile* f1_MergedESpec_1 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_2.c_str(),ImpactParaTag.c_str());
  TFile* f1_MergedESpec_2 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;

  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    cout<<"Read : "<<NameTem<<endl;
    f1_MergedESpec_1->cd();
    h2_y_PtA_Lab[0][iP] = (TH2D*) f1_MergedESpec_1->Get(NameTem);
    sprintf(NameTem,"h2_y_PtA_Lab_%s_%s_%s_%s",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_y_PtA_Lab[0][iP]->SetName(NameTem);
    double Err1 = 0; 
    Yield_1Region[0][iP] = GetYield(h2_y_PtA_Lab[0][iP],RapidityRange,PtARange,&Err1);
    Yield_Err_1Region[0][iP] = Err1;
    
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    f1_MergedESpec_2->cd();
    h2_y_PtA_Lab[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem);
    sprintf(NameTem,"h2_y_PtA_Lab_%s_%s_%s_%s",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_y_PtA_Lab[1][iP]->SetName(NameTem);
    double Err2 = 0;
    Yield_1Region[1][iP] = GetYield(h2_y_PtA_Lab[1][iP],RapidityRange,PtARange,&Err2);
    Yield_Err_1Region[1][iP] = Err2;
    
    R21_1Region[iP] = Yield_1Region[1][iP]/Yield_1Region[0][iP];
    double V1 = Yield_1Region[0][iP]; double V2 = Yield_1Region[1][iP];
    R21_Err_1Region[iP] = R21_1Region[iP]*Sqrt(Power(Err1/V1,2)+Power(Err2/V2,2));
  }
  
  //Draw the Pt/A~Rapidity
  sprintf(NameTem,"c1_Pt_A_Rapidity_%s_%s_%s",SystemTag_1.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Sys1->Divide(3,2);
  
  sprintf(NameTem,"c1_Pt_A_Rapidity_%s_%s_%s",SystemTag_2.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Sys2 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Sys2->Divide(3,2);
  
  for(int i=0;i<ParticleNum;i++)
  {
    c1_Pt_A_Rapidity_Sys1->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab[0][i]->Draw("colz");
    
    c1_Pt_A_Rapidity_Sys2->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab[1][i]->Draw("colz");
  }

  //the below is for calculating the alpha and Beta for 1 region.
  double Z[] = {1,1,1,2,2,2,3,3}; double N[] = {0,1,2,1,2,4,3,4};
  int IsValue = 1;
  for(int iP=0;iP<ParticleNum;iP++) 
  { 
    cout<<"R21_"<<ParticleName[iP]<<" : "<<R21_1Region[iP]<<endl;
    if(R21_1Region[iP]==0) { IsValue=0; } 
  }
  double Results[4] = {0,0,0,0}; //[Alpha, Alpha_Err, Beta, Beta_Err]

  if(IsValue==1)
  {
    Cal_Alpha_Beta(ParticleNum, Z, N, R21_1Region, R21_Err_1Region, Results);
  }

}

double GetYield(TH2D* h2_y_PtA_Lab,double* RapidityRange, double* PtARange,double* Err)
{
  TAxis* X_Axis = h2_y_PtA_Lab->GetXaxis();
  TAxis* Y_Axis = h2_y_PtA_Lab->GetYaxis();
  
  int XBinRange[2] = {X_Axis->FindBin(RapidityRange[0]),X_Axis->FindBin(RapidityRange[1])};
  int YBinRange[2] = {Y_Axis->FindBin(PtARange[0]),Y_Axis->FindBin(PtARange[1])};
  
  double Yield = 0;
  double Error = 0;
  
  for(int iX=XBinRange[0];iX<=XBinRange[1];iX++)
  {
    for(int iY=YBinRange[0];iY<=YBinRange[1];iY++)
    {
      double BinContent = h2_y_PtA_Lab->GetBinContent(iX,iY);
      double BinErr = h2_y_PtA_Lab->GetBinError(iX,iY);
      Yield += BinContent;
      if(BinContent>0) { Error += Power(BinErr,2); }
    }
  }
  
  if(Yield>0) { *Err = Sqrt(Error); }
return Yield;
}

TH2D* Rebin_ESpec(TH2D* h2_y_PtA_Lab, string NameTem, double* RapidityRange,double Rapidity_BinSize,double* PtARange,double PtA_BinSize)
{
  int XBinNum_Rebin = (RapidityRange[1]-RapidityRange[0])/Rapidity_BinSize;
  int YBinNum_Rebin = (PtARange[1]-PtARange[0])/PtA_BinSize;
  cout<<XBinNum_Rebin<<"  "<<YBinNum_Rebin<<endl;
  
  TH2D* h2_y_PtA_Lab_Rebin = new TH2D(NameTem.c_str(),";y/y_{BeamLab};Pt/A(MeV/c)",XBinNum_Rebin,RapidityRange[0],RapidityRange[1],YBinNum_Rebin,PtARange[0],PtARange[1]);
  TH2D* h2_y_PtA_Lab_Rebin_Err = (TH2D*) h2_y_PtA_Lab_Rebin->Clone((NameTem+"_Err").c_str());
  
  TAxis* X_Axis = h2_y_PtA_Lab->GetXaxis();
  TAxis* Y_Axis = h2_y_PtA_Lab->GetYaxis();
  int XBinNum = X_Axis->GetNbins();
  int YBinNum = Y_Axis->GetNbins();
  double XBinSize = X_Axis->GetBinCenter(2)-X_Axis->GetBinCenter(1);
  double YBinSize = Y_Axis->GetBinCenter(2)-Y_Axis->GetBinCenter(1);
  
  for(int iX=1;iX<=XBinNum;iX++)
  {
    for(int iY=1;iY<=YBinNum;iY++)
    {
      double Rapidity = X_Axis->GetBinCenter(iX)+XBinSize*gRandom->Uniform(-0.5,0.5);
      double PtA = Y_Axis->GetBinCenter(iY)+YBinSize*gRandom->Uniform(-0.5,0.5);
      double BinContent = h2_y_PtA_Lab->GetBinContent(iX,iY);
      double BinContentErr = h2_y_PtA_Lab->GetBinError(iX,iY);
      
      h2_y_PtA_Lab_Rebin->Fill(Rapidity,PtA,BinContent);
      h2_y_PtA_Lab_Rebin_Err->Fill(Rapidity,PtA,BinContentErr*BinContentErr);
    }
  }
  
  for(int iX=1;iX<=XBinNum;iX++)
  {
    for(int iY=1;iY<=YBinNum;iY++)
    {
      double BinErr = h2_y_PtA_Lab_Rebin_Err->GetBinContent(iX,iY);
      BinErr = Sqrt(BinErr);
      h2_y_PtA_Lab_Rebin->SetBinError(iX,iY,BinErr);
      double BinContent = h2_y_PtA_Lab_Rebin->GetBinContent(iX,iY);
      if(BinErr/BinContent>0.1) { h2_y_PtA_Lab_Rebin->SetBinContent(iX,iY,0); }
    }
  }
  
return h2_y_PtA_Lab_Rebin;
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

void Cal_Alpha_Beta(int ParticleNum, double* Z, double* N, double* R21, double* R21_Err, double* FitResults)
{
  double Z_Err[20] = {0}; double N_Err[20] = {0};
  double Ln_R21[20] = {0}; double Ln_R21_Err[20] = {0};
  for(int i=0;i<ParticleNum;i++)
  { 
    Z_Err[i] = 0; N_Err[i] = 0;
    Ln_R21[i] = Log(R21[i]); 
    Ln_R21_Err[i] = R21_Err[i]/R21[i];
  }
  int Index = 1;
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
  TF2* f2 = new TF2(NameTem,"[0]*x+[1]*y",-1,3,1,3);
  f2->SetParameters(0.3,-0.3);
  gr2->Fit(f2);
  TF2 *fit2 = (TF2*)gr2->FindObject(NameTem);
  double Results[2]; double ResultsErr[2];
  f2->GetParameters(Results);
  ResultsErr[0] = f2->GetParError(0); ResultsErr[1] = f2->GetParError(1);
  
  FitResults[0] = Results[0]; FitResults[1] = ResultsErr[0];
  FitResults[2] = Results[1]; FitResults[1] = ResultsErr[1];
  
  //the below is only for drawing the R21.
  TCanvas* c1_R21_N_Z = new TCanvas("c1_R21_N_Z","c1_R21_N_Z",1);
  c1_R21_N_Z->Divide(2,1);
  TGraphErrors* gr1_R21_N = new TGraphErrors(ParticleNum,N,R21,N_Err,R21_Err);
  TGraphErrors* gr1_R21_Z = new TGraphErrors(ParticleNum,Z,R21,Z_Err,R21_Err);
  c1_R21_N_Z->cd(1)->SetLogy(1);
  gr1_R21_N->Draw("A*");
  gr1_R21_N->GetXaxis()->SetLimits(0,6);
  gr1_R21_N->GetYaxis()->SetRangeUser(0.5,2.5);
  gr1_R21_N->GetXaxis()->SetTitle("N");
  gr1_R21_N->GetXaxis()->CenterTitle(1);
  gr1_R21_N->GetYaxis()->SetTitle("R21");
  gr1_R21_N->GetYaxis()->SetTitleOffset(1);
  gr1_R21_N->GetYaxis()->CenterTitle(1);
  
  
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d)",Results[0],Results[1],1);
  TF1* f1_R21_Z_H = new TF1("f1_R21_Z_H",NameTem,-0.2,2.2);
  f1_R21_Z_H->Draw("same");
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d)",Results[0],Results[1],2);
  TF1* f1_R21_Z_He = new TF1("f1_R21_Z_He",NameTem,0.8,4.2);
  f1_R21_Z_He->Draw("same");
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d)",Results[0],Results[1],3);
  TF1* f1_R21_Z_Li = new TF1("f1_R21_Z_Li",NameTem,2.8,5.2);
  f1_R21_Z_Li->Draw("same"); f1_R21_Z_Li->SetLineStyle(7);
  sprintf(NameTem,"TMath::Exp(x*%.3f+%.3f*%d)",Results[0],Results[1],4);
  TF1* f1_R21_Z_Be = new TF1("f1_R21_Z_Be",NameTem,2.8,6.2);
  f1_R21_Z_Be->Draw("same"); f1_R21_Z_Be->SetLineStyle(7);
  
  c1_R21_N_Z->cd(2)->SetLogy(1);
  gr1_R21_Z->Draw("A*");
  gr1_R21_Z->GetXaxis()->SetLimits(0,6);
  gr1_R21_Z->GetYaxis()->SetRangeUser(0.5,2.5);
  gr1_R21_Z->GetXaxis()->SetTitle("Z");
  gr1_R21_Z->GetXaxis()->CenterTitle(1);
  gr1_R21_Z->GetYaxis()->SetTitle("R21");
  gr1_R21_Z->GetYaxis()->CenterTitle(1);
  gr1_R21_Z->GetYaxis()->SetTitleOffset(1);
  
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x)",1,Results[0],Results[1]);
  TF1* f1_R21_N1_Z1 = new TF1("f1_R21_N1_Z1",NameTem,0.8,2.2);
  f1_R21_N1_Z1->Draw("same");
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x)",2,Results[0],Results[1]);
  TF1* f1_R21_N2_Z1 = new TF1("f1_R21_N2_Z1",NameTem,0.8,2.2);
  f1_R21_N2_Z1->Draw("same");
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x)",4,Results[0],Results[1]);
  TF1* f1_R21_N4_Z2 = new TF1("f1_R21_N4_Z2",NameTem,1.8,3.2);
  f1_R21_N4_Z2->Draw("same"); f1_R21_N4_Z2->SetLineStyle(7);
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x)",3,Results[0],Results[1]);
  TF1* f1_R21_N3_Z3 = new TF1("f1_R21_N3_Z3",NameTem,2.8,4.2);
  f1_R21_N3_Z3->Draw("same"); f1_R21_N3_Z3->SetLineStyle(7);
  sprintf(NameTem,"TMath::Exp(%.d*%.3f+%.3f*x)",5,Results[0],Results[1]);
  TF1* f1_R21_N5_Z3 = new TF1("f1_R21_N5_Z3",NameTem,2.8,4.2);
  f1_R21_N5_Z3->Draw("same"); f1_R21_N5_Z3->SetLineStyle(7);
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
  
  h1->SetLineWidth(2);
  h1->SetMarkerSize(1.5);
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
