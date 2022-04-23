#include "TMath.h"
using namespace TMath;

//this script will draw all the different particle together.
void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle, int MarkerStyle, int MarkerColor);
void Set_1D_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle);
TH1D* Select_Ekin_A(string SystemName, string ParticleName, TH2D* h2_Ekin_Theta_Lab,double* ThetaRange, double* EkinCut);
void RemovePoint_LargeRelError(TH1D* h1,double RelErr);

void Cal_Alpha_Beta(int Index, double* Ekin_BinEdge, TCanvas* c1_N, TCanvas* c1_Z, int ParticleNum, double* Z, double* N, double* R21, double* R21_Err,double* Alpha, double* Alpha_Err, double* Beta, double* Beta_Err, double* C, double* C_Err);

void _2_Cal_Isoscaling_Lab()
{
  gStyle->SetErrorX(0);
  
  string SystemTag_1 = "Ca40Ni58E140"; 
  string SystemTag_2 = "Ca48Ni64E140"; 
  
  double RelErr_Cut = 0.05;

  int RebinNum = 10;
  bool Is_R21Err_Applied=1;
  
  int iEff = 1;
  int ib = 0;
  double ThetaRange[2] = {45,55};
  
  const int EffOpt = 3;
  string EffCor_Tag[] = {"_noEff", "_GeoEff",  "_GeoReactionEff"};
  
  const int ParticleNum = 6;
  string ParticleName[]={"p","d","t","3He","4He","6He","6Li","7Li"};
  int LineColor[] = {1,1,1,2,2,2};
//  int MarkerStyle[] = {1,1,1,1,1,1};
  int MarkerStyle[] = {21,25,36,22,26,33};
  
  double EkinCut[ParticleNum][2] = {{20,198.0/1.0},{15,263.0/2.0},{12,312.0/3.0},{20,200},{18,200},{18,200}};
  //double EkinCut[ParticleNum][2] = {{20,100},{20,100},{20,100},{20,100},{20,100},{20,100}};
  
  const int ImpactParaNum = 3;
  string ImpactParaTag[] = {"bHat_0.00_0.40","bHat_0.40_0.70","bHat_0.70_1.10"};
  string ImpactPara_Select_Tag = "bHat";
  
  char NameTem[200];
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_1.c_str(),ImpactParaTag[ib].c_str());
  TFile* f1_MergedESpec_1 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;
  sprintf(NameTem,"./%s/f1_MergedData_%s.root",SystemTag_2.c_str(),ImpactParaTag[ib].c_str());
  TFile* f1_MergedESpec_2 = new TFile(NameTem,"read");
  cout<<"Open "<<NameTem<<endl;

  TH2D* h2_Ekin_Theta_Lab[2][ParticleNum];
  TH1D* h1_Ekin_Lab[2][ParticleNum];
  TH1D* h1_R21[ParticleNum];
  
  double Integral_Yield[2][ParticleNum];
  double Integral_YieldErr[2][ParticleNum];
  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    cout<<"Read : "<<NameTem<<endl;
    f1_MergedESpec_1->cd();
    h2_Ekin_Theta_Lab[0][iP] = (TH2D*) f1_MergedESpec_1->Get(NameTem);
    sprintf(NameTem,"h2_%s_%s_%s_Ekin_Theta_Lab%s_bhat",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
    h2_Ekin_Theta_Lab[0][iP]->SetName(NameTem);
    h1_Ekin_Lab[0][iP] = Select_Ekin_A(SystemTag_1, ParticleName[iP], h2_Ekin_Theta_Lab[0][iP], ThetaRange, EkinCut[iP]);
    Set_1D_HistoStyle(h1_Ekin_Lab[0][iP],"Ekin/A(MeV)","d^{2}M/dEd#Omega(MeV^{-1} #bullet Sr^{-1})",MarkerStyle[iP],2);
    
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    f1_MergedESpec_2->cd();
    h2_Ekin_Theta_Lab[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem);
    sprintf(NameTem,"h2_%s_%s_%s_Ekin_Theta_Lab%s_bhat",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
    h2_Ekin_Theta_Lab[1][iP]->SetName(NameTem);
    h1_Ekin_Lab[1][iP] = Select_Ekin_A(SystemTag_2, ParticleName[iP], h2_Ekin_Theta_Lab[1][iP], ThetaRange, EkinCut[iP]);
    Set_1D_HistoStyle(h1_Ekin_Lab[1][iP],"Ekin/A(MeV)","d^{2}M/dEd#Omega(MeV^{-1} #bullet Sr^{-1})",MarkerStyle[iP],1);
    
    //for calculate the integral yield
    int Bin1 = h1_Ekin_Lab[0][iP]->FindBin(EkinCut[iP][0]);
    int Bin2 = h1_Ekin_Lab[0][iP]->FindBin(EkinCut[iP][1]);
    double Err1 = 0; double Err2 = 0;
    Integral_Yield[0][iP] = h1_Ekin_Lab[0][iP]->IntegralAndError(Bin1,Bin2-1,Err1,"");
    Integral_YieldErr[0][iP] = Err1;
    
    Integral_Yield[1][iP] = h1_Ekin_Lab[1][iP]->IntegralAndError(Bin1,Bin2-1,Err2,"");
    Integral_YieldErr[1][iP] = Err2;
    
    //Rebin the energy spectrum 
    h1_Ekin_Lab[0][iP]->Rebin(RebinNum);
    h1_Ekin_Lab[0][iP]->Scale(1.0/RebinNum);
    
    h1_Ekin_Lab[1][iP]->Rebin(RebinNum);
    h1_Ekin_Lab[1][iP]->Scale(1.0/RebinNum);
  }
  
  //Draw the Ekin~Theta
  sprintf(NameTem,"c1_Ekin_Theta_Lab_%s_%s_%s",SystemTag_1.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
  TCanvas* c1_Ekin_Theta_Lab_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_Ekin_Theta_Lab_Sys1->Divide(3,2,0,0);
  
  sprintf(NameTem,"c1_Ekin_Theta_Lab_%s_%s_%s",SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str());
  TCanvas* c1_Ekin_Theta_Lab_Sys2 = new TCanvas(NameTem,NameTem,1);
  c1_Ekin_Theta_Lab_Sys2->Divide(3,2,0,0);
  for(int i=0;i<ParticleNum;i++)
  {
    c1_Ekin_Theta_Lab_Sys1->cd(i+1);//->SetGrid(1,1);
    h2_Ekin_Theta_Lab[0][i]->Draw("col");
    c1_Ekin_Theta_Lab_Sys2->cd(i+1);//->SetGrid(1,1);
    h2_Ekin_Theta_Lab[1][i]->Draw("col");
  }
  
  //Draw the Ekin spectrum
  sprintf(NameTem,"c1_Ekin_Lab_%s_%s_%s_%s_Theta_%.1f_%.1f", SystemTag_1.c_str(), SystemTag_2.c_str(), ImpactParaTag[ib].c_str(), EffCor_Tag[iEff].c_str(), ThetaRange[0], ThetaRange[1]);
  TCanvas* c1_Ekin_Lab = new TCanvas(NameTem,NameTem,1);
  c1_Ekin_Lab->Divide(3,2);
  
  for(int i=0;i<ParticleNum;i++)
  {
    c1_Ekin_Lab->cd(i+1)->SetLogy(1);
    h1_Ekin_Lab[0][i]->Draw();
    h1_Ekin_Lab[1][i]->Draw("same");
    if(i<3)
    { 
      h1_Ekin_Lab[0][i]->GetXaxis()->SetRangeUser(EkinCut[i][0],EkinCut[i][1]);
      h1_Ekin_Lab[0][i]->GetYaxis()->SetRangeUser(2E-4,4E-2);
    }
    else
    {
      h1_Ekin_Lab[0][i]->GetXaxis()->SetRangeUser(EkinCut[i][0],EkinCut[i][1]);
      h1_Ekin_Lab[0][i]->GetYaxis()->SetRangeUser(2E-6,2E-2);
    }
  }
  
   //Draw R21
  sprintf(NameTem,"c1_R21_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
  TCanvas* c1_R21 = new TCanvas(NameTem,NameTem,1);
  TLegend* Legend_R21 = new TLegend(0.6,0.6,0.8,0.8);
  for(int i=0;i<ParticleNum;i++)
  {
    c1_R21->cd()->SetGrid(1,1);
    sprintf(NameTem,"h1_R21_%s_%s_%s_%s_%s_y_%.1f_%.1f",ParticleName[i].c_str(),SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
    h1_R21[i] = (TH1D*) h1_Ekin_Lab[1][i]->Clone(NameTem);
    h1_R21[i]->Divide(h1_Ekin_Lab[0][i]);
    h1_R21[i]->GetYaxis()->SetTitle("R21");
    h1_R21[i]->GetXaxis()->SetTitle("Ekin/A(MeV)");
    
    RemovePoint_LargeRelError(h1_R21[i],RelErr_Cut);
    if(i==0)
    {
      h1_R21[i]->Draw();
      h1_R21[i]->GetYaxis()->SetTitle("R21");
      h1_R21[i]->GetYaxis()->SetRangeUser(0.6,2);
    }
    else if(i!=5) { h1_R21[i]->Draw("same"); }
    
    if(i<3) { h1_R21[i]->SetMarkerColor(1); h1_R21[i]->SetLineColor(1); }
    else { h1_R21[i]->SetMarkerColor(2); h1_R21[i]->SetLineColor(2); }
    
    if(i!=5) { Legend_R21->AddEntry(h1_R21[i],ParticleName[i].c_str(),"p"); }
  }
  Legend_R21->Draw();  
  
  //The below is for calculating the Alpha and Beta.
  double Z[] = {1,1,1,2,2,2}; double N[] = {0,1,2,1,2,4};
  int Ekin_Num = 0;
  double Ekin[100] = {0}; double Ekin_Err[100] = {0}; 
  double Alpha[100] = {0}; double Alpha_Err[100] = {0}; double Beta[100] = {0}; double Beta_Err[100] = {0};
  double C[100] = {0}; double C_Err[100] = {0};
  double Alpha_Plus_Beta[100] = {0}; double Alpha_Plus_Beta_Err[100] = {0};
  double Alpha_Subtract_Beta[100] = {0}; double Alpha_Subtract_Beta_Err[100] = {0};
  
  sprintf(NameTem,"c1_N_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
  TCanvas* c1_N = new TCanvas(NameTem,NameTem,1);
  c1_N->Divide(3,3);
  sprintf(NameTem,"c1_Z_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
  TCanvas* c1_Z = new TCanvas(NameTem,NameTem,1);
  c1_Z->Divide(3,3);
  
  //int FirstBinIndex = h1_R21[0]->FindFirstBinAbove(0);
  int XBinNum = h1_R21[0]->GetNbinsX();
  for(int i=1;i<=XBinNum;i++) 
  {
    int IsValue = 1; //if each particle have R21, then the alpha and beta will be calculated.
    
    double Ekin_R21[6] = {0}; double Ekin_R21_Err[6] = {0};
    for(int iP=0;iP<ParticleNum-1;iP++)
    {
      Ekin_R21[iP] = h1_R21[iP]->GetBinContent(i);
      if(Is_R21Err_Applied==1) { Ekin_R21_Err[iP] = h1_R21[iP]->GetBinError(i); }
      else { Ekin_R21_Err[iP] = 0.02*Ekin_R21[iP]; }
      if(Ekin_R21[iP]==0) { IsValue=0; }
    }
    
    if(IsValue==1) 
    {
      Ekin[Ekin_Num] = h1_R21[0]->GetBinCenter(i);
      cout<<endl<<endl<<" --> Ekin/A(MeV/c) : "<<Ekin[Ekin_Num]<<endl;
      cout.width(7);
      cout<<"      "<<"R21"<<"  "<<"R21_Err"<<endl;
      for(int iP=0;iP<ParticleNum-1;iP++)
      {
        cout.width(7);
        cout<<ParticleName[iP]<<"  "<<Ekin_R21[iP]<<"  "<<Ekin_R21_Err[iP]<<endl;
      }
      double Ekin_BinWidth = h1_R21[0]->GetBinWidth(i);
      double Ekin_BinEdge[2] = {Ekin[Ekin_Num]-0.5*Ekin_BinWidth,Ekin[Ekin_Num]+0.5*Ekin_BinWidth};
      
      Cal_Alpha_Beta(Ekin_Num, Ekin_BinEdge, c1_N, c1_Z, ParticleNum-1, Z, N, Ekin_R21, Ekin_R21_Err, Alpha, Alpha_Err, Beta, Beta_Err, C, C_Err);
      Alpha_Plus_Beta[Ekin_Num] = Alpha[Ekin_Num]+Beta[Ekin_Num];
      Alpha_Plus_Beta_Err[Ekin_Num] = Sqrt(Power(Alpha_Err[Ekin_Num],2)+Power(Beta_Err[Ekin_Num],2));
      Alpha_Subtract_Beta[Ekin_Num] = Alpha[Ekin_Num]-Beta[Ekin_Num];
      Alpha_Subtract_Beta_Err[Ekin_Num] = Sqrt(Power(Alpha_Err[Ekin_Num],2)+Power(Beta_Err[Ekin_Num],2));
      Ekin_Num++; 
    }
  }
  
  cout<<"Ekin Num : "<<Ekin_Num<<endl;
  TGraphErrors* gr1_Alpha = new TGraphErrors(Ekin_Num,Ekin,Alpha,Ekin_Err,Alpha_Err);
  TGraphErrors* gr1_Beta = new TGraphErrors(Ekin_Num,Ekin,Beta,Ekin_Err,Beta_Err);
  TGraphErrors* gr1_Alpha_Plus_Beta = new TGraphErrors(Ekin_Num,Ekin,Alpha_Plus_Beta,Ekin_Err,Alpha_Plus_Beta_Err);
  TGraphErrors* gr1_Alpha_Subtract_Beta = new TGraphErrors(Ekin_Num,Ekin,Alpha_Subtract_Beta,Ekin_Err,Alpha_Subtract_Beta_Err);
  
  TGraphErrors* gr1_C_Fit = new TGraphErrors(Ekin_Num,Ekin,C,Ekin_Err,C_Err);
  
  sprintf(NameTem,"c1_Alpha_Beta_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
  TCanvas* c1_Alpha_Beta = new TCanvas(NameTem,NameTem,1);
  
  gr1_Alpha->Draw("A*"); gr1_Alpha->SetMarkerStyle(20); gr1_Alpha->SetMarkerColor(2);
  gr1_Beta->Draw("same*"); gr1_Beta->SetMarkerStyle(21); gr1_Beta->SetMarkerColor(4);
  gr1_Alpha_Plus_Beta->Draw("same*"); gr1_Alpha_Plus_Beta->SetMarkerStyle(22); gr1_Alpha_Plus_Beta->SetMarkerColor(6);
  
  gr1_Alpha->GetYaxis()->SetRangeUser(-0.4,0.4);
  Set_1D_GraphStyle(gr1_Alpha,"Ekin/A(MeV)", "#Delta #mu/T");
  Set_1D_GraphStyle(gr1_Beta,"Ekin/A(MeV)", "#Delta #mu/T");
  gr1_Alpha->SetMarkerStyle(20); gr1_Alpha->SetMarkerColor(1); 
  gr1_Beta->SetMarkerStyle(24); gr1_Beta->SetMarkerColor(1);
  
  TLegend* legend_Alpha_Beta = new TLegend(0.6,0.6,0.8,0.8);
  legend_Alpha_Beta->SetNColumns(3);
  legend_Alpha_Beta->AddEntry(gr1_Alpha,"#Delta#mu_{n}/T","p");
  legend_Alpha_Beta->AddEntry(gr1_Beta,"#Delta#mu_{p}/T","p");
  legend_Alpha_Beta->AddEntry(gr1_Alpha_Plus_Beta,"(#Delta#mu_{n}+#Delta#mu_{p})/T","p");
  legend_Alpha_Beta->Draw();
  
  //the below is for calculating the R21 from the alpha and beta.
  c1_R21->cd();
  TGraph* gr1_R21_Alpha_Beta[ParticleNum] = {0};
  for(int iP=0;iP<ParticleNum-1;iP++)
  {
    double R21_Particle[100] = {0};
    for(int i=0;i<Ekin_Num;i++)
    {
      R21_Particle[i] = Exp((Alpha[i]*N[iP]+Beta[i]*Z[iP]+C[i]));
    }
    gr1_R21_Alpha_Beta[iP] = new TGraph(Ekin_Num,Ekin,R21_Particle);
    gr1_R21_Alpha_Beta[iP]->Draw("sameL");
    gr1_R21_Alpha_Beta[iP]->SetLineColor(LineColor[iP]);
    gr1_R21_Alpha_Beta[iP]->SetLineWidth(2);
    gr1_R21_Alpha_Beta[iP]->SetLineStyle(2);
  }
  
  sprintf(NameTem,"c1_C_Fit_%s_%s_%s_%s_y_%.1f_%.1f",SystemTag_1.c_str(), SystemTag_2.c_str(), ImpactParaTag[ib].c_str(), EffCor_Tag[iEff].c_str(), ThetaRange[0], ThetaRange[1]);
  TCanvas* c1_C_Fit = new TCanvas(NameTem,NameTem,1);
  gr1_C_Fit->Draw("AL*");
 
  //The R21 of integral yield
  double Integral_R21[ParticleNum];
  double Ekin_BinEdge[2] = {0,200};
  double Integral_R21_Err[ParticleNum];
  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    Integral_R21[iP] = Integral_Yield[1][iP]/Integral_Yield[0][iP];
    Integral_R21_Err[iP] = Integral_R21[iP]*Sqrt(Power(Integral_YieldErr[1][iP]/Integral_Yield[1][iP],2)+Power(Integral_YieldErr[0][iP]/Integral_Yield[0][iP],2));
    
    cout<<ParticleName[iP]<<" : "<< Integral_Yield[0][iP]<<" +/- "<<Integral_YieldErr[0][iP]
                          <<"   "<<Integral_Yield[1][iP]<<" +/- "<<Integral_YieldErr[1][iP]
                          <<" R21 = "<<Integral_R21[iP]<<" +/- "<<Integral_R21_Err[iP]<<endl;
  }
  
  sprintf(NameTem,"c1_Integral_N_%s_%s_%s_%s_y_%.1f_%.1f", SystemTag_1.c_str(), SystemTag_2.c_str(), ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
  TCanvas* c1_Integral_N = new TCanvas(NameTem,NameTem,1);
  c1_Integral_N->Divide(1,1);
  
  sprintf(NameTem,"c1_Integral_Z_%s_%s_%s_%s_y_%.1f_%.1f", SystemTag_1.c_str(), SystemTag_2.c_str(), ImpactParaTag[ib].c_str(),EffCor_Tag[iEff].c_str(),ThetaRange[0],ThetaRange[1]);
  TCanvas* c1_Integral_Z = new TCanvas(NameTem,NameTem,1);
  c1_Integral_Z->Divide(1,1);
  Cal_Alpha_Beta(0, Ekin_BinEdge, c1_Integral_N, c1_Integral_Z, ParticleNum-1, Z, N, Integral_R21, Integral_R21_Err, Alpha, Alpha_Err, Beta, Beta_Err, C, C_Err);

/*
  //the R21 of sub-region
  double EnergyRange[2] = {}; int EnergyNum = ;
  double ThetaRange_SubRegion[2] = {}; int ThetaNum = ;
  
  for(int iT=0;iT<ThetaNum;iT++)
  {
    double Theta1 = ; double Theta2 = ;
    for(int iE=0;iE<EnergyNum;iE++)
    {
      double E1 = ; double E2 = ;
      
      for(int iP=0;iP<ParticleNum;iP++)
      {
        
      }
    }
  }
 */
}

TH1D* Select_Ekin_A(string SystemName, string ParticleName, TH2D* h2_Ekin_Theta_Lab,double* ThetaRange,double* EkinCut)
{
  int YBin1 = h2_Ekin_Theta_Lab->GetYaxis()->FindBin(ThetaRange[0]);
  int YBin2 = h2_Ekin_Theta_Lab->GetYaxis()->FindBin(ThetaRange[1]);
  
  char NameTem[200];
  sprintf(NameTem,"h1_Ekin_Lab_%s_%s_%.1f_%.1f",SystemName.c_str(), ParticleName.c_str(), ThetaRange[0],ThetaRange[1]);
  
  TH1D* h1_Ekin = (TH1D*) h2_Ekin_Theta_Lab->ProjectionX(NameTem,YBin1,YBin2);
  double SolidAngle = TwoPi()*(Cos(DegToRad()*ThetaRange[0])-Cos(DegToRad()*ThetaRange[1]));
  h1_Ekin->Scale(1.0/SolidAngle); //the rapidity normalization.
  Set_1D_HistoStyle(h1_Ekin,"Ekin/A(MeV)","d^{2}M/d_{E}d_{#Omega}(MeV^{-1} #bullet Sr^{-1})",20,1);
  cout<<" Produce " << h1_Ekin->GetName()<<endl;
  
  //make a energy cut
  int XBinNum = h1_Ekin->GetNbinsX();
  for(int i=1;i<=XBinNum;i++)
  {
    double BinCenter = h1_Ekin->GetBinCenter(i);
    if(BinCenter<EkinCut[0] || BinCenter>EkinCut[1])
    {
      h1_Ekin->SetBinContent(i,0);
      h1_Ekin->SetBinError(i,0);
    }
  }
  
return h1_Ekin;
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

void Cal_Alpha_Beta(int Index, double* Ekin_BinEdge, TCanvas* c1_N, TCanvas* c1_Z, int ParticleNum, double* Z, double* N, double* R21, double* R21_Err,double* Alpha, double* Alpha_Err, double* Beta, double* Beta_Err, double* C, double* C_Err)
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
  
  TGraphErrors* gr1_R21_N_H = new TGraphErrors(3,N,R21,N_Err,R21_Err);
  TGraphErrors* gr1_R21_Z_H = new TGraphErrors(3,Z,R21,Z_Err,R21_Err);
  
  double N_He[3]; double N_Err_He[3]; 
  double Z_He[3]; double Z_Err_He[3];
  double R21_He[3]; double R21_Err_He[3];
  
  for(int i=0;i<2;i++)
  {
    N_He[i] = N[3+i]; N_Err_He[i] = N_Err[3+i];
    Z_He[i] = Z[3+i]; Z_Err_He[i] = Z_Err[3+i];
    R21_He[i] = R21[3+i]; R21_Err_He[i] = R21_Err[3+i];
  }
  
  TGraphErrors* gr1_R21_N_He = new TGraphErrors(2,N_He,R21_He,N_Err_He,R21_Err_He);
  TGraphErrors* gr1_R21_Z_He = new TGraphErrors(2,Z_He,R21_He,Z_Err_He,R21_Err_He);
  
  
  c1_N->cd(Index+1)->SetLogy(1);
  gr1_R21_N_H->Draw("A*");
  gr1_R21_N_H->GetXaxis()->SetLimits(0,6);
  gr1_R21_N_H->GetYaxis()->SetRangeUser(0.5,2);
  gr1_R21_N_He->Draw("same*");
  Set_1D_GraphStyle(gr1_R21_N_H,"N","R21");
  Set_1D_GraphStyle(gr1_R21_N_He,"N","R21");
  gr1_R21_N_H->SetMarkerStyle(25);
  gr1_R21_N_He->SetMarkerStyle(25);

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.08);
  sprintf(NameTem,"E_{kin}#in[%.1f,%.1f)",Ekin_BinEdge[0],Ekin_BinEdge[1]);
  latex->DrawLatex(2,0.7,NameTem);
  
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
  gr1_R21_Z_H->Draw("A*"); 
  gr1_R21_Z_H->GetXaxis()->SetLimits(0,6);
  gr1_R21_Z_H->GetYaxis()->SetRangeUser(0.5,2);
  gr1_R21_Z_He->Draw("same*");
  Set_1D_GraphStyle(gr1_R21_Z_H,"Z","R21");
  Set_1D_GraphStyle(gr1_R21_Z_He,"Z","R21");
  gr1_R21_Z_H->SetMarkerStyle(25);
  gr1_R21_Z_He->SetMarkerStyle(25);

  sprintf(NameTem,"E_{kin}#in[%.1f,%.1f)",Ekin_BinEdge[0],Ekin_BinEdge[1]);
  latex->DrawLatex(2,0.7,NameTem);
  
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
