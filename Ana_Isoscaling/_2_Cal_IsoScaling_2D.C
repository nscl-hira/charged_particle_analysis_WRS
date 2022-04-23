#include "TMath.h"
using namespace TMath;

//this script will draw all the different particle together.
void Set_1D_HistoStyle(TH1D* h1, string XTitle, string YTitle);
void Set_2D_HistoStyle(TH2D* h2, string XTitle, string YTitle);
void Set_GraphStyle(TGraphErrors* gr1,string XTitle, string YTitle);
TH1D* GetTemperature(string NameTem,TH1D** h1_PtA_Lab);
TH2D* Rebin_ESpec(TH2D* h2_y_PtA_Lab, string NameTem, double* RapidityRange,double Rapidity_BinSize,double* PtARange,double PtA_BinSize);
void RemovePoint_LargeRelError(TH1D* h1,double RelErr);
void Cal_Alpha_Beta(int ParticleNum, double* Z, double* N, double* R21, double* R21_Err, double* FitResults);
void Cal_Alpha_Beta(int Index, TCanvas* c1_N, TCanvas* c1_Z, int ParticleNum, double* Z, double* N, double* R21, double* R21_Err,double* Alpha, double* Alpha_Err, double* Beta, double* Beta_Err, double* C, double* C_Err);

void Cal_Alpha_Beta_Avg(TH2D* h2_Rapidity_PtA_Alpha,double* Avg_PtASelect,double* Alpha_Avg,double* Alpha_Avg_Err);

void _2_Cal_IsoScaling_2D()
{
  gStyle->SetErrorX(0);
  
  string SystemTag_1 = "Ca40Ni58E140";
  string SystemTag_2 = "Ca48Ni64E140";
  
  string ThetaTag = "Noise_Theta_30_75";
  string EffCor_Tag = "GeoEff";
  string ImpactParaTag = "bHat_0.00_0.40";
  
  double RelErr_Cut = 0.1;
  double PtA_ShowRange[2] = {100,600};
  
  double RapidityRange[2] = {0,1.0}; double Rapidity_BinSize = 0.1;
  double PtARange[2] = {0,600.0}; double PtA_BinSize = 25;
  
  const int ParticleNum = 6;
  string ParticleName[]={"p","d","t","3He","4He","6He","6Li","7Li"};
  int LineColor[] = {1,1,1,2,2,2};
  int MarkerStyle[] = {20,24,21,20,24,21};
  
  TH2D* h2_y_PtA_Lab_Rebin[2][ParticleNum];
  TH2D* h2_y_PtA_Lab_Rebin_R21[ParticleNum];
  TH2D* h2_y_PtA_Lab[2][ParticleNum];
  
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
    sprintf(NameTem,"h2_y_PtA_Lab_%s_%s_%s_%s_Rebin",SystemTag_1.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_y_PtA_Lab_Rebin[0][iP] = Rebin_ESpec(h2_y_PtA_Lab[0][iP],NameTem,RapidityRange,Rapidity_BinSize,PtARange,PtA_BinSize);
    
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    f1_MergedESpec_2->cd();
    h2_y_PtA_Lab[1][iP] = (TH2D*) f1_MergedESpec_2->Get(NameTem);
    sprintf(NameTem,"h2_y_PtA_Lab_%s_%s_%s_%s",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_y_PtA_Lab[1][iP]->SetName(NameTem);
    sprintf(NameTem,"h2_y_PtA_Lab_%s_%s_%s_%s_Rebin",SystemTag_2.c_str(),ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_y_PtA_Lab_Rebin[1][iP] = Rebin_ESpec(h2_y_PtA_Lab[1][iP],NameTem,RapidityRange,Rapidity_BinSize,PtARange,PtA_BinSize);
    
    sprintf(NameTem,"h2_y_PtA_Lab_R21_%s_%s_%s_Rebin",ParticleName[iP].c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
    h2_y_PtA_Lab_Rebin_R21[iP] = (TH2D*) h2_y_PtA_Lab_Rebin[1][iP]->Clone(NameTem);
    h2_y_PtA_Lab_Rebin_R21[iP]->Divide(h2_y_PtA_Lab_Rebin[0][iP]);
  }
  
  //Draw the Pt/A~Rapidity
  sprintf(NameTem,"c1_Pt_A_Rapidity_%s_%s_%s",SystemTag_1.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Sys1->Divide(3,2,0,0);
  
  sprintf(NameTem,"c1_Pt_A_Rapidity_Rebin_%s_%s_%s",SystemTag_1.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Rebin_Sys1 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Rebin_Sys1->Divide(3,2);
  
  sprintf(NameTem,"c1_Pt_A_Rapidity_%s_%s_%s",SystemTag_2.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Sys2 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Sys2->Divide(3,2);
  
  sprintf(NameTem,"c1_Pt_A_Rapidity_Rebin_%s_%s_%s",SystemTag_2.c_str(),ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Rebin_Sys2 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Rebin_Sys2->Divide(3,2);
  
  sprintf(NameTem,"c1_Pt_A_Rapidity_Rebin_R21_%s_%s",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Rebin_R21 = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Rebin_R21->Divide(3,2);
  
  for(int i=0;i<ParticleNum;i++)
  {
    c1_Pt_A_Rapidity_Sys1->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab[0][i]->Draw("col");
    
    c1_Pt_A_Rapidity_Rebin_Sys1->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab_Rebin[0][i]->Draw("colz");
    
    c1_Pt_A_Rapidity_Sys2->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab[1][i]->Draw("colz");
    
    c1_Pt_A_Rapidity_Rebin_Sys2->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab_Rebin[1][i]->Draw("colz");
    
    c1_Pt_A_Rapidity_Rebin_R21->cd(i+1)->SetGrid(1,1);
    h2_y_PtA_Lab_Rebin_R21[i]->Draw("colz");
  }

  //the below is for calculating the alpha and Beta for each pixel.
  sprintf(NameTem,"h2_y_PtA_Lab_Alpha_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TH2D* h2_Rapidity_PtA_Alpha = (TH2D*)h2_y_PtA_Lab_Rebin[0][0]->Clone(NameTem);
  h2_Rapidity_PtA_Alpha->Reset();
  Set_2D_HistoStyle(h2_Rapidity_PtA_Alpha, "", "");
  
  sprintf(NameTem,"h2_y_PtA_Lab_Beta_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TH2D* h2_Rapidity_PtA_Beta = (TH2D*)h2_y_PtA_Lab_Rebin[0][0]->Clone(NameTem);
  h2_Rapidity_PtA_Beta->Reset();
  Set_2D_HistoStyle(h2_Rapidity_PtA_Beta, "", "");
  
  sprintf(NameTem,"h2_y_PtA_Lab_Alpha_Plus_Beta_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TH2D* h2_Rapidity_PtA_Alpha_Plus_Beta = (TH2D*)h2_y_PtA_Lab_Rebin[0][0]->Clone(NameTem);
  h2_Rapidity_PtA_Alpha_Plus_Beta->Reset();
  Set_2D_HistoStyle(h2_Rapidity_PtA_Alpha_Plus_Beta, "", "");
  
  double Z[] = {1,1,1,2,2,2}; double N[] = {0,1,2,1,2,4};
  TAxis* X_Axis = h2_Rapidity_PtA_Alpha->GetXaxis();
  TAxis* Y_Axis = h2_Rapidity_PtA_Alpha->GetYaxis();
  int XBinNum = X_Axis->GetNbins();
  int YBinNum = Y_Axis->GetNbins();
  
  double Alpha_Limit[2] = {99,-99}; double Beta_Limit[2] = {99,-99}; double Alpha_Plus_Beta_Limit[2] = {99,-99};
  int PointNum = 0;
  double Alpha[1000]; double AlphaErr[1000];
  double Beta[1000];  double BetaErr[1000];
  
  for(int iX=1;iX<=XBinNum;iX++)
  {
    for(int iY=1;iY<=YBinNum;iY++)
    {
      int IsValue = 1; //if each particle have R21, then the alpha and beta will be calculated.
      double R21[6] = {0}; double R21_Err[6] = {0};
      
      for(int iP=0;iP<ParticleNum-1;iP++)
      {
        R21[iP] = h2_y_PtA_Lab_Rebin_R21[iP]->GetBinContent(iX,iY);
        R21_Err[iP] = h2_y_PtA_Lab_Rebin_R21[iP]->GetBinError(iX,iY);
        if(R21[iP]==0) { IsValue=0; }
      }
      
      h2_Rapidity_PtA_Alpha->SetBinContent(iX,iY,-0);
      h2_Rapidity_PtA_Beta->SetBinContent(iX,iY,-0);
      h2_Rapidity_PtA_Alpha_Plus_Beta->SetBinContent(iX,iY,-0);
      
      if(IsValue==1) 
      {
        cout<<endl<<endl<<" iX: "<<iX<<" iY: "<<iY<<"("<<X_Axis->GetBinCenter(iX)<<","<<Y_Axis->GetBinCenter(iY)<<")"<<endl;
        cout.width(7);
        cout<<"      "<<"R21"<<"  "<<"R21_Err"<<endl;
        for(int iP=0;iP<ParticleNum-1;iP++)
        {
          cout.width(7);
          cout<<ParticleName[iP]<<"  "<<R21[iP]<<"  "<<R21_Err[iP]<<endl;
        }
        double Results[4] = {0,0,0,0}; //[Alpha, Alpha_Err, Beta, Beta_Err]
        Cal_Alpha_Beta(ParticleNum-1, Z, N, R21, R21_Err, Results);
        
        h2_Rapidity_PtA_Alpha->SetBinContent(iX,iY,Results[0]);
        h2_Rapidity_PtA_Alpha->SetBinError(iX,iY,Results[1]);
        h2_Rapidity_PtA_Beta->SetBinContent(iX,iY,Abs(Results[2]));
        h2_Rapidity_PtA_Beta->SetBinError(iX,iY,Abs(Results[3]));
        h2_Rapidity_PtA_Alpha_Plus_Beta->SetBinContent(iX,iY,Results[0]+Results[2]);
        h2_Rapidity_PtA_Alpha_Plus_Beta->SetBinError(iX,iY,Sqrt(Power(Results[1],2)+Power(Results[2],2)));
        
        double XCenter = X_Axis->GetBinCenter(iX); double YCenter = Y_Axis->GetBinCenter(iY);
        
        //if(XCenter>0.4 && XCenter<0.6)
        {
          Alpha[PointNum] = Results[0]; AlphaErr[PointNum] = Results[1];
          Beta[PointNum] = Results[2]; BetaErr[PointNum] = Results[3];
          PointNum++;
        }
        //select the limit of Alpha and beta for a better showing
        if(Results[0]>Alpha_Limit[1]) { Alpha_Limit[1] = Results[0]; }
        if(Results[0]<Alpha_Limit[0]) { Alpha_Limit[0] = Results[0]; }
        
        if(Results[2]>Beta_Limit[1]) { Beta_Limit[1] = Results[2]; }
        if(Results[2]<Beta_Limit[0]) { Beta_Limit[0] = Results[2]; }
        
        if((Results[0]+Results[2])>Alpha_Plus_Beta_Limit[1]) { Alpha_Plus_Beta_Limit[1] = (Results[0]+Results[2]); }
        if((Results[0]+Results[2])<Alpha_Plus_Beta_Limit[0]) { Alpha_Plus_Beta_Limit[0] = (Results[0]+Results[2]); } 
      }
    }
  }
  
  cout<<" Alpha_Limit : ["<<Alpha_Limit[0]<<","<<Alpha_Limit[1]<<"]"<<endl;
  cout<<" Beta_Limit : ["<<Beta_Limit[0]<<","<<Beta_Limit[1]<<"]"<<endl;
  cout<<" Alpha_Plus_Beta_Limit : ["<<Alpha_Plus_Beta_Limit[0]<<","<<Alpha_Plus_Beta_Limit[1]<<"]"<<endl;
  
  sprintf(NameTem,"c1_Pt_A_Rapidity_Alpha_Beta_%s_%s",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TCanvas* c1_Pt_A_Rapidity_Alpha_Beta = new TCanvas(NameTem,NameTem,1);
  c1_Pt_A_Rapidity_Alpha_Beta->Divide(2,2);
  
  //Alpha_Limit[0] = 0.1; Alpha_Limit[1] = 0.4;
  Beta_Limit[0] = 0.1; Beta_Limit[1] = 0.3;
  //Alpha_Plus_Beta_Limit[0] = 0.01; Alpha_Plus_Beta_Limit[1] = 0.10;
  
  double Delta = 0;//Alpha_Limit[1]-Alpha_Limit[0];
  c1_Pt_A_Rapidity_Alpha_Beta->cd(1); h2_Rapidity_PtA_Alpha->Draw("colz");
  h2_Rapidity_PtA_Alpha->GetZaxis()->SetRangeUser(Alpha_Limit[0]-0.1*Delta,Alpha_Limit[1]+0.1*Delta);
  h2_Rapidity_PtA_Alpha->GetYaxis()->SetRangeUser(PtA_ShowRange[0],PtA_ShowRange[1]);
  
  Delta = 0;//Beta_Limit[1]-Beta_Limit[0];
  c1_Pt_A_Rapidity_Alpha_Beta->cd(2); h2_Rapidity_PtA_Beta->Draw("colz");
  h2_Rapidity_PtA_Beta->GetZaxis()->SetRangeUser(Beta_Limit[0]-0.1*Delta,Beta_Limit[1]+0.1*Delta);
  h2_Rapidity_PtA_Beta->GetYaxis()->SetRangeUser(PtA_ShowRange[0],PtA_ShowRange[1]);
  
  Delta = 0;//Alpha_Plus_Beta_Limit[1] - Alpha_Plus_Beta_Limit[0];
  c1_Pt_A_Rapidity_Alpha_Beta->cd(3); h2_Rapidity_PtA_Alpha_Plus_Beta->Draw("colz");
  h2_Rapidity_PtA_Alpha_Plus_Beta->GetZaxis()->SetRangeUser(Alpha_Plus_Beta_Limit[0]-0.1*Delta,Alpha_Plus_Beta_Limit[1]+0.1*Delta);
  h2_Rapidity_PtA_Alpha_Plus_Beta->GetYaxis()->SetRangeUser(PtA_ShowRange[0],PtA_ShowRange[1]);
  
  //the below is for the correlation of alpha and beta.
  Delta = 0;
  c1_Pt_A_Rapidity_Alpha_Beta->cd(4); 
  sprintf(NameTem,"gr1_Corr_Alpha_Beta_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TGraphErrors* gr1_Corr_Alpha_Beta = new TGraphErrors(PointNum,Alpha,Beta,AlphaErr,BetaErr);
  gr1_Corr_Alpha_Beta->SetName(NameTem);
  gr1_Corr_Alpha_Beta->Draw("A*");
  gr1_Corr_Alpha_Beta->SetMarkerStyle(25);
  Set_GraphStyle(gr1_Corr_Alpha_Beta, "#alpha", "#beta");
  
  TF1* func1_pol1 = new TF1("func1_pol1","-x+[0]",-0.5,0.5);
  gr1_Corr_Alpha_Beta->Fit("func1_pol1","","",-0.5,0.5);
  
  TLine* RefLine = new TLine(-0.5,0.5,0.5,-0.5);
  RefLine->SetLineWidth(2);
  RefLine->SetLineColor(4);
  RefLine->SetLineStyle(2);
  RefLine->Draw("same");
  
  sprintf(NameTem,"./f1_AlphaBeta_%s_%s_%s.root",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag.c_str());
  TFile* f1_AlphaBeta_2D = new TFile(NameTem,"update");
  h2_Rapidity_PtA_Alpha->Write("",TObject::kOverwrite);
  h2_Rapidity_PtA_Beta->Write("",TObject::kOverwrite);
  f1_AlphaBeta_2D->Close();
  
  //the below is for calculating the avg of alpha and beta along the rapidity
  double Avg_PtASelect[2] = {200,300}; 
  double Alpha_Avg[100] = {0}; double Beta_Avg[100] = {0};
  double Alpha_Avg_Err[100] = {0}; double Beta_Avg_Err[100] = {0};
  
  //Cal_Alpha_Beta_Avg(h2_Rapidity_PtA_Alpha,Avg_PtASelect,Alpha_Avg,Alpha_Avg_Err);
  //Cal_Alpha_Beta_Avg(h2_Rapidity_PtA_Beta,Avg_PtASelect,Beta_Avg,Beta_Avg_Err);
  
  
}

void Cal_Alpha_Beta_Avg(TH2D* h2_Rapidity_PtA,double* PtASelect,double* Avg,double* Avg_Err)
{
  TAxis* X_Axis = h2_Rapidity_PtA->GetXaxis();
  TAxis* Y_Axis = h2_Rapidity_PtA->GetYaxis();
  int XBinNum = X_Axis->GetNbins();
  int YBinNum = Y_Axis->GetNbins();
  
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
      if(BinErr/BinContent>0.2) { h2_y_PtA_Lab_Rebin->SetBinContent(iX,iY,0); } //if the relative error is larger than 0.2, then set it to 0.
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
  FitResults[2] = Results[1]; FitResults[3] = ResultsErr[1];
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
  TF2* f2 = new TF2(NameTem,"[0]*x+[1]*y+[2]",-1,3,1,3);
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
