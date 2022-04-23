void _99_Study_AlphaBeta_Corr()
{
  gStyle->SetErrorX(0);
  
  string SystemTag_1 = "Ca40Sn112E56"; 
  string SystemTag_2 = "Ca48Sn124E56"; 
  string ImpactParaTag = "bHat_0.00_0.40";
  string EffCor_Tag = "GeoEff";
  
  char NameTem[200];
  sprintf(NameTem,"./f1_AlphaBeta_%s_%s_%s.root",SystemTag_1.c_str(),SystemTag_2.c_str(),ImpactParaTag.c_str());
  TFile* f1_AlphaBeta_2D = new TFile(NameTem,"update");
  sprintf(NameTem,"h2_y_PtA_Lab_Alpha_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TH2D* h2_Rapidity_PtA_Alpha = (TH2D*) f1_AlphaBeta_2D->Get(NameTem);
  sprintf(NameTem,"h2_y_PtA_Lab_Beta_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TH2D* h2_Rapidity_PtA_Beta = (TH2D*) f1_AlphaBeta_2D->Get(NameTem);
  
  TAxis* X_Axis = h2_Rapidity_PtA_Alpha->GetXaxis();
  TAxis* Y_Axis = h2_Rapidity_PtA_Alpha->GetYaxis();
  int XBinNum = X_Axis->GetNbins();
  int YBinNum = Y_Axis->GetNbins();
  
  int PointNum = 0;
  double Alpha[1000]; double AlphaErr[1000];
  double Beta[1000];  double BetaErr[1000];
  
  for(int iX=1;iX<=XBinNum;iX++)
  {
    for(int iY=1;iY<=YBinNum;iY++)
    {
      double XCenter = X_Axis->GetBinCenter(iX); double YCenter = Y_Axis->GetBinCenter(iY);
        
      //if(XCenter>0.4 && XCenter<0.6)
      {
        Alpha[PointNum] = Results[0]; AlphaErr[PointNum] = Results[1];
        Beta[PointNum] = Results[2]; BetaErr[PointNum] = Results[3];
        PointNum++;
      }
    }
  }
  
  TCanvas* c1_Alpha_Beta = new TCanvas(NameTem,NameTem,1);
  c1_Alpha_Beta->Divide(2,2);
  
  
  sprintf(NameTem,"gr1_Corr_Alpha_Beta_%s_%s_Rebin",ImpactParaTag.c_str(),EffCor_Tag.c_str());
  TGraphErrors* gr1_Corr_Alpha_Beta = new TGraphErrors(PointNum,Alpha,Beta,AlphaErr,BetaErr);
  gr1_Corr_Alpha_Beta->SetName(NameTem);
  gr1_Corr_Alpha_Beta->Draw("A*");
  gr1_Corr_Alpha_Beta->SetMarkerStyle(25);
  Set_GraphStyle(gr1_Corr_Alpha_Beta, "#alpha", "#beta");
  
}
