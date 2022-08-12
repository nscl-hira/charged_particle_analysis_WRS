#include "TMath.h"
using namespace TMath;

double Get_EvtNum_uBallDS_bHatCut(TH1D* h1_uBallDS_Multi,TH1D* h1_uBallMulti_bHat,double* bHat_Cut);
void ApplyGeoEff(Hira_GeoEff* Hira_GeoEfficiency,TH2D* h2_Ekin_Theta_Lab_tem);
void SetStyle_1DHisto(TH1D* h1);
TH2D* Get_Rapidity_PtA(TH2D* h2_Ekin_Theta_Lab, int ParticleA, double ParticleMass, double Rapidity_Beam_inLabFrame,double* ThetaCut, double* EkinCut);

void _0_Ana_MergeRun()
{
  bool IsDraw = 0;
  double bHat_Cut[2] = {0.0,0.4}; //[0,0.4)
  
  char NameTem[200];
  string DataPath = "./ESpec";
  string Hira_BadMap_FilePath = "./GeoEff";
  
  Hira_PosCali* HiraPos = new Hira_PosCali();
  HiraPos->Read_Hira_PixelAngle("./Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat"); //up to now, this pixel angle should not be changed.
  
  Hira_ReactionLost* Hira_ReactionLost_Corrector = new Hira_ReactionLost();
  //Hira_ReactionLost_Corrector->Draw_ReactionEff();
  
  if(HiraPos==0) { cout<<"HiraPos==0"<<endl; }
  if(Hira_ReactionLost_Corrector==0) { cout<<"ReactionLost_Corrector==0"<<endl; }
  
  string SystemTag = "Ca48Ni64E56"; 
  Exp_RunInfo* RunInfo = new Exp_RunInfo();
  RunInfo->Read_RunInfo("./config/RunInfo.data");
  RunInfo->PrintRunInfo(SystemTag);
  int RunNum = RunInfo->Get_RunNum(SystemTag);
  double Rapidity_Beam_inLabFrame = RunInfo->Get_BeamRapidity_Lab(SystemTag);
  
  double Mass_1u = 931.49410242;
  const int ParticleNum = 12;
  string ParticleName[ParticleNum] = {"p","d","t","3He","4He","6He","6Li","7Li","8Li","7Be","9Be","10Be"};
  int ParticleZ[ParticleNum] = {1,1,1,2,2,2,3,3,3,4,4,4};
  int ParticleN[ParticleNum] = {0,1,2,1,2,4,3,4,5,3,5,6};
  double ParticleMass[ParticleNum] = {938.272,1876.179,2809.514, //H
                      2809.496,3728.510,6.01888589*Mass_1u, //He
                      6.0151228874*Mass_1u,7.016003437*Mass_1u, 8.02248625*Mass_1u, //9.02679019*Mass_1u,//Li
                      7.01692872*Mass_1u, 9.01218307*Mass_1u,10.01353470*Mass_1u//Be
                      }; //unit: MeV/c2
  
  double ThetaCut[ParticleNum][2] = {{30,75},{30,75},{30,75},
                                     {30,75},{30,75},{30,75},
                                     {30,75},{30,75},{30,75},
                                     {30,75},{30,75},{30,75}};
                                     
  double EkinCut[ParticleNum][2] = {{20,198.0/1.0},{15,263.0/2.0},{12,312.0/3.0},
                                    {20,200},{18,200},{13,200},
                                    {22,200},{22,200},{22,200},
                                    {22,200},{22,200},{22,200}};
  
  TH2D* h2_Ekin_Theta_Lab_RawCount[ParticleNum];
  TH2D* h2_Ekin_Theta_Lab_Normalized[ParticleNum]; //Normalized with the total event number.
  TH2D* h2_Ekin_Theta_Lab_Normalized_GeoEff[ParticleNum]; //Normalized with the total event number and considering the geometry efficiency
  TH2D* h2_Rapidity_PtA_Lab_Normalized_GeoEff[ParticleNum]; 
  TH2D* h2_Rapidity_PtA_Lab_RawCount[ParticleNum];
  
  double ThetaRange_Checking[3][2] = {{20,44},{44,64},{64,80}};
  int Theta_LineColor[3] = {2,6,1};
  TH1D* h1_Ekin_Lab_GeoEff[200][3][ParticleNum]; //3 Theta_Lab:(theta<=44),(44<theta<=64),(theta>=64)
  
  sprintf(NameTem,"%s/%s/f1_MergedData_bHat_%.2f_%.2f.root",DataPath.c_str(),SystemTag.c_str(),bHat_Cut[0],bHat_Cut[1]);
  TFile* f1_MergedData = new TFile(NameTem,"recreate");
  f1_MergedData->cd();
  
  sprintf(NameTem,"c1_Ekin_Lab_QA_%s",SystemTag.c_str());
  TCanvas* c1_Ekin_Lab_QA = new TCanvas(NameTem,NameTem,1200,800);
  c1_Ekin_Lab_QA->Divide(3,2);
  
  for(int iP=0;iP<ParticleNum;iP++)
  {
    sprintf(NameTem,"h2_Ekin_Theta_Lab_RawCount_%s",ParticleName[iP].c_str());
    h2_Ekin_Theta_Lab_RawCount[iP] = new TH2D(NameTem,";Ekin_{lab}(MeV/A);Theta_{lab}(Deg.)",200,0,200,600,20,80);
    
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Normalized_%s",ParticleName[iP].c_str());
    h2_Ekin_Theta_Lab_Normalized[iP] = new TH2D(NameTem,";Ekin_{lab}(MeV/A);Theta_{lab}(Deg.)",200,0,200,600,20,80);
    
    sprintf(NameTem,"h2_Ekin_Theta_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    h2_Ekin_Theta_Lab_Normalized_GeoEff[iP] = new TH2D(NameTem,";Ekin_{lab}(MeV/A);Theta_{lab}(Deg.)",200,0,200,600,20,80);
    
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_RawCount_%s",ParticleName[iP].c_str());
    h2_Rapidity_PtA_Lab_RawCount[iP] = new TH2D(NameTem,";y/y_{Beam-Lab};Pt/A(MeV/c)",100,0,1,600,0,600);
    
    sprintf(NameTem,"h2_Rapidity_PtA_Lab_Normalized_GeoEff_%s",ParticleName[iP].c_str());
    h2_Rapidity_PtA_Lab_Normalized_GeoEff[iP] = new TH2D(NameTem,";y/y_{Beam-Lab};Pt/A(MeV/c)",100,0,1,600,0,600);
  }
  
  double EvtNum_uBallDS = 0;
  
  for(int iRun=0;iRun<RunNum;iRun++)
  {
    int RunNO = RunInfo->Get_RunNO(SystemTag,iRun);
    string BadMapVersion = RunInfo->Get_BadMapVersion(SystemTag,iRun);
    string Trigger = RunInfo->Get_Trigger(SystemTag,iRun);
    cout<<SystemTag<<"  "<<RunNO<<"  "<<BadMapVersion<<"  "<<Trigger<<endl;
    
    if(Trigger!="uBallDS_uBallHira_uBallNW") { continue; }
    
    string GeoEffFile = Hira_BadMap_FilePath+"/f1_GeoEff_"+BadMapVersion+".root";
    
    Hira_BadMap* Hira_BadMapper = new Hira_BadMap();
    Hira_BadMapper->Set_BadMapper_Version(Hira_BadMap_FilePath, BadMapVersion);
    
    Hira_GeoEff* Hira_GeoEfficiency = new Hira_GeoEff();
    Hira_GeoEfficiency->ReadGeoEffHistogram(GeoEffFile.c_str());
    
    if(Hira_BadMapper==0) { cout<<"Hira_BadMapper==0"<<endl; }
    if(Hira_GeoEfficiency==0) { cout<<"Hira_GeoEfficiency==0"<<endl; }
    
    sprintf(NameTem,"%s/%s/bHat_%.2f_%.2f/f1_Run%d.root",DataPath.c_str(),SystemTag.c_str(),bHat_Cut[0],bHat_Cut[1],RunNO);
    //cout<<" Process "<<NameTem<<endl;
    TFile* f1_Data = new TFile(NameTem,"read");
    
    for(int iP=0;iP<ParticleNum;iP++)
    {
      sprintf(NameTem,"h2_Ekin_Theta_Lab_%s",ParticleName[iP].c_str());
      TH2D* h2_Ekin_Theta_Lab_tem = (TH2D*) f1_Data->Get(NameTem);
      h2_Ekin_Theta_Lab_Normalized[iP]->Add(h2_Ekin_Theta_Lab_tem);
      h2_Ekin_Theta_Lab_RawCount[iP]->Add(h2_Ekin_Theta_Lab_tem);
      
      TH2D* h2_Rapidity_PtA_RawCount_tem = Get_Rapidity_PtA(h2_Ekin_Theta_Lab_tem,(ParticleZ[iP]+ParticleN[iP]),ParticleMass[iP],Rapidity_Beam_inLabFrame,ThetaCut[iP], EkinCut[iP]);
      h2_Rapidity_PtA_Lab_RawCount[iP]->Add(h2_Rapidity_PtA_RawCount_tem);
      h2_Rapidity_PtA_RawCount_tem->Delete();
      
      //Calculate the GeoEff correctted histogram, the add into the sum 2D histogram.
      ApplyGeoEff(Hira_GeoEfficiency,h2_Ekin_Theta_Lab_tem);
      h2_Ekin_Theta_Lab_Normalized_GeoEff[iP]->Add(h2_Ekin_Theta_Lab_tem);
      
      TH2D* h2_Rapidity_PtA_tem = Get_Rapidity_PtA(h2_Ekin_Theta_Lab_tem,(ParticleZ[iP]+ParticleN[iP]),ParticleMass[iP],Rapidity_Beam_inLabFrame,ThetaCut[iP], EkinCut[iP]);
      h2_Rapidity_PtA_Lab_Normalized_GeoEff[iP]->Add(h2_Rapidity_PtA_tem);
      h2_Rapidity_PtA_tem->Delete();
      
      //if(iP==0 || iP==4)
      {
        for(int iTheta=0;iTheta<3;iTheta++)
        {
          sprintf(NameTem,"h1_Ekin_Lab_GeoEff_%s_Run%d_ThetaLab_%.1f_%.1f",ParticleName[iP].c_str(),RunNO,ThetaRange_Checking[iTheta][0],ThetaRange_Checking[iTheta][1]);
          int Bin1 = (h2_Ekin_Theta_Lab_tem->GetYaxis())->FindBin(ThetaRange_Checking[iTheta][0]);
          int Bin2 = (h2_Ekin_Theta_Lab_tem->GetYaxis())->FindBin(ThetaRange_Checking[iTheta][1]);
          
          f1_MergedData->cd();
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP] = (TH1D*) h2_Ekin_Theta_Lab_tem->ProjectionX(NameTem,Bin1,Bin2);
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->SetLineColor(Theta_LineColor[iTheta]);
          //the below is for considering the solid angle.
          double delta_SolidAngle = TwoPi()*(Cos(ThetaRange_Checking[iTheta][0]*DegToRad())-Cos(ThetaRange_Checking[iTheta][1]*DegToRad()));
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Scale(1.0/delta_SolidAngle);
        }
      }
    }
    
    TH1D* h1_uBallMulti_uBallDS = (TH1D*) f1_Data->Get("h1_uBallMulti_uBallDS");
    TH1D* h1_uBallMulti_uBallDS_bHatCut = (TH1D*) f1_Data->Get("h1_uBallMulti_uBallDS_bHatCut");
    TH1D* h1_uBallMulti_bHat = (TH1D*) f1_Data->Get("h1_uBallMulti_bHat");
    
    double EvtNum_uBallDS_1Run = Get_EvtNum_uBallDS_bHatCut(h1_uBallMulti_uBallDS,h1_uBallMulti_bHat,bHat_Cut);
    cout<<"   --> Run: "<<RunNO<<"  "<<EvtNum_uBallDS_1Run<<"  "<<300.0*h1_uBallMulti_uBallDS_bHatCut->Integral(1,40)<<endl<<endl;
    
    for(int iTheta=0;iTheta<3;iTheta++) 
    {
      f1_MergedData->cd();
      double RebinValue = 5.0;
      for(int iP=0;iP<6;iP++)
      {
        h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Rebin(RebinValue);
        h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Scale(1.0/(EvtNum_uBallDS_1Run*RebinValue));
        if(iP==5) 
        {
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Rebin(RebinValue);
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Scale(1.0/RebinValue); 
        }
        c1_Ekin_Lab_QA->cd(iP+1)->SetLogy(1);
        if(iRun==0) 
        { 
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Draw();
          h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->GetYaxis()->SetRangeUser(1E-4,4E-2);
          if(iP==5)
          {
            h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->GetYaxis()->SetRangeUser(1E-6,1E-2);
          }
        }
        else { h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->Draw("same"); }
        SetStyle_1DHisto(h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]);
        h1_Ekin_Lab_GeoEff[iRun][iTheta][iP]->GetYaxis()->SetTitle("d^{2}M/dEd#Omega(Arb.)");
      }
      
      c1_Ekin_Lab_QA->Update();
    }
    
    EvtNum_uBallDS += EvtNum_uBallDS_1Run;
    
    f1_Data->Close();
    Hira_BadMapper->Delete();
    Hira_GeoEfficiency->Delete();
  }
  
  cout<<"EvtNum_uBallDS: "<<EvtNum_uBallDS<<endl;
  
  //Normalize the sum histogram.
  for(int iP=0;iP<ParticleNum;iP++)
  {
    h2_Ekin_Theta_Lab_Normalized[iP]->Scale(1.0/EvtNum_uBallDS);
    h2_Ekin_Theta_Lab_Normalized_GeoEff[iP]->Scale(1.0/EvtNum_uBallDS);
    h2_Rapidity_PtA_Lab_Normalized_GeoEff[iP]->Scale(1.0/EvtNum_uBallDS);
  }
  
  
  
  if(IsDraw==1)
  {
    TCanvas* c1_Ekin_Theta_Lab = new TCanvas("c1_Ekin_Theta_Lab","c1_Ekin_Theta_Lab",1);
    c1_Ekin_Theta_Lab->Divide(3,2);
    TCanvas* c1_Ekin_Theta_Lab_GeoEff = new TCanvas("c1_Ekin_Theta_Lab_GeoEff","c1_Ekin_Theta_Lab_GeoEff",1);
    c1_Ekin_Theta_Lab_GeoEff->Divide(3,2);
    TCanvas* c1_Rapidity_PtA_Lab_GeoEff = new TCanvas("c1_Rapidity_PtA_Lab_GeoEff","c1_Rapidity_PtA_Lab_GeoEff",1);
    c1_Rapidity_PtA_Lab_GeoEff->Divide(3,2);
    for(int iP=0;iP<ParticleNum/2;iP++)
    {
      c1_Ekin_Theta_Lab->cd(iP+1);
      h2_Ekin_Theta_Lab_Normalized[iP]->Draw("colz");
      
      c1_Ekin_Theta_Lab_GeoEff->cd(iP+1);
      h2_Ekin_Theta_Lab_Normalized_GeoEff[iP]->Draw("colz");
      
      c1_Rapidity_PtA_Lab_GeoEff->cd(iP+1);
      h2_Rapidity_PtA_Lab_Normalized_GeoEff[iP]->Draw("colz");
    }
  }
  
  f1_MergedData->cd();
  for(int iP=0;iP<ParticleNum;iP++)
  {
    h2_Ekin_Theta_Lab_Normalized[iP]->Write("",TObject::kOverwrite);
    h2_Ekin_Theta_Lab_RawCount[iP]->Write("",TObject::kOverwrite);
    h2_Ekin_Theta_Lab_Normalized_GeoEff[iP]->Write("",TObject::kOverwrite);
    h2_Rapidity_PtA_Lab_Normalized_GeoEff[iP]->Write("",TObject::kOverwrite);
    h2_Rapidity_PtA_Lab_RawCount[iP]->Write("",TObject::kOverwrite);
  }
}

double Get_EvtNum_uBallDS_bHatCut(TH1D* h1_uBallDS_Multi,TH1D* h1_uBallMulti_bHat,double* bHat_Cut)
{
  
  int BinIndex[2] = {40,1};
  if(bHat_Cut[0]==0) 
  {
    BinIndex[0] = 40;
    BinIndex[1] = h1_uBallMulti_bHat->FindLastBinAbove(bHat_Cut[1])+1;
  }
  else if(bHat_Cut[1]==1)
  {
    BinIndex[0] = h1_uBallMulti_bHat->FindLastBinAbove(bHat_Cut[0]);
    BinIndex[1] = 5; //now, 5 is set as minimum cut.
  }
  else
  {
    BinIndex[0] = h1_uBallMulti_bHat->FindLastBinAbove(bHat_Cut[0]);
    BinIndex[1] = h1_uBallMulti_bHat->FindLastBinAbove(bHat_Cut[1])+1;
  }
  
  int EvtNum = 0;
  EvtNum = h1_uBallDS_Multi->Integral(BinIndex[1],BinIndex[0]);
return 300.0*EvtNum;
}

//this function already include the error propagation
void ApplyGeoEff(Hira_GeoEff* Hira_GeoEfficiency,TH2D* h2_Ekin_Theta_Lab_tem)
{
  int BinXNum = h2_Ekin_Theta_Lab_tem->GetNbinsX();
  int BinYNum = h2_Ekin_Theta_Lab_tem->GetNbinsY();
  
  TAxis* Y_Axis = h2_Ekin_Theta_Lab_tem->GetYaxis();
  
  for(int iY=1;iY<BinYNum;iY++)
  {
    double BinCenter = Y_Axis->GetBinCenter(iY);
    double EffGeo = Hira_GeoEfficiency->Get_GeoEff(BinCenter); //Geometry Efficiency.
    if(EffGeo<0.001) { EffGeo=0.001; }
    for(int iX=1;iX<BinYNum;iX++)
    {
      double BinContent = h2_Ekin_Theta_Lab_tem->GetBinContent(iX,iY)/EffGeo;
      double BinErr = h2_Ekin_Theta_Lab_tem->GetBinError(iX,iY)/EffGeo;
      
      h2_Ekin_Theta_Lab_tem->SetBinContent(iX,iY,BinContent);
      h2_Ekin_Theta_Lab_tem->SetBinError(iX,iY,BinErr);
      
    }
  }
}

TH2D* Get_Rapidity_PtA(TH2D* h2_Ekin_Theta_Lab, int ParticleA, double ParticleMass, double Rapidity_Beam_inLabFrame,double* ThetaCut, double* EkinCut)
{
  int BinXNum = h2_Ekin_Theta_Lab->GetNbinsX();
  int BinYNum = h2_Ekin_Theta_Lab->GetNbinsY();
  
  TAxis* X_Axis = h2_Ekin_Theta_Lab->GetXaxis();
  TAxis* Y_Axis = h2_Ekin_Theta_Lab->GetYaxis();
  
  double XBinSize = X_Axis->GetBinCenter(2)-X_Axis->GetBinCenter(1);
  double YBinSize = Y_Axis->GetBinCenter(2)-Y_Axis->GetBinCenter(1);
  
  TH2D* h2_Rapidity_PtA_tem = new TH2D("h2_Rapidity_PtA_tem",";y/y_{Beam-Lab};Pt/A(MeV/c)",100,0,1,600,0,600);
  TAxis* X_Axis_Rapidity = h2_Rapidity_PtA_tem->GetXaxis();
  TAxis* Y_Axis_PtA = h2_Rapidity_PtA_tem->GetYaxis();
  
  //loop for the (Ekin,Theta_Lab)
  for(int iY=1;iY<BinYNum;iY++)
  {
    for(int iX=1;iX<BinXNum;iX++)
    {
      double Ekin = X_Axis->GetBinCenter(iX)+XBinSize*gRandom->Uniform(-0.5,0.5);
      double Theta = Y_Axis->GetBinCenter(iY)+YBinSize*gRandom->Uniform(-0.5,0.5);
      
      if(Theta<ThetaCut[0] || Theta>ThetaCut[1]) { continue; } //here is make a cut on Ekin/A
      if(Ekin<EkinCut[0] || Ekin>EkinCut[1]) { continue; }
      
      Ekin = Ekin*ParticleA;
      
      double BinContent_FromEkin = h2_Ekin_Theta_Lab->GetBinContent(iX,iY);
      double BinErr_FromEkin = h2_Ekin_Theta_Lab->GetBinError(iX,iY);
      
      if(BinContent_FromEkin<=0) { continue; }
      
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
      double BinError_Pre = h2_Rapidity_PtA_tem->GetBinError(Rapidity_Bin,PtA_Bin);
      
      h2_Rapidity_PtA_tem->Fill(Rapidity_Lab,PtA_Lab,BinContent_FromEkin);
      h2_Rapidity_PtA_tem->SetBinError(Rapidity_Bin,PtA_Bin,Sqrt(BinErr_FromEkin*BinErr_FromEkin+BinError_Pre*BinError_Pre));
      
      double BinContent_Updated = h2_Rapidity_PtA_tem->GetBinContent(Rapidity_Bin,PtA_Bin);
      
      /*
      if(iX==80&&iY==300) 
      { cout<<ParticleA<<"  "<<BinContent_FromEkin<<"  "<<BinErr_FromEkin<<"  "<<BinContent_Updated
            <<"  "<<Sqrt(BinErr_FromEkin*BinErr_FromEkin+BinError_Pre*BinError_Pre)<<endl; }
      */
    }
  }
return h2_Rapidity_PtA_tem;
}

void SetStyle_1DHisto(TH1D* h1)
{
  h1->SetLineWidth(3);
  h1->GetXaxis()->CenterTitle(1);
  h1->GetYaxis()->CenterTitle(1);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetYaxis()->SetTitleOffset(1.0);
}
