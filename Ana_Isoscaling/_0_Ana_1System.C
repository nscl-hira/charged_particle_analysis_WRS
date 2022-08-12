#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "TMath.h"
using namespace TMath;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage);
void Draw_Multi(string SystemTag,int RunNO, double* bHat_Cut, TH1D* h1_uBallMulti_uBallDS,TH1D* h1_uBallMulti_uBallDS_bHatCut, TH1D* h1_uBallMulti_uBallHira,TH1D* h1_uBallMulti_uBallNW,TH1D* h1_uBallMulti_bHat,TH1D* h1_uBallMulti_fb);
void SetStyle_1DHisto(TH1D* h1);
int GetPID(int A,int Z,int ParticleNum, int* ParticleZ, int* ParticleN);

void _0_Ana_1System()
{
  bool IsDraw = false;
  string ESpecFolder = "ESpec";
  double bHat_Cut[2] = {0.0,0.4}; //[0,0.4)
  
  char NameTem[200];
  //string DataPath = "/usr/local/DataAna/CaliData_Update";
  string DataPath = "/mnt/analysis/e15190/Brown/rootfiles";
  string Hira_BadMap_FilePath = "./GeoEff";
  
  Hira_PosCali* HiraPos = new Hira_PosCali();
  HiraPos->Read_Hira_PixelAngle("./Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat"); //up to now, this pixel angle should not be changed.
  
  Hira_ReactionLost* Hira_ReactionLost_Corrector = new Hira_ReactionLost();
  //Hira_ReactionLost_Corrector->Draw_ReactionEff();
  
  if(HiraPos==0) { cout<<"HiraPos==0"<<endl; }
  if(Hira_ReactionLost_Corrector==0) { cout<<"ReactionLost_Corrector==0"<<endl; }
  
  string SystemTag = "Ca40Ni58E56";
  Exp_RunInfo* RunInfo = new Exp_RunInfo();
  RunInfo->Read_RunInfo("./config/RunInfo.data");
  RunInfo->PrintRunInfo(SystemTag);
  int RunNum = RunInfo->Get_RunNum(SystemTag);
  
  double Mass_1u = 931.49410242; //MeV
  const int ParticleNum = 12;
  string ParticleName[ParticleNum] = {"p","d","t","3He","4He","6He","6Li","7Li","8Li","7Be","9Be","10Be"};
  int ParticleZ[ParticleNum] = {1,1,1,2,2,2,3,3,3,4,4,4};
  int ParticleN[ParticleNum] = {0,1,2,1,2,4,3,4,5,3,5,6};
  double ParticleMass[ParticleNum] = {938.272,1876.179,2809.514, //H
                      2809.496,3728.510,6.01888589*Mass_1u, //He
                      6.0151228874*Mass_1u,7.016003437*Mass_1u, 8.02248625*Mass_1u, //9.02679019*Mass_1u,//Li
                      7.01692872*Mass_1u, 9.01218307*Mass_1u,10.01353470*Mass_1u//Be
                      }; //unit: MeV/c2
  
  TH1D* h1_uBallMulti_uBallDS = new TH1D("h1_uBallMulti_uBallDS",";uBall Multi;Count",40,0.5,40.5);
  TH1D* h1_uBallMulti_uBallDS_bHatCut = new TH1D("h1_uBallMulti_uBallDS_bHatCut",";uBall Multi;Count",40,0.5,40.5);
  TH1D* h1_uBallMulti_uBallHira = new TH1D("h1_uBallMulti_uBallHira",";uBall Multi;Count",40,0.5,40.5);
  TH1D* h1_uBallMulti_uBallNW = new TH1D("h1_uBallMulti_uBallNW",";uBall Multi;Count",40,0.5,40.5);
  TH1D* h1_uBallMulti_bHat = new TH1D("h1_uBallMulti_bHat",";uBall Multi;bHat",40,0.5,40.5);
  TH1D* h1_uBallMulti_fb = new TH1D("h1_uBallMulti_fb",";uBall Multi;fb(fm)",40,0.5,40.5);
  
  SetStyle_1DHisto(h1_uBallMulti_uBallDS); h1_uBallMulti_uBallDS->SetLineColor(1); 
  SetStyle_1DHisto(h1_uBallMulti_uBallDS_bHatCut); h1_uBallMulti_uBallDS_bHatCut->SetLineColor(6); 
  SetStyle_1DHisto(h1_uBallMulti_uBallHira); h1_uBallMulti_uBallHira->SetLineColor(2);
  SetStyle_1DHisto(h1_uBallMulti_uBallNW); h1_uBallMulti_uBallNW->SetLineColor(4);
  SetStyle_1DHisto(h1_uBallMulti_bHat); 
  SetStyle_1DHisto(h1_uBallMulti_fb);
    
  TH2D* h2_Ekin_Theta_Lab[ParticleNum];
  for(int iP=0;iP<ParticleNum;iP++)
  {
    sprintf(NameTem,"h2_Ekin_Theta_Lab_%s",ParticleName[iP].c_str());
    h2_Ekin_Theta_Lab[iP] = new TH2D(NameTem,";Ekin_{lab}(MeV/A);Theta_{lab}(Deg.)",200,0,200,600,20,80);
  }
  
  for(int iRun=0;iRun<RunNum;iRun++)
  {
    int RunNO = RunInfo->Get_RunNO(SystemTag,iRun);
    string BadMapVersion = RunInfo->Get_BadMapVersion(SystemTag,iRun);
    string Trigger = RunInfo->Get_Trigger(SystemTag,iRun);
    cout<<SystemTag<<"  "<<RunNO<<"  "<<BadMapVersion<<"  "<<Trigger<<endl;
    
    if(!(Trigger=="uBallDS_uBallHira_uBallNW" || Trigger=="uBallDS_uBallHira")) { continue; }
    
    sprintf(NameTem,"%s/CalibratedData_%d.root",DataPath.c_str(),RunNO);
    cout<<" Process "<<NameTem<<endl;
    
    string GeoEffFile = Hira_BadMap_FilePath+"/f1_GeoEff_"+BadMapVersion+".root";
    
    Hira_BadMap* Hira_BadMapper = new Hira_BadMap();
    Hira_BadMapper->Set_BadMapper_Version(Hira_BadMap_FilePath, BadMapVersion);
  
    Hira_GeoEff* Hira_GeoEfficiency = new Hira_GeoEff();
    Hira_GeoEfficiency->ReadGeoEffHistogram(GeoEffFile.c_str());
    
    if(Hira_BadMapper==0) { cout<<"Hira_BadMapper==0"<<endl; }
    if(Hira_GeoEfficiency==0) { cout<<"Hira_GeoEfficiency==0"<<endl; }
    
    //the below is for processing this run
    TChain* t1_Data = new TChain("E15190");
    t1_Data->AddFile(NameTem);
    
    Int_t Hira_fmulti;
    double uBall_fb; double uBall_fbhat; int uBall_fmulti;
    double fTheta[200]; double fPhi[200]; 
    int fnumtel[200]; int fnumstripf[200]; int fnumstripb[200];
    double fEnergySifCal[200]; double fEnergySibCal[200]; 
    double fEnergycsifCal[200];
    double fKinEnergy[200];
    int fZId[200]; int fAId[200];
    double fZ[200]; double fA[200];
    double TDCTriggers_uBallDS_TRG;
    double TDCTriggers_uBallHira_TRG;
    double TDCTriggers_uBallNW_TRG;
    
    t1_Data->SetMakeClass(1);
    t1_Data->SetBranchStatus("*",0);
    t1_Data->SetBranchAddress("uBall.fb",&uBall_fb);
    t1_Data->SetBranchAddress("uBall.fbhat",&uBall_fbhat);
    t1_Data->SetBranchAddress("uBall.fmulti",&uBall_fmulti);
    t1_Data->SetBranchAddress("HiRA.fmulti",&Hira_fmulti);
    t1_Data->SetBranchAddress("HiRA.fnumtel",fnumtel);
    t1_Data->SetBranchAddress("HiRA.fnumstripf",fnumstripf);
    t1_Data->SetBranchAddress("HiRA.fnumstripb",fnumstripb);
    t1_Data->SetBranchAddress("HiRA.fTheta",fTheta); //Attention: This fTheta is not used, this is calculated by the HiraPos.
    t1_Data->SetBranchAddress("HiRA.fPhi",fPhi);  //Attention: This fTheta is not used, this is calculated by the HiraPos.
    
    t1_Data->SetBranchAddress("HiRA.fKinEnergy",fKinEnergy);
  //  t1_Data->SetBranchAddress("HiRA.fZId",fZId);
  //  t1_Data->SetBranchAddress("HiRA.fAId",fAId);
    t1_Data->SetBranchAddress("HiRA.fZ",fZ); //20200425, Kuan's calibrated data does not include fZId and fAId, so I get it by floor(*+0.5);
    t1_Data->SetBranchAddress("HiRA.fA",fA);
    t1_Data->SetBranchAddress("TDCTriggers.uBall_DS_TRG",&TDCTriggers_uBallDS_TRG);
    t1_Data->SetBranchAddress("TDCTriggers.uBallHiRA_TRG",&TDCTriggers_uBallHira_TRG);
    t1_Data->SetBranchAddress("TDCTriggers.uBallNW_TRG",&TDCTriggers_uBallNW_TRG);
    
    int EvtNum = t1_Data->GetEntries();
    cout<<" EvtNum : "<<EvtNum<<endl;
    int EvtNum_1 = 0;
    
    h1_uBallMulti_uBallDS->Reset();
    h1_uBallMulti_uBallHira->Reset();
    h1_uBallMulti_uBallNW->Reset();
    h1_uBallMulti_bHat->Reset();
    h1_uBallMulti_fb->Reset();
    h1_uBallMulti_uBallDS_bHatCut->Reset();
    for(int iP=0;iP<ParticleNum;iP++) { h2_Ekin_Theta_Lab[iP]->Reset(); }
    
    for(int iEvt = 0;iEvt<EvtNum;iEvt++)
    {
      if(iEvt%100000==0) { printProgress (1.0*iEvt/EvtNum); }
      t1_Data->GetEntry(iEvt);
      
      if(TDCTriggers_uBallDS_TRG>-9990) { h1_uBallMulti_uBallDS->Fill(uBall_fmulti); }
      if(TDCTriggers_uBallHira_TRG>-9990) { h1_uBallMulti_uBallHira->Fill(uBall_fmulti); }
      if(TDCTriggers_uBallNW_TRG>-9990) { h1_uBallMulti_uBallNW->Fill(uBall_fmulti); }
      
      if(uBall_fmulti<5) { continue; } 
      
      //if(iEvt<100000) 
      {
        int BinIndex = h1_uBallMulti_bHat->FindBin(uBall_fmulti);
        h1_uBallMulti_bHat->SetBinContent(BinIndex,uBall_fbhat);
        BinIndex = h1_uBallMulti_fb->FindBin(uBall_fmulti);
        h1_uBallMulti_fb->SetBinContent(BinIndex,uBall_fb);
      }
      
      if(uBall_fbhat<bHat_Cut[0] || uBall_fbhat>=bHat_Cut[1]) { continue; }
      if(TDCTriggers_uBallDS_TRG>-9990) { h1_uBallMulti_uBallDS_bHatCut->Fill(uBall_fmulti); }
      //if(TDCTriggers_uBallHira_TRG<-9990 && Hira_fmulti>0) { cout<<iEvt<<"  "<<TDCTriggers_uBallDS_TRG<<"  "<<TDCTriggers_uBallHira_TRG<<"  "<<TDCTriggers_uBallNW_TRG<<"  "<<Hira_fmulti<<endl; EvtNum_1++;}
      
      for(int iP=0;iP<Hira_fmulti;iP++)
      {
        bool IsBad_Hira = Hira_BadMapper->IsBad_Hira(fnumtel[iP]);
        bool IsBad_CsI = Hira_BadMapper->IsBad_CsI(fnumtel[iP],fnumstripf[iP],fnumstripb[iP]);
        bool IsBad_StripX = Hira_BadMapper->IsBad_StripX(fnumtel[iP],fnumstripf[iP]);
        bool IsBad_StripY = Hira_BadMapper->IsBad_StripY(fnumtel[iP],fnumstripb[iP]);
        
        fZId[iP] = floor(fZ[iP]+0.5);
        fAId[iP] = floor(fA[iP]+0.5);
        
        if(IsBad_Hira==0 && IsBad_CsI==0 && IsBad_StripX==0 && IsBad_StripY==0)
        {
          int PID = GetPID(fAId[iP],fZId[iP],ParticleNum,ParticleZ,ParticleN);
          if(PID==-1) { continue; }
          
          //if(fAId[iP]==1 && fZId[iP]==1) //select proton
          {
            double P_Mag_Lab = Sqrt(fKinEnergy[iP]*fKinEnergy[iP]+2*fKinEnergy[iP]*ParticleMass[PID]);
            TVector3 P_3D_Lab(0,0,P_Mag_Lab);
            P_3D_Lab.SetTheta(HiraPos->GetTheta(fnumtel[iP],fnumstripf[iP],fnumstripb[iP])*DegToRad());
            P_3D_Lab.SetPhi(HiraPos->GetPhi(fnumtel[iP],fnumstripf[iP],fnumstripb[iP])*DegToRad());
            
            double Theta_Lab = P_3D_Lab.Theta()*RadToDeg();
            double Ekin_Lab = fKinEnergy[iP];
            
            h2_Ekin_Theta_Lab[PID]->Fill(Ekin_Lab/fAId[iP],Theta_Lab);
            
            /*
            double Pt_Lab = P_3D_Lab.Pt();//get the transverse momentum.
            double Pt_A_Lab = Pt_Lab/fAId[iP];
            
            TLorentzVector P_4D_CM(P_3D_Lab,ParticleMass[0]+fKinEnergy[iP]);
            double EffGeo = Hira_GeoEfficiency->Get_GeoEff(P_3D_Lab.Theta()*RadToDeg()); //Geometry Efficiency.
            double ReactionEff = Hira_ReactionLost_Corrector->Get_ReactionLost_CorEff(fAId[iP], fAId[iP], fKinEnergy[iP]); //Energy lost efficiency.
            */
            
          }
        }
      }//particles in the Hira.
      
    }//events of this run
    
    //cout<<"EvtNum_1: "<<EvtNum_1<<endl;
    t1_Data->Delete();
    if(IsDraw==1)
    {
      Draw_Multi(SystemTag,RunNO,bHat_Cut,
                    h1_uBallMulti_uBallDS,h1_uBallMulti_uBallDS_bHatCut,
                    h1_uBallMulti_uBallHira,h1_uBallMulti_uBallNW,
                    h1_uBallMulti_bHat,h1_uBallMulti_fb);
    }
    
    bool IsFolderExist = 1;
    //the below is for checking the existance of a folder
    sprintf(NameTem,"./%s",ESpecFolder.c_str());
    if(access(NameTem, 0)==-1) 
    { 
      mkdir(NameTem,S_IRWXU);
      sprintf(NameTem,"./%s/%s",ESpecFolder.c_str(),SystemTag.c_str());
      mkdir(NameTem,S_IRWXU);
      sprintf(NameTem,"./%s/%s/bHat_%.2f_%.2f",ESpecFolder.c_str(),SystemTag.c_str(),bHat_Cut[0],bHat_Cut[1]);
      mkdir(NameTem,S_IRWXU);
    }
    else
    {
      sprintf(NameTem,"./%s/%s",ESpecFolder.c_str(),SystemTag.c_str());
      if(access(NameTem, 0)==-1) 
      { 
        mkdir(NameTem,S_IRWXU);
        sprintf(NameTem,"./%s/%s/bHat_%.2f_%.2f",ESpecFolder.c_str(),SystemTag.c_str(),bHat_Cut[0],bHat_Cut[1]);
        mkdir(NameTem,S_IRWXU);
      }
      else
      {
        sprintf(NameTem,"./%s/%s/bHat_%.2f_%.2f",ESpecFolder.c_str(),SystemTag.c_str(),bHat_Cut[0],bHat_Cut[1]);
        if(access(NameTem, 0)==-1) { mkdir(NameTem,S_IRWXU); }
      }
    }
    
    sprintf(NameTem,"./%s/%s/bHat_%.2f_%.2f/f1_Run%d.root",ESpecFolder.c_str(),SystemTag.c_str(),bHat_Cut[0],bHat_Cut[1],RunNO);
    cout<<"Writting into : "<<NameTem<<endl;
    TFile* f1_AnaResults = new TFile(NameTem,"update");
    h1_uBallMulti_uBallDS->Write("",TObject::kOverwrite);
    h1_uBallMulti_uBallDS_bHatCut->Write("",TObject::kOverwrite);
    h1_uBallMulti_uBallHira->Write("",TObject::kOverwrite);
    h1_uBallMulti_uBallNW->Write("",TObject::kOverwrite);
    h1_uBallMulti_bHat->Write("",TObject::kOverwrite);
    h1_uBallMulti_fb->Write("",TObject::kOverwrite);
    for(int iP=0;iP<ParticleNum;iP++) { h2_Ekin_Theta_Lab[iP]->Write("",TObject::kOverwrite); }
    f1_AnaResults->Close();
    
  }//all runs of this system
  
}

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
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

int GetPID(int A, int Z, int ParticleNum, int* ParticleZ, int* ParticleN)
{
  int PID = -1;
  for(int i=0;i<ParticleNum;i++)
  {
    if(Z==ParticleZ[i] && A == (ParticleZ[i]+ParticleN[i])) { PID = i; break; }
  }
  
return PID;
}

void Draw_Multi(string SystemTag,int RunNo,double* bHat_Cut,TH1D* h1_uBallMulti_uBallDS,TH1D* h1_uBallMulti_uBallDS_bHatCut, TH1D* h1_uBallMulti_uBallHira,TH1D* h1_uBallMulti_uBallNW,TH1D* h1_uBallMulti_bHat,TH1D* h1_uBallMulti_fb)
{
  char NameTem[200];
  sprintf(NameTem,"c1_MultiDis_DiffTrigger_%s_run%d",SystemTag.c_str(),RunNo);
  TCanvas* c1_Multi_DiffTrigger = new TCanvas(NameTem,NameTem,1200,800);
  c1_Multi_DiffTrigger->Divide(4,2);
  
  sprintf(NameTem,"h1_uBallMulti_uBallDS_run%d",RunNo);
  TH1D* h1_uBallMulti_uBallDS_tem = (TH1D* )h1_uBallMulti_uBallDS->Clone(NameTem);
  sprintf(NameTem,"h1_uBallMulti_uBallHira_run%d",RunNo);
  TH1D* h1_uBallMulti_uBallHira_tem = (TH1D*) h1_uBallMulti_uBallHira->Clone(NameTem);
  sprintf(NameTem,"h1_uBallMulti_uBallNW_run%d",RunNo);
  TH1D* h1_uBallMulti_uBallNW_tem = (TH1D*) h1_uBallMulti_uBallNW->Clone(NameTem);
  c1_Multi_DiffTrigger->cd(1); h1_uBallMulti_uBallDS_tem->Draw();
  c1_Multi_DiffTrigger->cd(2); h1_uBallMulti_uBallHira_tem->Draw();
  c1_Multi_DiffTrigger->cd(3); h1_uBallMulti_uBallNW_tem->Draw();
  c1_Multi_DiffTrigger->cd(4)->SetGrid(1,1); 
  
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
  
  double BinCenter[2] = {h1_uBallMulti_bHat->GetBinCenter(BinIndex[0]),h1_uBallMulti_bHat->GetBinCenter(BinIndex[1])};
  cout<<endl<<"bHat: ["<<bHat_Cut[0]<<","<<bHat_Cut[1]
            <<") ; BinIndex:[ "<<BinIndex[0]<<","<<BinIndex[1]
            <<"]; BinCenter: ["<<BinCenter[0]<<","<<BinCenter[1]<<"]"<<endl;
  
  //Attention: in root, TH1::Integral(Bin1, Bin2), this integral include the content of bin1 and bin2, means[bin1, bin2]
   
  sprintf(NameTem,"h1_uBallMulti_uBallDS_Normalized_%s_run%d",SystemTag.c_str(),RunNo);
  TH1D* h1_uBallMulti_uBallDS_Normalized = (TH1D*) h1_uBallMulti_uBallDS->Clone(NameTem);
  h1_uBallMulti_uBallDS_Normalized->Scale(1.0/h1_uBallMulti_uBallDS_Normalized->Integral(BinIndex[1],BinIndex[0]));
  h1_uBallMulti_uBallDS_Normalized->Draw();
  
  sprintf(NameTem,"h1_uBallMulti_uBallHira_Normalized_%s_run%d",SystemTag.c_str(),RunNo);
  TH1D* h1_uBallMulti_uBallHira_Normalized = (TH1D*) h1_uBallMulti_uBallHira->Clone(NameTem);
  h1_uBallMulti_uBallHira_Normalized->Scale(1.0/h1_uBallMulti_uBallHira_Normalized->Integral(BinIndex[1],BinIndex[0]));
  h1_uBallMulti_uBallHira_Normalized->Draw("same");
  
  sprintf(NameTem,"h1_uBallMulti_uBallNW_Normalized_%s_run%d",SystemTag.c_str(),RunNo);
  TH1D* h1_uBallMulti_uBallNW_Normalized = (TH1D*) h1_uBallMulti_uBallNW->Clone(NameTem);
  h1_uBallMulti_uBallNW_Normalized->Scale(1.0/h1_uBallMulti_uBallNW_Normalized->Integral(BinIndex[1],BinIndex[0]));
  h1_uBallMulti_uBallNW_Normalized->Draw("same");
  
  TLegend* legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->AddEntry(h1_uBallMulti_uBallDS,"uBallDS","lp");
  legend->AddEntry(h1_uBallMulti_uBallHira,"uBallHira","lp");
  legend->AddEntry(h1_uBallMulti_uBallNW,"uBallNW","lp");
  legend->Draw();
  
  c1_Multi_DiffTrigger->cd(5)->SetGrid(1,1);
  h1_uBallMulti_bHat->Draw();
  c1_Multi_DiffTrigger->cd(6)->SetGrid(1,1);
  h1_uBallMulti_fb->Draw();
  
  c1_Multi_DiffTrigger->cd(7)->SetGrid(1,1);
  sprintf(NameTem,"h1_uBallMulti_uBallDS_bHatCut_%.1f_%.1f_run%d",bHat_Cut[0],bHat_Cut[1],RunNo);
  TH1D* h1_uBallMulti_uBallDS_bHatCut_tem = (TH1D* )h1_uBallMulti_uBallDS_bHatCut->Clone(NameTem);
  h1_uBallMulti_uBallDS_bHatCut_tem->Draw();
  
  c1_Multi_DiffTrigger->cd(8)->SetGrid(1,1);
  TH1D* h1_uBallMulti_uBallDS_cummulative = (TH1D*)h1_uBallMulti_uBallDS_Normalized->GetCumulative(0,"_cumulative");
  TH1D* h1_uBallMulti_uBallHira_cummulative = (TH1D*)h1_uBallMulti_uBallHira_Normalized->GetCumulative(0,"_cumulative");
  TH1D* h1_uBallMulti_uBallNW_cummulative = (TH1D*)h1_uBallMulti_uBallNW_Normalized->GetCumulative(0,"_cumulative");
  h1_uBallMulti_uBallDS_cummulative->SetMaximum(1.0);
  h1_uBallMulti_uBallHira_cummulative->SetMaximum(1.0);
  h1_uBallMulti_uBallNW_cummulative->SetMaximum(1.0);
  h1_uBallMulti_uBallDS_cummulative->GetXaxis()->SetRangeUser(BinIndex[1]-2,BinIndex[0]+1);
  h1_uBallMulti_uBallDS_cummulative->Draw("hist");
  h1_uBallMulti_uBallHira_cummulative->Draw("same");
  h1_uBallMulti_uBallNW_cummulative->Draw("same");
  c1_Multi_DiffTrigger->Update();
}
