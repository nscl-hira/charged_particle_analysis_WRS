#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include "HiRA_GeoEff.hh"
ClassImp(HiRA_GeoEff);

HiRA_GeoEff::HiRA_GeoEff()
{
  h2_WholeHira_Theta_Phi_Lab = 0;
  h2_noBadMap_Theta_Phi_Lab = 0;
  h2_BadMap_Theta_Phi_Lab = 0;
  h1_BadMap_Theta_Lab_HitCount = 0;
  h1_BadMap_Theta_Lab_Eff = 0;
  h1_noBadMap_Theta_Lab_HitCount = 0;
  h1_noBadMap_Theta_Lab_Eff = 0;
  h1_HiraHit_Multi = 0;
  Hira_BadMapper = 0;
  f1_Results = 0;
  IsApplyPixelAngle = 1;
  
  //the below is the default setting of theta range and bin size.
  ThetaRange[0] = 20; ThetaRange[1] = 80; //unit:deg.
  ThetaBinSize = 0.1; //unit:deg.
}

HiRA_GeoEff::~HiRA_GeoEff()
{;}

#ifdef NPTOOL
void HiRA_GeoEff::ReadSimData(string SimFileName,string RootFileNameForStore)
{
  if(Hira_BadMapper==0) { cout<<"Attention: the Hira_BadMapper==0 !!!, you need to set the BadMapper!"<<endl; }
  TFile* f1 = new TFile(SimFileName.c_str(),"read");
  TTree* t1 = (TTree*) f1->Get("SimulatedTree");
  
  THiraData* DetData= new THiraData(); //attention: this class is defined in the NPTool, so the NPTool's environment has to be set before study the GeoEff.
  TInitialConditions* SrcData = new TInitialConditions();
  
  t1->SetBranchAddress("Hira",&DetData);
  t1->SetBranchAddress("InitialConditions",&SrcData);
  
  h2_WholeHira_Theta_Phi_Lab = 0;
  h2_BadMap_Theta_Phi_Lab = 0;
  h2_noBadMap_Theta_Phi_Lab = 0;
  h1_BadMap_Theta_Lab_HitCount = 0;
  h1_BadMap_Theta_Lab_Eff = 0;
  h1_noBadMap_Theta_Lab_HitCount = 0;
  h1_noBadMap_Theta_Lab_Eff = 0;
  
  h2_WholeHira_Theta_Phi_Lab = new TH2D("h2_WholeHira_Theta_Phi_Lab",";#theta_{Lab}(Deg.);#phi_{Lab}(Deg.)",720,0,180,1440,0,360);
  h2_BadMap_Theta_Phi_Lab = new TH2D("h2_BadMap_Theta_Phi_Lab",";#theta_{Lab}(Deg.);#phi_{Lab}(Deg.)",600,20,80,1000,150,250);
  h2_noBadMap_Theta_Phi_Lab = new TH2D("h2_noBadMap_Theta_Phi_Lab",";#theta_{Lab}(Deg.);#phi_{Lab}(Deg.)",600,20,80,1000,150,250);
  h1_BadMap_Theta_Lab_HitCount = new TH1D("h1_BadMap_Theta_Lab_HitCount",";#theta_{Lab}(Deg.);count",(int)((ThetaRange[1]-ThetaRange[0])/ThetaBinSize),ThetaRange[0],ThetaRange[1]);
  h1_noBadMap_Theta_Lab_HitCount = new TH1D("h1_noBadMap_Theta_Lab_HitCount",";#theta_{Lab}(Deg.);count",(int)((ThetaRange[1]-ThetaRange[0])/ThetaBinSize),ThetaRange[0],ThetaRange[1]);
  h1_BadMap_Theta_Lab_Eff = (TH1D*) h1_BadMap_Theta_Lab_HitCount->Clone("h1_BadMap_Theta_Lab_Eff");
  h1_noBadMap_Theta_Lab_Eff = (TH1D*) h1_noBadMap_Theta_Lab_HitCount->Clone("h1_noBadMap_Theta_Lab_Eff");
  
  h1_BadMap_Theta_Lab_Eff->Reset();
  h1_BadMap_Theta_Lab_Eff->GetYaxis()->SetTitle("Eff/d#theta(Prob./Deg.)");
  h1_noBadMap_Theta_Lab_Eff->Reset();
  h1_noBadMap_Theta_Lab_Eff->GetYaxis()->SetTitle("Eff/d#theta(Prob./Deg.)");
  
  h1_HiraHit_Multi = new TH1D("h1_HiraHit_Multi",";HiraHit_Multi;Count",20,0,20);
  
  int EvtNum = t1->GetEntries();
  cout<<" EvtNum : "<<EvtNum<<endl;
  int TotalSimEvtNum =  EvtNum;
  
  for(int iEvt=0;iEvt<EvtNum;iEvt++)
  {
    t1->GetEntry(iEvt);
    if(iEvt%100000==0) { printProgress (1.0*iEvt/EvtNum); }
    
    int Multi_Det = DetData->GetHiraStripXEMult();
    h1_HiraHit_Multi->Fill(Multi_Det);
    //if(Multi_Det<=0) {continue;}
    
    int Multi_Src = SrcData->GetParticleMultiplicity();
    for(int iSrc=0;iSrc<Multi_Src;iSrc++)
    {
      double Px = SrcData->GetMomentumDirectionX(iSrc);
      double Py = SrcData->GetMomentumDirectionY(iSrc);
      double Pz = SrcData->GetMomentumDirectionZ(iSrc);
      
      TVector3 P_Lab(Px,Py,Pz);
      
      double Theta = P_Lab.Theta()*RadToDeg(); //real theta, no digitized
      double Phi = P_Lab.Phi()*RadToDeg(); //real phi, no digitized
      
      if(Phi<0) { Phi+=360; }
      h2_WholeHira_Theta_Phi_Lab->Fill(Theta,Phi);
      
      if(DetData->GetHiraStripXTMult()==0 || DetData->GetHiraStripYEMult()==0) { continue; }
      int HiraIndex = DetData->GetHiraStripXEDetectorNbr(0);
      int X_StripIndex = DetData->GetHiraStripXEStripNbr(0);
      int Y_StripIndex = DetData->GetHiraStripYEStripNbr(0);
      // in NPTool, the X_StripIndex:0->31 is from -x to +x, but in the data analysis, it is from +x -> -x, so here for using the experimental dead silicon strip index, I flip it.
      X_StripIndex = 31-X_StripIndex;
      bool IsBad_Hira = Hira_BadMapper->IsBad_Hira(HiraIndex);
      bool IsBad_CsI = Hira_BadMapper->IsBad_CsI(HiraIndex,X_StripIndex,Y_StripIndex); //Attention: the X_StripIndex already be flipped.
      bool IsBad_StripX = Hira_BadMapper->IsBad_StripX(HiraIndex,X_StripIndex);
      bool IsBad_StripY = Hira_BadMapper->IsBad_StripY(HiraIndex,Y_StripIndex);
      
      double Theta_Digitized = HiraPos->GetTheta(HiraIndex, X_StripIndex, Y_StripIndex);
      double Phi_Digitized = HiraPos->GetPhi(HiraIndex, X_StripIndex, Y_StripIndex);;
      if(IsApplyPixelAngle==1)
      { 
        Theta = Theta_Digitized; 
        Phi = Phi_Digitized;
      }
      //with BadStripMap Filter.
      h2_noBadMap_Theta_Phi_Lab->Fill(Theta,Phi);
      h1_noBadMap_Theta_Lab_HitCount->Fill(Theta);
      if(IsBad_Hira==0 && IsBad_CsI==0 && IsBad_StripX==0 && IsBad_StripY==0)
      {
        h2_BadMap_Theta_Phi_Lab->Fill(Theta,Phi);
        h1_BadMap_Theta_Lab_HitCount->Fill(Theta);
      }
    }
  }
  Cal_GeoEff(h1_BadMap_Theta_Lab_HitCount,h1_BadMap_Theta_Lab_Eff,TotalSimEvtNum);
  Cal_GeoEff(h1_noBadMap_Theta_Lab_HitCount,h1_noBadMap_Theta_Lab_Eff,TotalSimEvtNum);
  
  f1_Results = new TFile(RootFileNameForStore.c_str(),"update");
  h2_WholeHira_Theta_Phi_Lab->Write("",TObject::kOverwrite);
  h2_BadMap_Theta_Phi_Lab->Write("",TObject::kOverwrite);
  h2_noBadMap_Theta_Phi_Lab->Write("",TObject::kOverwrite);
  h1_BadMap_Theta_Lab_HitCount->Write("",TObject::kOverwrite);
  h1_BadMap_Theta_Lab_Eff->Write("",TObject::kOverwrite);
  h1_noBadMap_Theta_Lab_HitCount->Write("",TObject::kOverwrite);
  h1_noBadMap_Theta_Lab_Eff->Write("",TObject::kOverwrite);
  h1_HiraHit_Multi->Write("",TObject::kOverwrite);
  f1_Results->Close();
  Draw_Info();
}
#endif

void HiRA_GeoEff::Cal_GeoEff(TH1D* h1_Count,TH1D* h1_GeoEff,int TotalSimEvtNum)
{
  int BinNum = h1_Count->GetNbinsX();
  double BinSize = h1_Count->GetBinCenter(2)-h1_Count->GetBinCenter(1);
  double Phi_Eff = 360.0/(250.0-150.0);
  double Theta_Eff = 1.0*TotalSimEvtNum/(80.0-20.0);
  
  for(int i=1;i<=BinNum;i++)
  {
    double BinContent = h1_Count->GetBinContent(i);
    double BinContent_RelErr = 0;
    if(BinContent!=0) { BinContent_RelErr = Sqrt(BinContent)/BinContent; }
    double GeoEff = BinContent/(Phi_Eff*Theta_Eff*BinSize);
    h1_GeoEff->SetBinContent(i,GeoEff);
    h1_GeoEff->SetBinError(i,GeoEff*BinContent_RelErr); //the relative error is same with the ;
  }
}

void HiRA_GeoEff::Draw_Info()
{
  string Hira_BadMap_Version = "";
  if(Hira_BadMapper!=0)
  {
    Hira_BadMap_Version = Hira_BadMapper->Get_BadMapper_Version();
  }
  
  char NameTem[200];
  if(h2_WholeHira_Theta_Phi_Lab!=0) 
  { 
    sprintf(NameTem,"c1_%s_WholeHira_Theta_Phi_Lab",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h2_WholeHira_Theta_Phi_Lab->Draw("colz"); 
  }
  if(h2_BadMap_Theta_Phi_Lab!=0) 
  { 
    sprintf(NameTem,"c1_%s_BadMap_Theta_Phi_Lab",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h2_BadMap_Theta_Phi_Lab->Draw("colz"); 
  }
  if(h2_noBadMap_Theta_Phi_Lab!=0) 
  { 
    sprintf(NameTem,"c1_%s_noBadMap_Theta_Phi_Lab",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h2_noBadMap_Theta_Phi_Lab->Draw("colz");
  }
  if(h1_BadMap_Theta_Lab_HitCount!=0) 
  {
    sprintf(NameTem,"c1_%s_BadMap_Theta_Lab_HitCount",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h1_BadMap_Theta_Lab_HitCount->Draw("hist"); 
  }
  if(h1_BadMap_Theta_Lab_Eff!=0) 
  { 
    sprintf(NameTem,"c1_%s_BadMap_Theta_Lab_Eff",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h1_BadMap_Theta_Lab_Eff->Draw("hist"); 
  }
  if(h1_noBadMap_Theta_Lab_HitCount!=0) 
  {
    sprintf(NameTem,"c1_%s_noBadMap_Theta_Lab_HitCount",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h1_noBadMap_Theta_Lab_HitCount->Draw("hist"); 
  }
  if(h1_noBadMap_Theta_Lab_Eff!=0) 
  { 
    sprintf(NameTem,"c1_%s_noBadMap_Theta_Lab_Eff",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h1_noBadMap_Theta_Lab_Eff->Draw("hist"); 
  }
  if(h1_HiraHit_Multi!=0) 
  { 
    sprintf(NameTem,"c1_%s_HiraHit_Multi",Hira_BadMap_Version.c_str());
    TCanvas* c1_tem = new TCanvas(NameTem,NameTem,1);
    h1_HiraHit_Multi->Draw("hist"); 
  }
}

void HiRA_GeoEff::ReadGeoEffHistogram(string RootFileName)
{
  f1_Results = new TFile(RootFileName.c_str(),"read");
  h1_BadMap_Theta_Lab_Eff = (TH1D*) f1_Results->Get("h1_BadMap_Theta_Lab_Eff");
  if(h1_BadMap_Theta_Lab_Eff!=0) { cout<<"Get the GeoEff Correction histogram!"<<endl; }
  h2_WholeHira_Theta_Phi_Lab = (TH2D*) f1_Results->Get("h2_WholeHira_Theta_Phi_Lab");
  h2_BadMap_Theta_Phi_Lab = (TH2D*) f1_Results->Get("h2_BadMap_Theta_Phi_Lab");
  h2_noBadMap_Theta_Phi_Lab = (TH2D*) f1_Results->Get("h2_noBadMap_Theta_Phi_Lab");
  h1_BadMap_Theta_Lab_HitCount = (TH1D*) f1_Results->Get("h1_BadMap_Theta_Lab_HitCount");
  h1_HiraHit_Multi = (TH1D*) f1_Results->Get("h1_HiraHit_Multi");
}

double HiRA_GeoEff::Get_GeoEff(double Theta_Lab)
{
  int BinIndex = h1_BadMap_Theta_Lab_Eff->FindBin(Theta_Lab);
  double Eff = h1_BadMap_Theta_Lab_Eff->GetBinContent(BinIndex);
return Eff;
}

void HiRA_GeoEff::printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}
