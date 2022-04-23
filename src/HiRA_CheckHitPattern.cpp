#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include "HiRA_CheckHitPattern.hh"
ClassImp(HiRA_CheckHitPattern);

HiRA_CheckHitPattern::HiRA_CheckHitPattern(string SystemTagTem, string RunTagTem, string Hira_BadMap_VersionTem)
{
  SystemTag = SystemTagTem;
  RunTag = RunTagTem;
  Hira_BadMap_Version = Hira_BadMap_VersionTem;
  
  h1_WholeHira_Multi_Dis = 0;
  for(int i=0;i<2;i++)
  {
    h2_WholeHira_Theta_Phi_Lab[i] = 0;
    for(int j=0;j<12;j++)
    {
      h2_1Hira_Theta_Phi_Lab[j][i] = 0;
      h1_1Hira_Theta_HitCount[j][i] = 0;
      h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[j][i] = 0;
    }
  }
  IsActive_BadMap = 0;
  IsHistoInitialized = 0;
  Is_HiraPos_Applied=1;
}

HiRA_CheckHitPattern::~HiRA_CheckHitPattern()
{;}

void HiRA_CheckHitPattern::SetAnaTag(string SystemTagTem, string RunTagTem, string Hira_BadMap_VersionTem)
{
  SystemTag = SystemTagTem;
  RunTag = RunTagTem;
  Hira_BadMap_Version = Hira_BadMap_VersionTem;
}

void HiRA_CheckHitPattern::InitialHisto()
{
  if(IsHistoInitialized==0)
  {
    h1_WholeHira_Multi_Dis = new TH1D("h1_WholeHira_Multi_Dis",";Multi;Count",20,0,20);
    string ConfigTag[2] = {"noBadMap","withBadMap"};
    char NameTem[200];
    for(int i=0;i<2;i++)
    {
      h2_WholeHira_Theta_Phi_Lab[i] = new TH2D(("h2_WholeHira_Theta_Phi_Lab"+ConfigTag[i]).c_str(),";#theta_{Lab}(Deg.);#phi_{Lab}(Deg.)",600,20,80,1000,150,250);;
      for(int j=0;j<12;j++)
      {
        sprintf(NameTem,"h2_Hira%d_Theta_Phi_Lab_%s",j,ConfigTag[i].c_str());
        h2_1Hira_Theta_Phi_Lab[j][i] = new TH2D(NameTem,";#theta_{Lab}(Deg.);#phi_{Lab}(Deg.)",600,20,80,1000,150,250);
        sprintf(NameTem,"h1_Hira%d_Theta_HitCount_%s",j,ConfigTag[i].c_str());
        h1_1Hira_Theta_HitCount[j][i] = new TH1D(NameTem,";#theta_{Lab};Count",120,20,80);
        sprintf(NameTem,"h1_Hira%d_Theta_HitCount_NormalizedWithPixelNum_%s",j,ConfigTag[i].c_str());
        h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[j][i] = new TH1D(NameTem,";#theta_{Lab};Count",120,20,80);
      }
    }
  }
  else
  {
    h1_WholeHira_Multi_Dis->Reset();
    for(int i=0;i<2;i++)
    {
      h2_WholeHira_Theta_Phi_Lab[i]->Reset();
      for(int j=0;j<12;j++)
      {
        h2_1Hira_Theta_Phi_Lab[j][i]->Reset();
        h1_1Hira_Theta_HitCount[j][i]->Reset();
        h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[j][i]->Reset();
      }
    }
  }
  IsHistoInitialized = 1;
}

void HiRA_CheckHitPattern::ReadExpData(int FileNum, string ExpFileName[],string RootFilePathForStore)
{
  if(Hira_BadMapper==0) { cout<<"Attention: the Hira_BadMapper==0 !!!, you need to set the BadMapper!"<<endl; }
  
  //Initial the histogram for storing.
  InitialHisto();
  
  char NameTem[200];
  //reading the data.
  TChain* t1_Data = new TChain("E15190");
  int FileNum_Existed = 0;
  for(int i=0;i<FileNum;i++)
  {
    cout<<"Reading "<<ExpFileName[i]<<" : ";
    ifstream infile(ExpFileName[i]);
    if(!infile.good()) { cout<< " Not Existed! "<<endl; continue; }
    else { cout<< " Existed! "<<endl; infile.close(); }
    
    t1_Data->AddFile(ExpFileName[i].c_str());
    FileNum_Existed++;
  }
  cout<<" ---> Existed File Num: "<<FileNum_Existed<<endl;
  
  Int_t fmulti;
  double fTheta[200]; double fPhi[200]; 
  int fnumtel[200]; int fnumstripf[200]; int fnumstripb[200];
  
  t1_Data->SetMakeClass(1);
  t1_Data->SetBranchStatus("*",0);
  t1_Data->SetBranchAddress("HiRA.fmulti",&fmulti);
  t1_Data->SetBranchAddress("HiRA.fnumtel",fnumtel);
  t1_Data->SetBranchAddress("HiRA.fnumstripf",fnumstripf);
  t1_Data->SetBranchAddress("HiRA.fnumstripb",fnumstripb);
  t1_Data->SetBranchAddress("HiRA.fTheta",fTheta);
  t1_Data->SetBranchAddress("HiRA.fPhi",fPhi);
  
  int EvtNum = t1_Data->GetEntries();
  cout<<" EvtNum : "<<EvtNum<<endl;
  
  for(int iEvt = 0;iEvt<EvtNum;iEvt++)
  {
    //if(iEvt%1000000==0) { cout<<"iEvt: "<<iEvt<<endl; }
    if(iEvt%100000==0) { printProgress (1.0*iEvt/EvtNum); }
    t1_Data->GetEntry(iEvt);
    h1_WholeHira_Multi_Dis->Fill(fmulti);
    
    for(int iP=0;iP<fmulti;iP++)
    {
      bool IsBad_Hira = Hira_BadMapper->IsBad_Hira(fnumtel[iP]); 
      bool IsBad_CsI = Hira_BadMapper->IsBad_CsI(fnumtel[iP],fnumstripf[iP],fnumstripb[iP]);
      bool IsBad_StripX = Hira_BadMapper->IsBad_StripX(fnumtel[iP],fnumstripf[iP]);
      bool IsBad_StripY = Hira_BadMapper->IsBad_StripY(fnumtel[iP],fnumstripb[iP]);
      
      double Theta = fTheta[iP]*RadToDeg();
      double Phi = fPhi[iP]*RadToDeg();
      
      if(Is_HiraPos_Applied==1)
      {
        Theta = HiraPos->GetTheta(fnumtel[iP],fnumstripf[iP],fnumstripb[iP]);
        Phi = HiraPos->GetPhi(fnumtel[iP],fnumstripf[iP],fnumstripb[iP]);
      }
      
      h2_WholeHira_Theta_Phi_Lab[0]->Fill(Theta,Phi);
      h1_1Hira_Theta_HitCount[fnumtel[iP]][0]->Fill(Theta);
      h2_1Hira_Theta_Phi_Lab[fnumtel[iP]][0]->Fill(Theta,Phi);
      
      if(IsActive_BadMap==1 && IsBad_Hira==0 && IsBad_CsI==0 && IsBad_StripX==0 && IsBad_StripY==0)
      {
        h2_WholeHira_Theta_Phi_Lab[1]->Fill(Theta,Phi);
        h1_1Hira_Theta_HitCount[fnumtel[iP]][1]->Fill(Theta);
        h2_1Hira_Theta_Phi_Lab[fnumtel[iP]][1]->Fill(Theta,Phi);
      }
    }
  }
  
  //the below is for calculating the normalized hit distribution along the theta direction.
  for(int i=0;i<12;i++)
  {
    GetNormalized_CountNum(h1_1Hira_Theta_HitCount[i][0],h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][0],h2_1Hira_Theta_Phi_Lab[i][0]);
    GetNormalized_CountNum(h1_1Hira_Theta_HitCount[i][1],h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][1],h2_1Hira_Theta_Phi_Lab[i][1]);
  }
  
  //store the histogram for each run.
  string FileFroStore = RootFilePathForStore+"/f1_HitPattern_"+SystemTag+"_"+RunTag+"_"+Hira_BadMap_Version+".root";
  
  TFile* f1_Hira_HitPattern_Checking = new TFile(FileFroStore.c_str(),"update");
  h1_WholeHira_Multi_Dis->Write("",TObject::kOverwrite);
  h2_WholeHira_Theta_Phi_Lab[0]->Write("",TObject::kOverwrite);
  h2_WholeHira_Theta_Phi_Lab[1]->Write("",TObject::kOverwrite);
  for(int i=0;i<12;i++)
  {
    h1_1Hira_Theta_HitCount[i][0]->Write("",TObject::kOverwrite);
    h1_1Hira_Theta_HitCount[i][1]->Write("",TObject::kOverwrite);
    h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][0]->Write("",TObject::kOverwrite);
    h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][1]->Write("",TObject::kOverwrite);
    h2_1Hira_Theta_Phi_Lab[i][0]->Write("",TObject::kOverwrite);
    h2_1Hira_Theta_Phi_Lab[i][1]->Write("",TObject::kOverwrite);
  }
  f1_Hira_HitPattern_Checking->Close(); //this root file must be closed for the next run analysis.
}

void HiRA_CheckHitPattern::GetNormalized_CountNum(TH1D* h1_Count,TH1D* h1_Count_Normalized,TH2D* h2_Pixel_Dis)
{
  int XBinNum = h1_Count->GetNbinsX();
  double BinSize = h1_Count->GetBinCenter(2)-h1_Count->GetBinCenter(1);
  for(int iBin=1;iBin<=XBinNum;iBin++)
  {
    int PixelNum = 0;
    if(h1_Count->GetBinContent(iBin)>0)
    {
      double X0 = h1_Count->GetBinCenter(iBin)-0.5*BinSize;
      double X1 = h1_Count->GetBinCenter(iBin)+0.5*BinSize;
      //cout<<X0<<"  "<<X1<<endl;
      int Pixel_BinX_Num = h2_Pixel_Dis->GetNbinsX();
      int Pixel_BinY_Num = h2_Pixel_Dis->GetNbinsY();
      
      for(int i=1;i<=Pixel_BinX_Num;i++)
      {
        for(int j=1;j<=Pixel_BinY_Num;j++)
        {
          double Pixel_Center =  0;
          Pixel_Center = h2_Pixel_Dis->GetXaxis()->GetBinCenter(i);
          
          if( h2_Pixel_Dis->GetBinContent(i,j)>1 && Pixel_Center>=X0 && Pixel_Center<X1)
          { 
            //cout<<Pixel_Center_X<<"  ";
            PixelNum++; 
          }
        }
      }
      //cout<<endl;
    }
    //cout<<"PixelNum : "<<PixelNum<<endl;
    if(PixelNum>=1)
    {
      h1_Count_Normalized->SetBinContent(iBin,h1_Count->GetBinContent(iBin)/PixelNum);
    }
  }
}

void HiRA_CheckHitPattern::Draw_Info()
{
  if(h1_WholeHira_Multi_Dis!=0)
  { 
    TCanvas* c1_WholeHira_Multi_Dis = new TCanvas("c1_WholeHira_Multi_Dis","c1_WholeHira_Multi_Dis",1);
    h1_WholeHira_Multi_Dis->Draw("hist");
  }
  if(h2_WholeHira_Theta_Phi_Lab[0]!=0 && h2_WholeHira_Theta_Phi_Lab[1]!=0)
  {
    TCanvas* c1_WholeHira_Theta_Phi_Lab = new TCanvas("c1_WholeHira_Theta_Phi_Lab","c1_WholeHira_Theta_Phi_Lab",1);
    c1_WholeHira_Theta_Phi_Lab->Divide(2,1);
    c1_WholeHira_Theta_Phi_Lab->cd(1); h2_WholeHira_Theta_Phi_Lab[0]->Draw("colz");
    c1_WholeHira_Theta_Phi_Lab->cd(2); h2_WholeHira_Theta_Phi_Lab[1]->Draw("colz");
  }
  
  TCanvas* c1_1Hira_Theta_HitCount = new TCanvas("c1_1Hira_Theta_HitCount","c1_1Hira_Theta_HitCount",1);
  c1_1Hira_Theta_HitCount->Divide(4,3);
  TCanvas* c1_1Hira_Theta_HitCount_NormalizedWithPixelNum = new TCanvas("c1_1Hira_Theta_HitCount_NormalizedWithPixelNum","c1_1Hira_Theta_HitCount_NormalizedWithPixelNum",1);
  c1_1Hira_Theta_HitCount_NormalizedWithPixelNum->Divide(4,3);
  
  for(int i=0;i<12;i++)
  {
    c1_1Hira_Theta_HitCount->cd(i+1);
    h1_1Hira_Theta_HitCount[i][0] -> Draw("hist");
    h1_1Hira_Theta_HitCount[i][0] -> SetLineWidth(2); h1_1Hira_Theta_HitCount[i][0] -> SetLineColor(1);
    h1_1Hira_Theta_HitCount[i][1] -> Draw("histsame");
    h1_1Hira_Theta_HitCount[i][0] -> SetLineWidth(2); h1_1Hira_Theta_HitCount[i][1] -> SetLineColor(2);
    c1_1Hira_Theta_HitCount_NormalizedWithPixelNum->cd(i+1);
    h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][0] ->Draw("hist");
    h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][0] ->SetLineWidth(2); h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][0] ->SetLineColor(1);
    h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][1] ->Draw("histsame");
    h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][1] ->SetLineWidth(2); h1_1Hira_Theta_HitCount_NormalizedWithPixelNum[i][1] ->SetLineColor(2);
  }
  
}

void HiRA_CheckHitPattern::printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}
