#include "HiRA_ReactionLost.hh"
ClassImp(HiRA_ReactionLost);

HiRA_ReactionLost::HiRA_ReactionLost()
{
  Initial_ReactionLost_CorEff();
}

HiRA_ReactionLost::~HiRA_ReactionLost()
{;}

void HiRA_ReactionLost::Initial_ReactionLost_CorEff()
{
  ReactionLost_Cor_ParticleNum = 5; //currently, Sean only check the reaction lost correction up to 4He, the heavier particle will be taken as 4He.
  int ZTem[] = {1,1,1,2,2}; int ATem[] = {1,2,3,3,4};
  double A[] = {-1.772780E-4,-3.124610E-4,-2.480730E-4,-9.824110E-5,-8.952670E-5};
  double B[] = {-1.120190E-5,-6.187160E-6,-4.943390E-6,-7.446480E-7,-5.984310E-7};
  
  string ParticleName[] = {"P","D","T","3He","4He"};
  char NameTem[200];
  int LineColor_DiffPID[] = {1,2,3,4,6,7,8,9};
  
  for(int iPID = 0;iPID<ReactionLost_Cor_ParticleNum;iPID++)
  {
    ReactionLost_Cor_ParticleName[iPID] = ParticleName[iPID];
    ReactionLost_Cor_Z[iPID] = ZTem[iPID];
    ReactionLost_Cor_A[iPID] = ATem[iPID];
    sprintf(NameTem,"f1_ReactionLost_CorEff_%s",ParticleName[iPID].c_str());
    f1_ReactionLost_CorEff[iPID] = new TF1(NameTem,"TMath::Exp([0]*x*x+[1]*x)",0,500);
    f1_ReactionLost_CorEff[iPID]->SetParameter(0,B[iPID]);
    f1_ReactionLost_CorEff[iPID]->SetParameter(1,A[iPID]);
    f1_ReactionLost_CorEff[iPID]->SetLineColor(LineColor_DiffPID[iPID]);
    f1_ReactionLost_CorEff[iPID]->SetLineWidth(2);
  }
}

void HiRA_ReactionLost::Draw_ReactionEff()
{
  TCanvas* c1_ReactionEff = new TCanvas("c1_ReactionEff","c1_ReactionEff",1);
  c1_ReactionEff->cd();
  TLegend* legend_ReactionEff = new TLegend(0.6,0.6,0.8,0.8);
  
  for(int iPID = 0;iPID<ReactionLost_Cor_ParticleNum;iPID++)
  {
    if(iPID!=0) { f1_ReactionLost_CorEff[iPID] -> Draw("same"); }
    else { f1_ReactionLost_CorEff[iPID] ->Draw(); }
    legend_ReactionEff->AddEntry(f1_ReactionLost_CorEff[iPID],ReactionLost_Cor_ParticleName[iPID].c_str(),"l");
    f1_ReactionLost_CorEff[iPID]->GetXaxis()->SetTitle("Ekin(MeV)");
    f1_ReactionLost_CorEff[iPID]->GetYaxis()->SetTitle("Eff_{Reaction}");
  }
  legend_ReactionEff->Draw();
}

double HiRA_ReactionLost::Get_ReactionLost_CorEff(int Z, int A, double Ekin_Total_Lab)
{
  int Index = -1;
  for(int iPID=0;iPID<ReactionLost_Cor_ParticleNum;iPID++)
  {
    if(Z==ReactionLost_Cor_Z[iPID] && A==ReactionLost_Cor_A[iPID])
    {
      Index = iPID;
      break;
    }
  }
  if(Index!=-1) { return f1_ReactionLost_CorEff[Index]->Eval(Ekin_Total_Lab); }
  else { return f1_ReactionLost_CorEff[ReactionLost_Cor_ParticleNum-1]->Eval(Ekin_Total_Lab); } //the others are set to same with the last particle( always is alpha ).
}

TH1D* HiRA_ReactionLost::CorrESpec_perA(string HistoName, TH1D* h1_ESpec_Origin, int Z, int A)
{
  TH1D* h1_Corr = (TH1D*) h1_ESpec_Origin->Clone(HistoName.c_str());
  
  int XBinNum = h1_Corr->GetNbinsX();
  for(int iX=1;iX<=XBinNum;iX++)
  {
    double BinCenter = h1_Corr->GetBinCenter(iX);
    double BinErr = h1_Corr->GetBinError(iX);
    double BinContent = h1_Corr->GetBinContent(iX);
    
    double CorrEff = Get_ReactionLost_CorEff(Z,A,BinCenter*A);
    if(CorrEff==0) { continue; } 
    h1_Corr->SetBinContent(iX,BinContent/CorrEff);
    h1_Corr->SetBinError(iX,BinErr/CorrEff);
  }
  return h1_Corr;
}
