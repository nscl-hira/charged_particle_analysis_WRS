#include "TMath.h"
using namespace TMath;

void Cal_PixelAngle()
{
  TVector3 BeamSpotPos(0,0,0);
  string BeamSpotPos_Tag = "BeamPos_0_0_0";
  bool IsDrawWithSim = 1;
  
  ofstream outfile(("PixelAngle_"+BeamSpotPos_Tag+".dat").c_str());
  outfile<<"Tele EF EB theta(degree) phi(degree)"<<endl;
  TH2D* h2_PixelAngle = new TH2D(("h2_PixelAngle_"+BeamSpotPos_Tag).c_str(),";#theta_{Lab}(Deg.);#phi_{Lab}(Deg.)",600,20,80,1000,150,250);
  
  //the below is for calculating the basic vector of local coordinate.
  TVector3 LocalVectorUnit[12][3];
  TVector3 LocalFrame_Origin[12];
  
  ifstream infile("romerArmData.detector");
  for(int i=0;i<12;i++)
  {
    string DetTag; int HiraIndex;
    infile>>DetTag>>HiraIndex;
    if(HiraIndex!=i) { cout<<"Read error!"<<endl; }
    string PointTag; double PosX[4]; double PosY[4]; double PosZ[4];
    TVector3 PosTem[4];
    for(int j=0;j<4;j++)
    {
      infile>>PointTag>>PosX[j]>>PosY[j]>>PosZ[j];
      PosTem[j].SetXYZ(PosX[j],PosY[j],PosZ[j]);
    }
    
    TVector3 w = 1.0/4.0*(PosTem[0]+PosTem[1]+PosTem[2]+PosTem[3]);
    TVector3 v = PosTem[0]-PosTem[2];
    TVector3 u = v.Cross(w);
    v = w.Cross(u);
    
    LocalFrame_Origin[i] = w;
    LocalVectorUnit[i][0] = u.Unit();
    LocalVectorUnit[i][1] = v.Unit();
    LocalVectorUnit[i][2] = w.Unit();
  }
  
  //Loop all the strip of Hira telescope and write the global theta phi in a txt file as a map.
  char NameTem[200];
  for(int iHira=0;iHira<12;iHira++)
  {
    for(int iX=0;iX<32;iX++)
    {
      for(int iY=0;iY<32;iY++)
      {
         double posX = 0.5*64-0.5*2.0-iX*2.0;
	 double posY = (iY + 0.5 )*2.0 - 0.5*64;
	 double posZ = 7+6+0.75; //unit: mm, dist-Mylar-to-1stSi, dist-1stSi-to-2ndSi, thickness of 2ndSi.
	 
	 // PosGlobal is the position on the 2nd silicon with (0,0,0) as the original point
	 TVector3 PosGlobal = LocalFrame_Origin[iHira]
	                    + posX*LocalVectorUnit[iHira][0]
	                    + posY*LocalVectorUnit[iHira][1]
	                    + posZ*LocalVectorUnit[iHira][2];
         
         PosGlobal = PosGlobal-BeamSpotPos;
         
         double Theta_Global = PosGlobal.Theta()*RadToDeg(); 
         double Phi_Global = PosGlobal.Phi()*RadToDeg();
         if(Phi_Global<0) { Phi_Global += 360; }
         sprintf(NameTem,"%d   %d   %d   %.4f   %.4f",iHira,iX,iY,Theta_Global,Phi_Global);
         outfile<<NameTem<<endl;
         h2_PixelAngle->Fill(Theta_Global,Phi_Global);
      }
    }
    
  }
  
  
  TCanvas* c1_PixelAngle = new TCanvas(("c1_PixelAngle_"+BeamSpotPos_Tag).c_str(),("c1_PixelAngle_"+BeamSpotPos_Tag).c_str(),1);
  
  if(IsDrawWithSim==1)
  {
    TFile* f1_G4_Sim = new TFile("../GeoEff/f1_GeoEff_V1.root","update");
    //TH2D* h2_Sim_Theta_Phi = (TH2D*) f1_G4_Sim->Get("h2_noBadMap_Theta_Phi_Lab");
    TH2D* h2_Sim_Theta_Phi = (TH2D*) f1_G4_Sim->Get("h2_BadMap_Theta_Phi_Lab");
    h2_Sim_Theta_Phi->Draw("col");
    h2_PixelAngle->Draw("same");
  }
  else
  { h2_PixelAngle->Draw(); }
  
  TFile* f1_PixelAngle = new TFile(("f1_PixelAngle_"+BeamSpotPos_Tag+".root").c_str(),"update");
  f1_PixelAngle->cd();
  h2_PixelAngle->Write("",TObject::kOverwrite);
  f1_PixelAngle->Close();
}
