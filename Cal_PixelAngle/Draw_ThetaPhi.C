void Draw_ThetaPhi()
{
  string Energy = "140";
  char NameTem[200];
  sprintf(NameTem,"pixelAngle_%sMeV.dat",Energy.c_str());
  
  ifstream infile(NameTem);
  infile.ignore(99,'\n');
  
  sprintf(NameTem,"h2_Angle_Pixel_%sMeV",Energy.c_str());
  TH2D* h2_Angle_Pixel = new TH2D(NameTem,";Theta(Deg.);Phi(Deg.)",600,20,80,1000,150,250);
  int TeleIndex; int EF; int EB; double Theta; double Phi;
  while(infile)
  {
    infile>>TeleIndex>>EF>>EB>>Theta>>Phi;
    h2_Angle_Pixel->Fill(Theta,Phi);
  }
  
  TCanvas* c1_Angle_Pixel = new TCanvas("c1_Angle_Pixel","c1_Angle_Pixel",1);
  h2_Angle_Pixel->Draw("colz");
  TFile* f1_HitAngleDis_Hira = new TFile("f1_HitPattern_Hira.root","update");
  h2_Angle_Pixel->Write("",TObject::kOverwrite);
}
