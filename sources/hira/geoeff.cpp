#include "geoeff.hh"

geoeff::geoeff()
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

    // the below is the default setting of theta range and bin size.
    ThetaRange[0] = 20;
    ThetaRange[1] = 80; // unit:deg.
    ThetaBinSize = 0.1; // unit:deg.
}

void geoeff::Cal_GeoEff(TH1D *h1_Count, TH1D *h1_GeoEff, int TotalSimEvtNum)
{
    int BinNum = h1_Count->GetNbinsX();
    double BinSize = h1_Count->GetBinCenter(2) - h1_Count->GetBinCenter(1);
    double Phi_Eff = 360.0 / (250.0 - 150.0);
    double Theta_Eff = 1.0 * TotalSimEvtNum / (80.0 - 20.0);

    for (int i = 1; i <= BinNum; i++)
    {
        double BinContent = h1_Count->GetBinContent(i);
        double BinContent_RelErr = 0;
        if (BinContent != 0)
        {
            BinContent_RelErr = TMath::Sqrt(BinContent) / BinContent;
        }
        double GeoEff = BinContent / (Phi_Eff * Theta_Eff * BinSize);
        h1_GeoEff->SetBinContent(i, GeoEff);
        h1_GeoEff->SetBinError(i, GeoEff * BinContent_RelErr); // the relative error is same with the ;
    }
}

void geoeff::ReadGeoEffHistogram(const std::string &RootFileName)
{
    f1_Results = new TFile(RootFileName.c_str(), "read");
    h1_BadMap_Theta_Lab_Eff = (TH1D *)f1_Results->Get("h1_BadMap_Theta_Lab_Eff");
    if (h1_BadMap_Theta_Lab_Eff != 0)
    {
        std::cout << "Get the GeoEff Correction histogram!" << std::endl;
    }
    h2_WholeHira_Theta_Phi_Lab = (TH2D *)f1_Results->Get("h2_WholeHira_Theta_Phi_Lab");
    h2_BadMap_Theta_Phi_Lab = (TH2D *)f1_Results->Get("h2_BadMap_Theta_Phi_Lab");
    h2_noBadMap_Theta_Phi_Lab = (TH2D *)f1_Results->Get("h2_noBadMap_Theta_Phi_Lab");
    h1_BadMap_Theta_Lab_HitCount = (TH1D *)f1_Results->Get("h1_BadMap_Theta_Lab_HitCount");
    h1_HiraHit_Multi = (TH1D *)f1_Results->Get("h1_HiraHit_Multi");
}

double geoeff::Get_GeoEff(const double &Theta_Lab)
{
    int BinIndex = h1_BadMap_Theta_Lab_Eff->FindBin(Theta_Lab);
    double Eff = h1_BadMap_Theta_Lab_Eff->GetBinContent(BinIndex);
    return Eff;
}
