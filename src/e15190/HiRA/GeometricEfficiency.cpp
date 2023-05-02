#include "GeometricEfficiency.hh"

GeometricEfficiency::GeometricEfficiency()
{
    f1_Results = 0;

    ThetaRange[0] = 20;
    ThetaRange[1] = 80; // unit:deg.
    ThetaBinSize = 0.1; // unit:deg.
}

void GeometricEfficiency::Cal_GeometricEfficiency(TH1D *h1_Count, TH1D *h1_GeometricEfficiency, int TotalSimEvtNum)
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
        double GeometricEfficiency = BinContent / (Phi_Eff * Theta_Eff * BinSize);
        h1_GeometricEfficiency->SetBinContent(i, GeometricEfficiency);
        h1_GeometricEfficiency->SetBinError(i, GeometricEfficiency * BinContent_RelErr); // the relative error is same with the ;
    }
}

void GeometricEfficiency::ReadGeometricEfficiencyHistogram(const std::string &RootFileName)
{
    fs::path path = RootFileName;
    if (!fs::exists(path))
    {
        std::string msg = RootFileName + " does not exists.";
        throw std::invalid_argument(msg.c_str());
    }
    f1_Results = new TFile(RootFileName.c_str(), "READ");

    for (auto &name : this->h1_names)
    {
        this->h1_collection[name] = (TH1D *)f1_Results->Get(name.c_str());
    }
    for (auto &name : this->h2_names)
    {
        this->h2_collection[name] = (TH2D *)f1_Results->Get(name.c_str());
    }
}

double GeometricEfficiency::Get_GeometricEfficiency(const double &Theta_Lab)
{
    int BinIndex = h1_collection["h1_BadMap_Theta_Lab_Eff"]->FindBin(Theta_Lab);
    double Eff = h1_collection["h1_BadMap_Theta_Lab_Eff"]->GetBinContent(BinIndex);

    if (Eff > 0 && Eff < 0.001)
    {
        Eff = 0.001;
    }
    return Eff;
}
