#ifndef GeometricEfficiency_hh
#define GeometricEfficiency_hh

#include <iostream>
#include <fstream>
#include <array>
#include <filesystem>
namespace fs = std::filesystem;

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"

class GeometricEfficiency
{
protected:
    std::vector<std::string> h1_names = {"h1_BadMap_Theta_Lab_HitCount", "h1_BadMap_Theta_Lab_Eff", "h1_noBadMap_Theta_Lab_HitCount", "h1_noBadMap_Theta_Lab_Eff", "h1_HiraHit_Multi"};

    std::vector<std::string> h2_names = {"h2_WholeHira_Theta_Phi_Lab", "h2_noBadMap_Theta_Phi_Lab", "h2_BadMap_Theta_Phi_Lab"};

public:
    GeometricEfficiency()
    {
        ThetaRange[0] = 20;
        ThetaRange[1] = 80; // unit:deg.
        ThetaBinSize = 0.1; // unit:deg.
    }
    ~GeometricEfficiency(){};

    void SetTheta(const std::array<double, 2> Range, const double &BinSize)
    {
        ThetaRange[0] = Range[0];
        ThetaRange[1] = Range[1];
        ThetaBinSize = BinSize;
    }

    void Cal_GeometricEfficiency(TH1D *h1_Count, TH1D *h1_GeometricEfficiency, int TotalSimEvtNum)
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

    // the below is for get the GeometricEfficiency from the Histogram.
    void ReadGeometricEfficiencyHistogram(const std::string &RootFileName)
    {
        fs::path path = RootFileName;
        if (!fs::exists(path))
        {
            std::string msg = RootFileName + " does not exists.";
            throw std::invalid_argument(msg.c_str());
            std::exit(1);
        }
        TFile* f1_Results = new TFile(RootFileName.c_str(), "READ");

        std::vector<std::string> containing_histograms;
        TList* list_of_keys = f1_Results->GetListOfKeys();
        for (TObject* key : *list_of_keys){
            std::string name_str = key->GetName();
            containing_histograms.push_back(name_str);
        }

        for (auto &name : this->h1_names)
        {
            if (std::find(containing_histograms.begin(), containing_histograms.end(), name) == containing_histograms.end()) {
                std::string msg = name + " does not exists in " + RootFileName;
                std::cerr << msg << std::endl;
                continue;
            }
            this->h1_collection[name] = (TH1D *)f1_Results->Get(name.c_str());
            this->h1_collection[name]->SetDirectory(0);
        }
        for (auto &name : this->h2_names)
        {
            if (std::find(containing_histograms.begin(), containing_histograms.end(), name) == containing_histograms.end()) {
                std::string msg = name + " does not exists in " + RootFileName;
                std::cerr << msg << std::endl;
                continue;
            }
            this->h2_collection[name] = (TH2D *)f1_Results->Get(name.c_str());
            this->h2_collection[name]->SetDirectory(0);
        }
    }

    double Get_GeometricEfficiency(const double &theta_lab_deg)
    {
        int BinIndex = h1_collection["h1_BadMap_Theta_Lab_Eff"]->FindBin(theta_lab_deg);
        double Eff = h1_collection["h1_BadMap_Theta_Lab_Eff"]->GetBinContent(BinIndex);

        if (Eff > 0 && Eff < 0.001)
        {
            Eff = 0.001;
        }
        return Eff;
    }

    std::map<std::string, TH2D *> h2_collection;
    std::map<std::string, TH1D *> h1_collection;
    std::array<double, 2> ThetaRange;
    double ThetaBinSize;
};

#endif
