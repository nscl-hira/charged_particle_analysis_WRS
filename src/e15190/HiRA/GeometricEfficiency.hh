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

class GeometricEfficiency
{
protected:
    std::vector<std::string> h1_names = {"h1_BadMap_Theta_Lab_HitCount", "h1_BadMap_Theta_Lab_Eff", "h1_noBadMap_Theta_Lab_HitCount", "h1_noBadMap_Theta_Lab_Eff", "h1_HiraHit_Multi"};

    std::vector<std::string> h2_names = {"h2_WholeHira_Theta_Phi_Lab", "h2_noBadMap_Theta_Phi_Lab", "h2_BadMap_Theta_Phi_Lab"};

public:
    GeometricEfficiency();
    ~GeometricEfficiency(){};

    void SetTheta(const std::array<double, 2> Range, const double &BinSize)
    {
        ThetaRange[0] = Range[0];
        ThetaRange[1] = Range[1];
        ThetaBinSize = BinSize;
    }

    void Cal_GeometricEfficiency(TH1D *h1_Count, TH1D *h1_GeometricEfficiency, int TotalSimEvtNum);

    // the below is for get the GeometricEfficiency from the Histogram.
    void ReadGeometricEfficiencyHistogram(const std::string &RootFileName);
    double Get_GeometricEfficiency(const double &Theta_Lab); // unit: Deg.
    TH1D *Get_ThetaEff_Histogram() { return this->h1_collection["h1_BadMap_Theta_Lab_Eff"]; }

private:
    TFile *f1_Results;
    std::map<std::string, TH2D *> h2_collection;
    std::map<std::string, TH1D *> h1_collection;
    std::array<double, 2> ThetaRange;
    double ThetaBinSize;
};

#endif
