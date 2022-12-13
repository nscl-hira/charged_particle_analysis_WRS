#include "data_reduction.hh"
int main(int argc, char **argv)
{
    angles *PositionCalibration = new angles("../database/Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat");
    // std::string path_data = "/data/HiRA_Cali/48Ca64Ni_140MeVu/CalibratedData_4023.root";
    std::string path_data = argv[1];
    std::string path_out = argv[2];

    std::string tr_name = "E15190";

    TChain *in_chain = input_chain(path_data, tr_name);

    TFile *outfile = new TFile(path_out.c_str(), "RECREATE");
    outfile->cd();
    TTree *out_tr = get_tree();

    for (int ievt = 0; ievt < in_chain->GetEntries(); ievt++)
    {
        in_chain->GetEntry(ievt);
        for (int n = 0; n < structure.Hira_multi; n++)
        {
            structure.Hira_theta[n] = PositionCalibration->GetTheta(int(structure.Hira_telescope[n]), int(structure.Hira_stripf[n]), int(structure.Hira_stripb[n]));
            structure.Hira_phi[n] = PositionCalibration->GetPhi(int(structure.Hira_telescope[n]), int(structure.Hira_stripf[n]), int(structure.Hira_stripb[n]));
        }
        out_tr->Fill();
    }
    outfile->Write();
    outfile->Close();

    return 0;
}