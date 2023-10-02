#include "anal.hh"
#include <argparse/argparse.hpp>
#include <utilities/ame.hpp>
#include <utilities/physics.hpp>
#include <utilities/particle.hpp>
#include <progressbar/progressbar.hpp>
#include "reader.hpp"
#include "hira.hpp"

argparse::ArgumentParser get_arguments(int argc, char *argv[]);
std::map<int, double> Get_ImpactParameter_Map(const std::string &path_bimp);
bool detected(const particle &ptcl, HiRA *&hira_detector);
void fill_histograms(const std::string &mode, const particle &ptcl, const double &ybeam, const double &weight);
void normalize_histograms(const double &norm);
void save_results(const std::string &filename);
int main(int argc, char *argv[])
{
    auto program = get_arguments(argc, argv);
    auto ame_mass_table = get_ame_mass_table();

    std::string reaction = program.get<std::string>("--reaction");
    auto bmulti_map = Get_ImpactParameter_Map((PATH_BMULTI_MAP / Form("%s.dat", reaction.c_str())).string());

    Reader *reader = new Reader(
        reaction,
        DIR_DATA.string(),
        PATH_RUNINFO.string()
        //
    );

    HiRA *hira_detector = new HiRA();
    hira_detector->Initialize_LaserAngles(PATH_PIXEL_ANGLES.string());
    hira_detector->Initialize_ReactionLost();
    define_histograms();

    int beamA, beamZ, targetA, targetZ, beam_energy;
    read_reaction(program.get<std::string>("--reaction"), beamA, beamZ, targetA, targetZ, beam_energy);
    double beam_mass = ame_mass_table[{beamZ, beamA}];
    double target_mass = ame_mass_table[{targetZ, targetA}];
    double betacms = get_reaction_beta(beam_mass, target_mass, beam_energy, beamA);
    double rapidity_beam = get_beam_rapidity(beam_mass, beam_energy, beamA);

    std::cout << "beam mass: " << beam_mass << std::endl;
    std::cout << "target mass: " << target_mass << std::endl;
    std::cout << "beam energy: " << beam_energy << std::endl;
    std::cout << "beta cms: " << betacms << std::endl;
    std::cout << "rapidity beam: " << rapidity_beam << std::endl;

    double norm = 0.;
    long analyzed_nevents = 0;
    bool completed = 0;
    for (int iRun = 0; iRun < reader->GetRuns(); iRun++)
    {
        if (completed)
        {
            break;
        }

        int RunIndex = reader->mRunIndex[iRun];
        if (program.get<std::string>("--mode") == "nruns")
        {
            if (std::find(program.get<std::vector<int>>("--run-indices").begin(), program.get<std::vector<int>>("--run-indices").end(), RunIndex) == program.get<std::vector<int>>("--run-indices").end())
            {
                continue;
            }
        }

        std::string badmap_version = reader->mBadMapVersion[iRun];
        std::string trigger = reader->mTriggerCondition[iRun];

        std::cout << "RunIndex: " << RunIndex << std::endl;
        std::cout << "badmap version: " << badmap_version << std::endl;
        std::cout << "trigger: " << trigger << std::endl;

        if (!(trigger == "uBallDS_uBallHira_uBallNW" || trigger == "uBallDS_uBallHira"))
        {
            continue;
        }

        hira_detector->Initialize_StripMap(DIR_GEOEFF, badmap_version);
        fs::path path_geoeff = DIR_GEOEFF / ("f1_GeoEff_" + badmap_version + ".root");
        hira_detector->Initialize_GeometricEfficiency(path_geoeff);

        reader->Initialize_Chain(iRun);
        ProgressBar pbar(reader->GetEntries(iRun), Form("RunIndex: %i", RunIndex));

        for (int ievt = 0; ievt < reader->GetEntries(iRun); ievt++)
        {
            reader->GetEntry(ievt);
            analyzed_nevents += 1;
            if (program.get<std::string>("--mode") == "nevents" && analyzed_nevents > program.get<long>("--nevents"))
            {
                completed = 1;
                break;
            }
            double tdc_trigger_uball_ds = reader->GetTDC_Trig_Uball_DS();
            double tdc_trigger_uball_hira = reader->GetTDC_Trig_Uball_HiRA();
            double tdc_trigger_uball_nw = reader->GetTDC_Trig_Uball_NW();

            int uball_multi = reader->GetUball_Multi();
            int hira_multi = reader->GetHiRA_Multi();

            int *hira_A = reader->GetHiRA_A();
            int *hira_Z = reader->GetHiRA_Z();
            double *hira_kinergy = reader->GetHiRA_Kinergy();

            int *hira_numcsi = reader->GetHiRA_NumCSI();
            int *hira_numtel = reader->GetHiRA_NumTel();
            int *hira_numstripf = reader->GetHiRA_NumStripF();
            int *hira_numstripb = reader->GetHiRA_NumStripB();

            // cut on reduced-b, multi
            std::vector<int> cut_on_multi = program.get<std::vector<int>>("--cut-on-uball");
            std::vector<double> cut_on_b = program.get<std::vector<double>>("--cut-on-reduced-b");
            if (uball_multi < cut_on_multi[0] || uball_multi > cut_on_multi[1] || bmulti_map[uball_multi] > cut_on_b[1] || bmulti_map[uball_multi] < cut_on_b[0])
            {
                continue;
            }

            if (tdc_trigger_uball_ds > -9990)
            {
                norm += 1.;
            }

            for (int idx = 0; idx < hira_multi; idx++)
            {
                if (!hira_detector->IsGoodStrip(hira_numtel[idx], hira_numstripf[idx], hira_numstripb[idx]))
                {
                    continue;
                }

                if (ame_mass_table.count({hira_Z[idx], hira_A[idx]}) == 0)
                {
                    continue;
                }

                double mass = ame_mass_table[{hira_Z[idx], hira_A[idx]}];
                double theta = hira_detector->GetTheta(hira_numtel[idx], hira_numstripf[idx], hira_numstripb[idx]) * TMath::DegToRad();
                double phi = hira_detector->GetPhi(hira_numtel[idx], hira_numstripf[idx], hira_numstripb[idx]) * TMath::DegToRad();

                double pmag = TMath::Sqrt(pow(hira_kinergy[idx] + mass, 2.) - pow(mass, 2.));
                double px = pmag * TMath::Sin(theta) * TMath::Cos(phi);
                double py = pmag * TMath::Sin(theta) * TMath::Sin(phi);
                double pz = pmag * TMath::Cos(theta);

                particle ptcl = {
                    hira_A[idx] - hira_Z[idx],
                    hira_Z[idx],
                    px / hira_A[idx],
                    py / hira_A[idx],
                    pz / hira_A[idx],
                    mass,
                    "lab"};
                ptcl.initialize(betacms);

                if (detected(ptcl, hira_detector))
                {
                    double geo = hira_detector->Get_GeometricEfficiency(ptcl.theta_lab_deg);
                    double reaction = hira_detector->Get_ReactionLost_CorEff(ptcl.z, ptcl.a, ptcl.kinergy_lab);
                    fill_histograms("raw", ptcl, rapidity_beam, 1.);

                    if (geo > 0.)
                    {
                        fill_histograms("geo", ptcl, rapidity_beam, 1. / geo);
                    }

                    if (reaction > 0. && geo > 0.)
                    {
                        fill_histograms("georeaction", ptcl, rapidity_beam, 1. / geo / reaction);
                    }
                }
            }
            pbar.Update();
        }
        pbar.Finish();
        reader->Clear_Chain(iRun);
    }
    normalize_histograms(norm * 300.);
    save_results(program.get<std::string>("--output-file"));
}

void normalize_histograms(const double &norm)
{
    for (auto &mode : modes)
    {
        for (auto &pn : particles)
        {
            h2_kinergy_theta_lab[mode][pn]->Scale(1. / norm);
            h2_pt_rapidity_lab[mode][pn]->Scale(1. / norm);
        }
    }
    return;
}

void save_results(const std::string &filename)
{
    TFile *output_file = new TFile(filename.c_str(), "RECREATE");
    output_file->cd();
    for (auto &mode : modes)
    {
        for (auto &pn : particles)
        {
            h2_kinergy_theta_lab[mode][pn]->Write();
            h2_pt_rapidity_lab[mode][pn]->Write();
        }
    }
    output_file->Write();
    output_file->Close();
    return;
}

bool detected(const particle &ptcl, HiRA *&hira_detector)
{
    return (
        hira_detector->AcceptTheta({ptcl.z, ptcl.a}, ptcl.theta_lab_deg) &&
        hira_detector->AcceptKinergy({ptcl.z, ptcl.a}, ptcl.kinergy_lab / ptcl.a));
}

void fill_histograms(const std::string &mode, const particle &ptcl, const double &ybeam, const double &weight)
{
    if (nz2particle.count({ptcl.n, ptcl.z}) == 0)
    {
        return;
    }

    std::string name = nz2particle[{ptcl.n, ptcl.z}];
    h2_kinergy_theta_lab[mode][name]->Fill(ptcl.kinergy_lab / ptcl.a, ptcl.theta_lab_deg, weight);
    h2_pt_rapidity_lab[mode][name]->Fill(ptcl.rapidity_lab / ybeam, ptcl.pmag_trans / ptcl.a, weight);
    return;
}

argparse::ArgumentParser get_arguments(int argc, char *argv[])
{
    argparse::ArgumentParser program("analysis", "0.1.0");
    program.add_argument("-r", "--reaction", "reaction name")
        .default_value(std::string{"Ca40Ni58E56"})
        .action([](const std::string &value)
                {
            static const std::vector<std::string> choices = { 
                "Ca40Ni58E56", "Ca40Ni58E140", "Ca48Ni64E56",  "Ca48Ni64E140",
                "Ca40Sn112E56", "Ca40Sn112E140", "Ca48Sn124E56", "Ca48Sn124E140"
            };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "Ca40Ni58E56" }; })
        .required();

    program.add_argument("--mode", "analysis mode (full, nevents, nruns)")
        .default_value(std::string{"full"})
        .action([](const std::string &value)
                {
            static const std::vector<std::string> choices = { "full", "nevents", "nruns"};
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "full" }; });

    program.add_argument("--run-indices", "run indices to be analyzed, only used if `mode` is `nruns`")
        .default_value(std::vector<int>{})
        .nargs('+')
        .scan<'i', int>()
        .required();

    program.add_argument("--nevents", "only analyze `nevents`, only used if `mode` is `nevents`")
        .default_value(DBL_MAX)
        .scan<'i', long>()
        .required();

    program.add_argument("-o", "--output-file", "output file name")
        .default_value(std::string{"output.root"})
        .required();

    program.add_argument("-c", "--cut-on-uball", "cut on uball multiplicity")
        .nargs(2)
        .default_value(std::vector<int>{0, 128})
        .scan<'i', int>()
        .required();

    program.add_argument("-b", "--cut-on-reduced-b", "cut on reduced impact parameter, currently not implemented")
        .nargs(2)
        .default_value(std::vector<double>{0., 0.4})
        .scan<'g', double>()
        .required();

    program.add_argument("-d", "--debug", "debug mode")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--disable-progress-bar", "disable progress bar")
        .default_value(true)
        .implicit_value(false);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    return program;
}

std::map<int, double> Get_ImpactParameter_Map(const std::string &path_bimp)
{
    std::map<int, double> ReducedImpactParameter;
    std::ifstream infile(path_bimp.c_str());
    int multiplicity;
    double reduced_bimp;
    double buffer;
    infile.ignore(99, '\n');

    while (infile >> multiplicity >> reduced_bimp)
    {
        for (int _ = 0; _ < 3; _++)
        {
            infile >> buffer;
        }
        ReducedImpactParameter[multiplicity] = reduced_bimp;
    }

    int MAX_MULTIPLICITY = 100;
    for (int i = multiplicity; i <= MAX_MULTIPLICITY; i++)
    {
        ReducedImpactParameter[i] = reduced_bimp;
    }
    return ReducedImpactParameter;
}
