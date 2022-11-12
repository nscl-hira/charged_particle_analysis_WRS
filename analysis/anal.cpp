#include "anal.hh"

// #define DEBUG

const double uballDS = 300.;
std::string path_runinfo = PROJECT_DIR / "database/config/RunInfo.data";
std::string path_badmap = PROJECT_DIR / "database/GeoEff";
std::string path_angles = PROJECT_DIR / "database/Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat";

struct manager
{
    std::string reaction;
    std::array<double, 2> impact_parameter;
    fs::path dir_data, path_out, path_runinfo, path_badmap, path_angles;

    double betacms, beam_rapidity;
    double normalization = 0.;
    long total_events = 0;
    long passed_events = 0;

    runinfo *m_runinfo;
    eventcut m_eventcut;
    // hira m_hira;
    Hira *m_hira;

    std::map<int, RootReader *> readers;
    std::map<int, int> runIndex;
    std::map<int, std::string> badmaps;
    std::map<int, std::string> triggers;

    histograms hist_raw;
    histograms hist_geoeff;
    histograms hist_alleff;

    void init();
    void run();
    void finish();
};

int main(int argc, char **argv)
{
    std::string reaction = argv[1];
    std::string dir_data = argv[2];
    std::string path_out = argv[3];
    double bmin = std::stod(argv[4]);
    double bmax = std::stod(argv[5]);
    if (argc >= 7)
    {
        path_runinfo = argv[6];
        path_badmap = argv[7];
        path_angles = argv[8];
    }

    std::cout << "reaction : " << reaction << std::endl;
    std::cout << "data directory : " << dir_data << std::endl;
    std::cout << "output path : " << path_out << std::endl;
    std::cout << "impact parameter : " << bmin << "\t" << bmax << std::endl;
    if (argc >= 7)
    {
        std::cout << "path RunInfo : " << path_runinfo << std::endl;
        std::cout << "path Badmap : " << path_badmap << std::endl;
        std::cout << "path hira_angles : " << path_angles << std::endl;
    }

    manager manager = {reaction, {bmin, bmax}, dir_data, path_out, path_runinfo, path_badmap, path_angles};
    manager.init();
    manager.run();
    manager.finish();
}

void manager::init()
{
    this->m_runinfo = new runinfo(this->path_runinfo);
    this->betacms = this->m_runinfo->GetBetaCMS(this->reaction);
    this->beam_rapidity = this->m_runinfo->GetBeamRapidity(this->reaction);

    std::vector<branch> branches = {
        {"TDCTriggers.uBall_DS_TRG", "double"},
        {"TDCTriggers.uBallHiRA_TRG", "double"},
        {"TDCTriggers.uBallNW_TRG", "double"},

        {"uBall.fb", "double"},
        {"uBall.fbhat", "double"},
        {"uBall.fmulti", "int"},
        {"uBall.fnumring", "int[uBall.fmulti]"},
        {"uBall.fnumdet", "int[uBall.fmulti]"},

        {"HiRA.fmulti", "int"},
        {"HiRA.fnumtel", "int[HiRA.fmulti]"},
        {"HiRA.fnumcsi", "int[HiRA.fmulti]"},
        {"HiRA.fnumstripf", "int[HiRA.fmulti]"},
        {"HiRA.fnumstripb", "int[HiRA.fmulti]"},
        {"HiRA.fTheta", "double[HiRA.fmulti]"},
        {"HiRA.fPhi", "double[HiRA.fmulti]"},
        {"HiRA.fThetaCMS", "double[HiRA.fmulti]"},
        {"HiRA.fKinEnergy", "double[HiRA.fmulti]"},
        {"HiRA.fKinEnergyCMS", "double[HiRA.fmulti]"},
        {"HiRA.fMomentum", "double[HiRA.fmulti]"},
        {"HiRA.fMomentumCMS", "double[HiRA.fmulti]"},
        {"HiRA.fZId", "int[HiRA.fmulti]"},
        {"HiRA.fAId", "int[HiRA.fmulti]"},
    };

    for (auto &br : branches)
    {
        br.autofill();
    }

    int NumberOfRun = this->m_runinfo->GetRunNumber(this->reaction);
    for (int iRun = 0; iRun < NumberOfRun; iRun++)
    {
        int RunIndex = this->m_runinfo->GetRunIndex(this->reaction, iRun);
        std::string path_file = Form("%s/CalibratedData_%d.root", this->dir_data.c_str(), RunIndex);
        this->readers[iRun] = new RootReader("E15190");
        this->readers[iRun]->add_file(path_file.c_str());
        this->readers[iRun]->set_branches(branches);

        this->runIndex[iRun] = RunIndex;
        this->badmaps[iRun] = this->m_runinfo->GetBadMapVersion(this->reaction, iRun);
        this->triggers[iRun] = this->m_runinfo->GetTriggerCondition(this->reaction, iRun);
    }

    this->m_hira = new Hira();
    this->m_hira->Init_Angle(this->path_angles);
    // m_hira->Init)

    // m_hira = {this->path_badmap, this->path_angles};
    // m_hira.init();

    this->m_eventcut = {{1, 80}, {5, 80}, this->impact_parameter};
    this->hist_raw = {"raw", this->betacms, this->beam_rapidity};
    this->hist_geoeff = {"geoeff", this->betacms, this->beam_rapidity};
    this->hist_alleff = {"alleff", this->betacms, this->beam_rapidity};

    this->hist_raw.init();
    this->hist_geoeff.init();
    this->hist_alleff.init();
}

void manager::run()
{
    for (auto &[iRun, iReader] : this->readers)
    {
        int RunIndex = this->runIndex[iRun];
        std::string badmap_version = this->badmaps[iRun];
        std::string Trigger = this->triggers[iRun];

        std::cout << this->reaction << "  " << RunIndex << "  " << badmap_version << "  " << Trigger << std::endl;

        if (!(Trigger == "uBallDS_uBallHira_uBallNW" || Trigger == "uBallDS_uBallHira"))
        {
            continue;
        }

        this->m_hira->Init_Badmap(this->path_badmap, badmap_version);

        fs::path path_geoeff = this->path_badmap / ("f1_GeoEff_" + badmap_version + ".root");
        this->m_hira->Init_GeoEff(path_geoeff);

        // this->m_hira.m_badmap = new badmap(this->path_badmap, badmap_version);
        // fs::path path_geoeff = this->path_badmap / ("f1_GeoEff_" + badmap_version + ".root");
        // this->m_hira.m_geoeff = new geoeff();
        // this->m_hira.m_geoeff->ReadGeoEffHistogram(path_geoeff.c_str());

        // total num of event in 1 file
        int EventNum = iReader->tree->GetEntries();
        this->total_events += EventNum;
        std::cout << EventNum << std::endl;
        for (int ievt = 0; ievt < EventNum; ievt++)
        {
            if (ievt % 100000 == 0 && ievt > 0)
            {
                std::cout << ievt << "events processed..." << std::endl;
            }

            double TDCTriggers_uBallDS_TRG, TDCTriggers_uBallHira_TRG, TDCTriggers_uBallNW_TRG;
            int Hira_fmulti, uBall_fmulti;
            double uBall_fbhat;
            int *fnumtel, *fnumstripf, *fnumstripb, *fnumcsi;
            double *fMomentum;
            int *fZId, *fAId;

            std::map<std::string, std::any> map = iReader->get_entry(ievt);
            try
            {
                TDCTriggers_uBallDS_TRG = std::any_cast<double>(map["TDCTriggers.uBall_DS_TRG"]);
                TDCTriggers_uBallHira_TRG = std::any_cast<double>(map["TDCTriggers.uBallHiRA_TRG"]);
                TDCTriggers_uBallNW_TRG = std::any_cast<double>(map["TDCTriggers.uBallNW_TRG"]);

                uBall_fmulti = std::any_cast<int>(map["uBall.fmulti"]);
                uBall_fbhat = std::any_cast<double>(map["uBall.fbhat"]);

                Hira_fmulti = std::any_cast<int>(map["HiRA.fmulti"]);
                fnumtel = std::any_cast<int *>(map["HiRA.fnumtel"]);
                fnumstripf = std::any_cast<int *>(map["HiRA.fnumstripf"]);
                fnumstripb = std::any_cast<int *>(map["HiRA.fnumstripb"]);
                fnumcsi = std::any_cast<int *>(map["HiRA.fnumcsi"]);
                fZId = std::any_cast<int *>(map["HiRA.fZId"]);
                fAId = std::any_cast<int *>(map["HiRA.fAId"]);
                fMomentum = std::any_cast<double *>(map["HiRA.fMomentum"]);
            }
            catch (const std::bad_any_cast &e)
            {
                std::cout << e.what() << '\n';
            }

#ifndef DEBUG
            event event = {Hira_fmulti, uBall_fmulti, uBall_fbhat};
            if (!this->m_eventcut.pass(event))
            {
                continue;
            }
            passed_events++;
            if (TDCTriggers_uBallDS_TRG > -9990)
            {
                this->normalization += 1.;
            }

            for (int n = 0; n < Hira_fmulti; n++)
            {
                // double thetalab = this->m_hira.m_angles->GetTheta(fnumtel[n], fnumstripf[n], fnumstripb[n]) * TMath::DegToRad();
                // double phi = this->m_hira.m_angles->GetPhi(fnumtel[n], fnumstripf[n], fnumstripb[n]) * TMath::DegToRad();

                double thetalab = this->m_hira->GetTheta(fnumtel[n], fnumstripf[n], fnumstripb[n]) * TMath::DegToRad();

                double phi = this->m_hira->GetPhi(fnumtel[n], fnumstripf[n], fnumstripb[n]) * TMath::DegToRad();

                double px = fMomentum[n] * TMath::Sin(thetalab) * TMath::Cos(phi);
                double py = fMomentum[n] * TMath::Sin(thetalab) * TMath::Sin(phi);
                double pz = fMomentum[n] * TMath::Cos(thetalab);

                particle particle = {fAId[n] - fZId[n], fZId[n], px, py, pz};
                particle.init(this->betacms);

                if (this->m_hira->pass(particle))
                {

                    if (!this->m_hira->pass_badmap(fnumtel[n], fnumstripf[n], fnumstripb[n]))
                    {
                        continue;
                    }

                    // double EffGeo = this->m_hira.m_geoeff->Get_GeoEff(particle.thetalab);
                    // double ReactionEff = this->m_hira.m_reactionlost->Get_ReactionLost_CorEff(particle.zid, particle.aid, particle.ekinlab);

                    double EffGeo = this->m_hira->Get_GeoEff(particle.thetalab);
                    double ReactionEff = this->m_hira->Get_ReactionLost_CorEff(particle.zid, particle.aid, particle.ekinlab);

                    if (EffGeo > 0 && EffGeo < 0.001)
                    {
                        EffGeo = 0.001;
                    }
                    if (ReactionEff > 0 && ReactionEff < 0.001)
                    {
                        ReactionEff = 0.001;
                    }

                    this->hist_raw.fill(particle);
                    if (EffGeo > 0)
                    {
                        this->hist_geoeff.fill(particle, 1. / EffGeo);
                    }
                    if (EffGeo > 0 && ReactionEff > 0)
                    {
                        this->hist_alleff.fill(particle, 1. / EffGeo / ReactionEff);
                    }
                }
            }
#endif
        }
        // this->m_hira.reset();
        delete iReader;
    }
}

void manager::finish()
{
    TFile *outf = new TFile(this->path_out.c_str(), "RECREATE");
    this->normalization *= uballDS;
    this->hist_raw.normalize(1. / this->normalization);
    this->hist_geoeff.normalize(1. / this->normalization);
    this->hist_alleff.normalize(1. / this->normalization);

    this->hist_raw.write();
    this->hist_geoeff.write();
    this->hist_alleff.write();
    outf->Write();
    outf->Close();

    std::cout << this->passed_events << " / " << this->total_events << std::endl;
}

// void get_path_io(string &path_in, string &path_out, const string &this->reaction, const string &dir_data, const string &dir_out, const double &b1, const double &b2)
// {

//     std::vector<std::string> res;
//     std::regex regex_nuc("[A-Z][a-z][0-9]+");
//     std::regex regex_en("E[0-9]+");

//     sregex_iterator current_match(this->reaction.begin(), this->reaction.end(), regex_nuc);
//     sregex_iterator last_match;

//     while (current_match != last_match)
//     {
//         smatch match = *current_match;
//         res.push_back(string(match.str()));
//         current_match++;
//     }

//     std::string beam = res[0];
//     std::string target = res[1];

//     std::smatch match;
//     std::regex_search(this->reaction, match, regex_en);
//     std::string energy = match[0];

//     std::string beam_nuc = beam.substr(0, 2);
//     std::string target_nuc = target.substr(0, 2);
//     int beam_a = stoi(beam.substr(2, beam.length() - 2));
//     int target_a = stoi(target.substr(2, target.length() - 2));
//     int beam_e = stoi(energy.substr(1, energy.length() - 1));

//     path_in = Form("%s/%i%s%i%s_%iMeVu", dir_data.c_str(), beam_a, beam_nuc.c_str(), target_a, target_nuc.c_str(), beam_e);

//     path_out = Form("%s/%s_bhatmin%.1f_bhatmax%.1f.root", dir_out.c_str(), this->reaction.c_str(), b1, b2);
//     return;
// }
