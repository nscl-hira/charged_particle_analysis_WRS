#include "anal.hh"
#include "HistogramManager.hh"

AME ame;
std::map<int, double> Get_ImpactParameter_Map(const std::string &path_bimp);

int main(int argc, char *argv[])
{
    ArgumentParser argparser(argc, argv);

    auto ReducedImpactParameter = Get_ImpactParameter_Map(argparser.path_bimp);
    Reader *reader = new Reader(
        argparser.reaction,
        argparser.dir_data.string(),
        argparser.path_runinfo.string());

    HistogramManager *hist_raw = new HistogramManager("raw");
    HistogramManager *hist_geoeff = new HistogramManager("geoeff");
    HistogramManager *hist_alleff = new HistogramManager("alleff");

    double beam_mass = ame.GetMass(reader->GetBeamZ(), reader->GetBeamA());
    double target_mass = ame.GetMass(reader->GetTargetZ(), reader->GetTargetA());
    double betacms = Physics::GetReactionBeta(beam_mass, target_mass, reader->GetBeamEnergy(), reader->GetBeamA());
    double rapidity_beam = Physics::GetBeamRapidity(beam_mass, target_mass, reader->GetBeamEnergy(), reader->GetBeamA());

    HiRA *HiRA_detector = new HiRA();
    HiRA_detector->Initialize_LaserAngles(argparser.path_pixel_angles);

    int Runs = reader->GetRuns();

    double Normalization = 0;
    int NumberOfPassedEvents = 0;

    for (int iRun = 0; iRun < Runs; iRun++)
    {
        int RunIndex = reader->mRunIndex[iRun];
        std::string badmap_version = reader->mBadMapVersion[iRun];
        std::string trigger = reader->mTriggerCondition[iRun];

        if (!(trigger == "uBallDS_uBallHira_uBallNW" || trigger == "uBallDS_uBallHira"))
        {
            return 0;
        }

        HiRA_detector->Initialize_StripMap(argparser.dir_geoeff, badmap_version);

        fs::path path_geoeff = argparser.dir_geoeff / ("f1_GeoEff_" + badmap_version + ".root");
        HiRA_detector->Initialize_GeometricEfficiency(path_geoeff);

        reader->Initialize_Chain(iRun);
        for (int ievt = 0; ievt < reader->GetEntries(iRun); ievt++)
        {
            reader->GetEntry(ievt);

            double tdc_trigger_uball_ds = reader->GetTDC_Trig_Uball_DS();
            double tdc_trigger_uball_hira = reader->GetTDC_Trig_Uball_HiRA();
            double tdc_trigger_uball_nw = reader->GetTDC_Trig_Uball_NW();

            int uball_multi = reader->GetUball_Multi();
            int hira_multi = reader->GetHiRA_Multi();

            int *hira_A = reader->GetHiRA_A();
            int *hira_Z = reader->GetHiRA_Z();

            double *hira_pmag = reader->GetHiRA_Pmag();
            int *hira_numcsi = reader->GetHiRA_NumCSI();
            int *hira_numtel = reader->GetHiRA_NumTel();
            int *hira_numstripf = reader->GetHiRA_NumStripF();
            int *hira_numstripb = reader->GetHiRA_NumStripB();

            // event cut here
            if (
                ReducedImpactParameter[uball_multi] > argparser.cut_on_impact_parameter[1] || ReducedImpactParameter[uball_multi] < argparser.cut_on_impact_parameter[0] ||
                uball_multi < argparser.cut_on_multiplicity[0] ||
                uball_multi > argparser.cut_on_multiplicity[1])
            {
                continue;
            }

            if (uball_multi >= argparser.cut_on_multiplicity[0] && tdc_trigger_uball_ds > TDC_TRIGGER_THRESHOLD)
            {
                Normalization += 1.;
            }
            NumberOfPassedEvents += 1;
            for (int ip = 0; ip < hira_multi; ip++)
            {
                if (!HiRA_detector->IsGoodStrip(hira_numtel[ip], hira_numstripf[ip], hira_numstripb[ip]))
                {
                    continue;
                }
                if (hira_A[ip] < 0 || hira_Z[ip] < 0 || hira_A[ip] < hira_Z[ip])
                {
                    continue;
                }

                double thetalab = HiRA_detector->GetTheta(hira_numtel[ip], hira_numstripf[ip], hira_numstripb[ip]) * TMath::DegToRad();

                double phi = HiRA_detector->GetPhi(hira_numtel[ip], hira_numstripf[ip], hira_numstripb[ip]) * TMath::DegToRad();

                double px = hira_pmag[ip] * TMath::Sin(thetalab) * TMath::Cos(phi);
                double py = hira_pmag[ip] * TMath::Sin(thetalab) * TMath::Sin(phi);
                double pz = hira_pmag[ip] * TMath::Cos(thetalab);

                Particle particle = {hira_A[ip] - hira_Z[ip], hira_Z[ip], px, py, pz};
                particle.initialize(betacms, rapidity_beam, ame.GetMass(particle.Z, particle.N + particle.Z));

                if (HiRA_detector->Pass(particle))
                {
                    double EffGeo = HiRA_detector->Get_GeometricEfficiency(particle.theta_lab);
                    double ReactionEff = HiRA_detector->Get_ReactionLost_CorEff(particle.Z, particle.Z + particle.N, particle.kinergy_lab);

                    hist_raw->Fill(particle, 1.);

                    if (EffGeo > 0 && EffGeo < 0.001)
                    {
                        EffGeo = 0.001;
                        hist_geoeff->Fill(particle, 1. / EffGeo);
                    }
                    if (ReactionEff > 0 && ReactionEff < 0.001)
                    {
                        ReactionEff = 0.001;
                        hist_alleff->Fill(particle, 1. / EffGeo / ReactionEff);
                    }
                }
            }
        }
    }

    Normalization *= 300.;
    hist_raw->Normalize(1. / Normalization);
    hist_geoeff->Normalize(1. / Normalization);
    hist_alleff->Normalize(1. / Normalization);

    TFile *file = new TFile(argparser.output_file.c_str(), "RECREATE");
    hist_raw->Write();
    hist_geoeff->Write();
    hist_alleff->Write();
    file->Write();
    file->Close();
}

std::map<int, double> Get_ImpactParameter_Map(const std::string &path_bimp)
{
    std::map<int, double> ReducedImpactParameter;
    std::ifstream infile(path_bimp.c_str());
    int multiplicity;
    double reduced_bimp;
    double buffer;
    infile.ignore(99, '\n');

    while (infile >> multiplicity)
    {
        infile >> reduced_bimp;
        for (int _ = 0; _ < 3; _++)
        {
            infile >> buffer;
        }
        ReducedImpactParameter[multiplicity] = reduced_bimp;
    }
    return ReducedImpactParameter;
}
