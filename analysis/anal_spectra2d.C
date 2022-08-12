

struct particle {
    int aid, zid;
    double p, thetalab, phi;
    double mass, pcms, pzcms, rapiditylab, rapiditycms, thetacms, ekinlab, ekincms, pt; 
    int id; 
    std::string name;
    void autofill(const double& betacms);
};

double get_mass(const int& aid, const int& zid);
int get_pid(const int& aid, const int& zid);
void get_path_io(string& path_in, string& path_out, const string& system, const string& dir_data, const string& dir_out, const double& b1, const double& b2);

void anal_spectra2d( 
    std::string system = "Ca48Ni64E140",
    double b1= 0.,
    double b2= 0.4,
    std::string dir_out = "../output/spectra2d",
    std::string dir_data = "/data/HiRA_Cali",
    std::string path_runinfo = "../database/config/RunInfo.data", 
    std::string path_badmap = "../database/GeoEff", 
    std::string path_angles = "../database/Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat"
)
{
    std::string path_data, path_out;
    get_path_io(path_data, path_out, system, dir_data, dir_out, b1, b2);
    
    std::vector<std::string> particlename = { "p", "d", "t", "3He", "4He" };
    std::vector<std::string> analtag = { "raw", "geoeff", "alleff" };
    
	//eventcut
	std::array<int, 2> uBall_fmulti_cut = {5, 80};
	std::array<int, 2> Hira_fmulti_cut = {1, 80};

    std::vector<std::array<double, 2>> Ekinlabcut = {
        {20.0, 198.0}, {15.0, 263.0/2}, {12.0, 312/3.}, {20., 200.}, {18., 200.}
    };
    std::vector<std::array<double, 2>> thetalabcut = {
        {30.,75.}, {30.,75.}, {30.,75.}, {30.,75.}, {30.,75.},
    };
    
//==========================Define histogram======================================//    
    
    TH2D* h2_pta_yybeam_lab[analtag.size()][particlename.size()];
    TH2D* h2_ekin_theta_lab[analtag.size()][particlename.size()];
    TH2D* h2_ekin_theta_cms[analtag.size()][particlename.size()];
    
    for (unsigned int anal = 0; anal < analtag.size(); anal++) {
        for (unsigned int i = 0; i < particlename.size(); i++) {
            h2_pta_yybeam_lab[anal][i] = new TH2D(("h2_pta_yybeam_lab_"+analtag[anal]+"_"+particlename[i]).c_str(),("h2_pta_yybeam_lab_"+analtag[anal]+"_"+particlename[i]).c_str(), 100, 0, 1, 600, 0, 600);
            
            h2_ekin_theta_lab[anal][i] = new TH2D(("h2_ekin_theta_lab_"+analtag[anal]+"_"+particlename[i]).c_str(),("h2_ekin_theta_lab_"+analtag[anal]+"_"+particlename[i]).c_str(), 600, 20, 80, 300, 0, 300);
            
            h2_ekin_theta_cms[anal][i] = new TH2D(("h2_ekin_theta_cms_"+analtag[anal]+"_"+particlename[i]).c_str(),("h2_ekin_theta_cms_"+analtag[anal]+"_"+particlename[i]).c_str(), 1200, 20, 140, 300, 0, 300);
            
            
            h2_pta_yybeam_lab[anal][i]->Sumw2();
            h2_ekin_theta_lab[anal][i]->Sumw2();
            h2_ekin_theta_cms[anal][i]->Sumw2();
        }
    }    
        
    
//==========================finishing definitions======================================//        
    
	Exp_RunInfo* RunInfo = new Exp_RunInfo();
	RunInfo->Read_RunInfo(path_runinfo);
    
	int NumberOfRun = RunInfo->Get_RunNum(system);
	double betaCMS = RunInfo->Get_BetaZ_LabToCM(system);
	double ybeam = RunInfo->Get_BeamRapidity_Lab(system);
    
    std::cout << "NumberOfRun= " << NumberOfRun << std::endl;
    std::cout << "betaCMS= " << betaCMS << std::endl;
	std::cout << "ybeam lab= " << ybeam << std::endl;
	
    HiRA_ReactionLost* Hira_ReactionLost_Corrector = new HiRA_ReactionLost();
	HiRA_PosCali* HiraPos = new HiRA_PosCali();
    HiraPos->Read_Hira_PixelAngle(path_angles);
	
    int total_passed_events_uball = 0;
    int total_passed_events_hira = 0;
    
    int uball_triggered_event = 0;
    
    
	for (int iRun = 0; iRun < NumberOfRun; iRun++) {
		
		int RunIndex = RunInfo->Get_RunNO(system, iRun);
		string BadMapVersion = RunInfo->Get_BadMapVersion(system,iRun);
		string Trigger = RunInfo->Get_Trigger(system,iRun);
		std::cout << system << "  " << RunIndex << "  " << BadMapVersion << "  " << Trigger << std::endl;
        
        if(!(Trigger=="uBallDS_uBallHira_uBallNW" || Trigger=="uBallDS_uBallHira")) { continue; }
        
		HiRA_BadMap* Hira_BadMapper = new HiRA_BadMap();
		Hira_BadMapper->Set_BadMapper_Version(path_badmap, BadMapVersion);
		
		std::string GeoEffFile = path_badmap + "/f1_GeoEff_" + BadMapVersion + ".root";
		HiRA_GeoEff* Hira_GeoEfficiency = new HiRA_GeoEff();
		Hira_GeoEfficiency->ReadGeoEffHistogram(GeoEffFile.c_str());
        
        std::string path_file = Form("%s/CalibratedData_%d.root", path_data.c_str(), RunIndex);
		TChain* t1_Data = new TChain("E15190");
		t1_Data->AddFile(path_file.c_str());
		
		//Set Branch address	
		Int_t Hira_fmulti;Int_t uBall_fmulti;
		double uBall_fb; double uBall_fbhat; 
        int uBall_fnumring[100];
        int uBall_fnumdet[100];
		double fTheta[50]; double fPhi[50]; double fThetaCMS[50];
		int fnumtel[50]; int fnumstripf[50]; int fnumstripb[50]; int fnumcsi[50];
		double fEnergySifCal[50]; double fEnergySibCal[50];
		double fEnergycsifCal[50];
		double fKinEnergy[50];
		double fKinEnergyCMS[50];
		double fMomentum[50];
		double fMomentumCMS[50];
		int fZId[50]; int fAId[50];
		double TDCTriggers_uBallDS_TRG;
		double TDCTriggers_uBallHira_TRG;
		double TDCTriggers_uBallNW_TRG;
		
		//branch address point to number array 
		t1_Data->SetMakeClass(1);
		t1_Data->SetBranchStatus("*", 0);
		t1_Data->SetBranchAddress("uBall.fb", &uBall_fb);
		t1_Data->SetBranchAddress("uBall.fbhat", &uBall_fbhat);
		t1_Data->SetBranchAddress("uBall.fmulti", &uBall_fmulti);
   		t1_Data->SetBranchAddress("uBall.fnumring", &uBall_fnumring);
        t1_Data->SetBranchAddress("uBall.fnumdet", &uBall_fnumdet);
		t1_Data->SetBranchAddress("HiRA.fmulti", &Hira_fmulti);
		t1_Data->SetBranchAddress("HiRA.fnumtel", fnumtel);
		t1_Data->SetBranchAddress("HiRA.fnumcsi", fnumcsi);
		t1_Data->SetBranchAddress("HiRA.fnumstripf", fnumstripf);
		t1_Data->SetBranchAddress("HiRA.fnumstripb", fnumstripb);
		t1_Data->SetBranchAddress("HiRA.fTheta", fTheta);
		t1_Data->SetBranchAddress("HiRA.fPhi", fPhi);
		t1_Data->SetBranchAddress("HiRA.fThetaCMS", fThetaCMS);
		t1_Data->SetBranchAddress("HiRA.fKinEnergy", fKinEnergy);
		t1_Data->SetBranchAddress("HiRA.fKinEnergyCMS", fKinEnergyCMS);
		t1_Data->SetBranchAddress("HiRA.fMomentumCMS", fMomentumCMS);
		t1_Data->SetBranchAddress("HiRA.fMomentum", fMomentum);
		t1_Data->SetBranchAddress("HiRA.fZId", fZId);
		t1_Data->SetBranchAddress("HiRA.fAId", fAId);
		t1_Data->SetBranchAddress("TDCTriggers.uBall_DS_TRG", &TDCTriggers_uBallDS_TRG);
		t1_Data->SetBranchAddress("TDCTriggers.uBallHiRA_TRG", &TDCTriggers_uBallHira_TRG);
		t1_Data->SetBranchAddress("TDCTriggers.uBallNW_TRG", &TDCTriggers_uBallNW_TRG);
		
		//in one event:
		int EventNum = t1_Data->GetEntries();	//total num of event in 1 file
        
        for (int EventIndex = 0; EventIndex < EventNum; EventIndex++) {
			
            if (EventIndex%100000==0 && EventIndex>0){
                std::cout <<EventIndex << "events processed..." << std::endl;
            }
            
			t1_Data->GetEntry(EventIndex);
            
            if (uBall_fmulti >= uBall_fmulti_cut[0] && uBall_fbhat >= b1 && uBall_fbhat < b2) {
                if (TDCTriggers_uBallDS_TRG > -9990) {
                    uball_triggered_event ++;
                }
                
                total_passed_events_uball++;
                if (Hira_fmulti >= Hira_fmulti_cut[0] && Hira_fmulti < Hira_fmulti_cut[1]) {
                    total_passed_events_hira++;
                }
            }
            
			if (uBall_fmulti >= uBall_fmulti_cut[0] && Hira_fmulti >= Hira_fmulti_cut[0] && Hira_fmulti < Hira_fmulti_cut[1] && uBall_fbhat >= b1 && uBall_fbhat < b2) {

                for (int n = 0; n < Hira_fmulti; n++) {
                    
                    double thetalab = HiraPos->GetTheta(fnumtel[n], fnumstripf[n], fnumstripb[n])* DegToRad();
                    double philab  = HiraPos->GetPhi(fnumtel[n], fnumstripf[n], fnumstripb[n])* DegToRad();
                    
                    particle particle = {fAId[n], fZId[n], fMomentum[n], thetalab, philab};
                    particle.autofill(betaCMS);
                    
                    if (particle.id >= 5){continue;}

                    bool IsBad_Hira = Hira_BadMapper->IsBad_Hira(fnumtel[n]);
                    bool IsBad_CsI = Hira_BadMapper->IsBad_CsI(fnumtel[n], fnumstripf[n], fnumstripb[n]);
                    bool IsBad_StripX = Hira_BadMapper->IsBad_StripX(fnumtel[n], fnumstripf[n]);
                    bool IsBad_StripY = Hira_BadMapper->IsBad_StripY(fnumtel[n], fnumstripb[n]);

                    double EffGeo = Hira_GeoEfficiency->Get_GeoEff(particle.thetalab*RadToDeg()); 
                    double ReactionEff = Hira_ReactionLost_Corrector->Get_ReactionLost_CorEff(particle.zid, particle.aid, particle.ekinlab);
                    
                    if (EffGeo > 0 && EffGeo < 0.001){EffGeo = 0.001;}
                    if (ReactionEff > 0 && ReactionEff < 0.001){ReactionEff = 0.001;}
                    
                    bool BadMapCut = (IsBad_CsI == 0 && IsBad_Hira == 0 && IsBad_StripX == 0 && IsBad_StripY == 0);
                    
                    bool KineCut = (particle.ekinlab/particle.aid >= Ekinlabcut[particle.id][0] && particle.ekinlab/particle.aid <= Ekinlabcut[particle.id][1]);
                    
                    bool AngleCut = (particle.thetalab *TMath::RadToDeg() >= thetalabcut[particle.id][0] && particle.thetalab * TMath::RadToDeg() <= thetalabcut[particle.id][1]);
                    
                    bool EffCut = (EffGeo > 0 && ReactionEff > 0);
                    bool GeoEffCut = (EffGeo > 0);
                    bool ReactionEffCut = (ReactionEff > 0);
                    
                    if (BadMapCut && KineCut && AngleCut) {
                                
                        h2_pta_yybeam_lab[0][particle.id]->Fill(particle.rapiditylab/ybeam, particle.pt/particle.aid);
                        
                        h2_ekin_theta_lab[0][particle.id]->Fill(particle.thetalab * TMath::RadToDeg(), particle.ekinlab/particle.aid);
                        
                        h2_ekin_theta_cms[0][particle.id]->Fill(particle.thetacms * TMath::RadToDeg(), particle.ekincms/particle.aid);
                        
                        if (GeoEffCut) {
                            h2_pta_yybeam_lab[1][particle.id]->Fill(particle.rapiditylab/ybeam, particle.pt/particle.aid, 1./EffGeo);
                            
                            h2_ekin_theta_lab[1][particle.id]->Fill(particle.thetalab * TMath::RadToDeg(), particle.ekinlab/particle.aid, 1./EffGeo);
                        
                            h2_ekin_theta_cms[1][particle.id]->Fill(particle.thetacms * TMath::RadToDeg(), particle.ekincms/particle.aid, 1./EffGeo);
                        }
                        if (EffCut) {
                            h2_pta_yybeam_lab[2][particle.id]->Fill(particle.rapiditylab/ybeam, particle.pt/particle.aid, 1./EffGeo/ReactionEff);
                            
                            h2_ekin_theta_lab[2][particle.id]->Fill(particle.thetalab * TMath::RadToDeg(), particle.ekinlab/particle.aid, 1./EffGeo/ReactionEff);
                        
                            h2_ekin_theta_cms[2][particle.id]->Fill(particle.thetacms * TMath::RadToDeg(), particle.ekincms/particle.aid, 1./EffGeo/ReactionEff);
                        }
                    }
                }//each particle				
            }
        }//each event
        t1_Data->Delete();
    }//each run file
    
    TFile* outfile  = new TFile(path_out.c_str(),"RECREATE");
    
    TNamed* tn_uball = new TNamed("total_passed_events_uball", Form("%i", total_passed_events_uball));
    TNamed* tn_uball_triggered_event = new TNamed("uball_triggered_event", Form("%i", uball_triggered_event));
    TNamed* tn_hira = new TNamed("total_passed_events_hira", Form("%i", total_passed_events_hira));
    tn_uball_triggered_event->Write();
    tn_uball->Write();
    tn_hira->Write();
    
    for (unsigned int anal = 0; anal < analtag.size(); anal++) {
        for (unsigned int i = 0; i < particlename.size(); i++) {
            h2_pta_yybeam_lab[anal][i]->Write();
            h2_ekin_theta_lab[anal][i]->Write();
            h2_ekin_theta_cms[anal][i]->Write();    
        }
    }    

    
    outfile->Write();
    outfile->Write();
    
    
}

int get_pid(const int& aid, const int& zid) {
    std::vector<std::array<int,2>> array_az = {
        {1, 1}, {2, 1}, {3, 1}, {3, 2}, {4, 2}
    };
    
    std::array<int, 2> temp = {aid, zid};
    int pid = std::find(array_az.begin(), array_az.end(), temp) - array_az.begin();
    
    
    if (aid == 1 && zid == 1) {pid = 0;}
    else if (aid == 2 && zid == 1) {pid = 1;}
    else if (aid == 3 && zid == 1) {pid = 2;}
    else if (aid == 3 && zid == 2) {pid = 3;}
    else if (aid == 4 && zid == 2) {pid = 4;}
    else {pid = 5;}
    
    return pid;
}

double get_mass(const int& aid, const int& zid) {
    
    std::vector<double> masses = {
//         938.272, 1875.610, 2808.9182, 2808.3870, 3727.374
        938.272, 1876.179, 2809.514,  2809.496, 3728.510 //Rensheng numbers
    };
    
    std::vector<std::array<int,2>> array_az = {
        {1, 1}, {2, 1}, {3, 1}, {3, 2}, {4, 2}
    };
    
    std::array<int, 2> temp = {aid, zid};
    int pid = std::find(array_az.begin(), array_az.end(), temp) - array_az.begin();
    double mass = (pid < 5) ? masses[pid] : 100000 ;
    return mass;
}


void particle::autofill(const double& betacms) {
    
    this->id = get_pid(this->aid, this->zid);
    this->mass = get_mass(this->aid, this->zid);
    
    double px = this->p * TMath::Sin(this->thetalab) * TMath::Cos(this->phi);
    double py = this->p * TMath::Sin(this->thetalab) * TMath::Sin(this->phi);
    double pz = this->p * TMath::Cos(this->thetalab);
    
    this->pzcms = (1./TMath::Sqrt(1-betacms*betacms)) * (pz -betacms * TMath::Sqrt(pow(this->p, 2.) + pow(this->mass, 2.)));

    this->pt = TMath::Sqrt(pow(px, 2.) + pow(py, 2.));
    this->pcms = TMath::Sqrt(pow(this->pt, 2.) + pow(this->pzcms, 2.));
    
    this->ekinlab = TMath::Sqrt(pow(this->p, 2.) + pow(this->mass, 2.)) - this->mass;
    this->ekincms = TMath::Sqrt(pow(this->pcms, 2.) + pow(this->mass, 2.)) - this->mass;
        
    this->rapiditycms = 0.5*TMath::Log((this->ekincms + this->mass + this->pzcms) / (this->ekincms + this->mass - this->pzcms));
    
    this->rapiditylab = 0.5*TMath::Log((this->ekinlab + this->mass + pz) / (this->ekinlab + this->mass - pz));

    this->thetacms = TMath::ATan2(this->pt, this->pzcms);
    
}




void get_path_io(string& path_in, string& path_out, const string& system, const string& dir_data, const string& dir_out, const double& b1, const double& b2){
    
    std::vector<std::string> res;
    std::regex regex_nuc("[A-Z][a-z][0-9]+");
    std::regex regex_en("E[0-9]+");
    
    sregex_iterator current_match(system.begin(), system.end(), regex_nuc);
    sregex_iterator last_match;
    
    while (current_match != last_match){
        smatch match = *current_match;
        res.push_back(string(match.str()));
        current_match++;
    }
    
    std::string beam = res[0];
    std::string target = res[1];
    
    std::smatch match;
    std::regex_search(system, match, regex_en);
    std::string energy = match[0];
    
    std::string beam_nuc = beam.substr(0, 2);
    std::string target_nuc = target.substr(0, 2);
    int beam_a = stoi(beam.substr(2, beam.length() - 2));
    int target_a = stoi(target.substr(2, target.length() - 2));
    int beam_e = stoi(energy.substr(1, energy.length() - 1));
 
    
    path_in = Form("%s/%i%s%i%s_%iMeVu", dir_data.c_str(), beam_a, beam_nuc.c_str(), target_a, target_nuc.c_str(), beam_e);
    
    path_out = Form("%s/%s_bhatmin%.1f_bhatmax%.1f.root", dir_out.c_str(), system.c_str(), b1, b2);
    return;
}







