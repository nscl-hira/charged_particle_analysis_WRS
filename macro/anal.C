#include "TMath.h"
using namespace TMath;


void anal( 
            string OutputFileName = "Ca48Ni64E140_bhat0.0_0.4.root",
            string DataPath = "/data/HiRA_Cali/48Ca64Ni_140MeVu",
            string system = "Ca48Ni64E140", 
            string runinfo = "config/RunInfo.data", 
            string Hira_BadMap_FilePath = "GeoEff", 
            string AnglePath = "Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat",
            double b1= 0.,
            double b2= 0.4
        )
{

    
    string nametag[5]={"p","d","t","3He","4He"};
    string analtag[3]={"","_geoeff","_alleff"};
	//eventcut
	int uBall_fmulti_cut[2]={5,35};
	int Hira_fmulti_cut[2]={1,80};

	// punchthrough 
    double Ekinlabcut[5][2]={{20.0,198.0},{15.0,263.0/2},{12.0,312/3.},{20.,200.},{18.,200.}};
    double thetalabcut[2] = {30.,75.};
    
    
    
//==========================Define histogram======================================//    
    
    TH2D* h2_pta_yybeam_lab[3][5];
    TH2D* h2_ekin_theta_lab[3][5];
    TH2D* h2_ekin_theta_cms[3][5];

    for (int anal=0;anal<3;anal++)
    {
        for (int i=0;i<5;i++)
        {
            h2_pta_yybeam_lab[anal][i] = new TH2D(("h2_pta_yybeam_lab_"+analtag[anal]+"_"+nametag[i]).c_str(),("h2_pta_yybeam_lab_"+analtag[anal]+"_"+nametag[i]).c_str(), 100,0,1,600,0,600);
            
            h2_ekin_theta_lab[anal][i] = new TH2D(("h2_ekin_theta_lab_"+analtag[anal]+"_"+nametag[i]).c_str(),("h2_ekin_theta_lab_"+analtag[anal]+"_"+nametag[i]).c_str(), 600,20,80,300,0,300);
            
            h2_ekin_theta_cms[anal][i] = new TH2D(("h2_ekin_theta_cms_"+analtag[anal]+"_"+nametag[i]).c_str(),("h2_ekin_theta_cms_"+analtag[anal]+"_"+nametag[i]).c_str(), 800,60,140,300,0,300);
            
            h2_pta_yybeam_lab[anal][i]->Sumw2();
            h2_ekin_theta_lab[anal][i]->Sumw2();
            h2_ekin_theta_cms[anal][i]->Sumw2();
        }
    }    
        
    
//==========================finishing definitions======================================//        
    
	Exp_RunInfo* RunInfo = new Exp_RunInfo();
	RunInfo->Read_RunInfo(runinfo);
	int NumberOfRun = RunInfo->Get_RunNum(system);
	cout<<"NumberOfRun= "<<NumberOfRun<<endl;
	double betaCMS = RunInfo->Get_BetaZ_LabToCM(system);
	double ybeam = RunInfo->Get_BeamRapidity_Lab(system);
	cout<<"ybeam lab= "<<ybeam<<endl;
	
    
    HiRA_ReactionLost* Hira_ReactionLost_Corrector = new HiRA_ReactionLost();
	HiRA_PosCali* HiraPos = new HiRA_PosCali();
    HiraPos->Read_Hira_PixelAngle(AnglePath);
	char FileName[70]; 
	
    int TotalPassedEvents=0;
    int TotalPassedTracks=0;
    
    
	for (int iRun=0;iRun<NumberOfRun;iRun++)
	{
		
		int RunIndex = RunInfo->Get_RunNO(system, iRun);
		string BadMapVersion = RunInfo->Get_BadMapVersion(system,iRun);
		string Trigger = RunInfo->Get_Trigger(system,iRun);
		cout<<system<<"  "<<RunIndex<<"  "<<BadMapVersion<<"  "<<Trigger<<endl;
		
		HiRA_BadMap* Hira_BadMapper = new HiRA_BadMap();
		Hira_BadMapper->Set_BadMapper_Version(Hira_BadMap_FilePath, BadMapVersion);
		
		string GeoEffFile = Hira_BadMap_FilePath+"/f1_GeoEff_"+BadMapVersion+".root";
		HiRA_GeoEff* Hira_GeoEfficiency = new HiRA_GeoEff();
		Hira_GeoEfficiency->ReadGeoEffHistogram(GeoEffFile.c_str());
		
		sprintf(FileName,"%s/CalibratedData_%d.root",DataPath.c_str(),RunIndex);
		TChain* t1_Data = new TChain("E15190");
		t1_Data->AddFile(FileName);
		
		//Set Branch address	
		Int_t Hira_fmulti;Int_t uBall_fmulti;
		double uBall_fb; double uBall_fbhat; 
		double fTheta[20]; double fPhi[20]; double fThetaCMS[20];
		int fnumtel[20]; int fnumstripf[20]; int fnumstripb[20]; int fnumcsi[20];
		double fEnergySifCal[20]; double fEnergySibCal[20];
		double fEnergycsifCal[20];
		double fKinEnergy[20];
		double fKinEnergyCMS[20];
		double fMomentum[20];
		double fMomentumCMS[20];
		int fZId[20]; int fAId[20];
		double TDCTriggers_uBallDS_TRG;
		double TDCTriggers_uBallHira_TRG;
		double TDCTriggers_uBallNW_TRG;
		
		//branch address point to number array 
		t1_Data->SetMakeClass(1);
		//set branch not processed
		t1_Data->SetBranchStatus("*", 0);
		t1_Data->SetBranchAddress("uBall.fb", &uBall_fb);
		t1_Data->SetBranchAddress("uBall.fbhat", &uBall_fbhat);
		t1_Data->SetBranchAddress("uBall.fmulti", &uBall_fmulti);
		t1_Data->SetBranchAddress("HiRA.fmulti", &Hira_fmulti);
		t1_Data->SetBranchAddress("HiRA.fnumtel", fnumtel);
		t1_Data->SetBranchAddress("HiRA.fnumcsi", fnumcsi);
		t1_Data->SetBranchAddress("HiRA.fnumstripf", fnumstripf);
		t1_Data->SetBranchAddress("HiRA.fnumstripb", fnumstripb);
		t1_Data->SetBranchAddress("HiRA.fTheta", fTheta);
		t1_Data->SetBranchAddress("HiRA.fPhi", fPhi);
		t1_Data->SetBranchAddress("HiRA.fThetaCMS",fThetaCMS);
		t1_Data->SetBranchAddress("HiRA.fKinEnergy", fKinEnergy);
		t1_Data->SetBranchAddress("HiRA.fKinEnergyCMS",fKinEnergyCMS);
		t1_Data->SetBranchAddress("HiRA.fMomentumCMS",fMomentumCMS);
		t1_Data->SetBranchAddress("HiRA.fMomentum",fMomentum);
		
		t1_Data->SetBranchAddress("HiRA.fZId", fZId);
		t1_Data->SetBranchAddress("HiRA.fAId", fAId);

		t1_Data->SetBranchAddress("TDCTriggers.uBall_DS_TRG", &TDCTriggers_uBallDS_TRG);
		t1_Data->SetBranchAddress("TDCTriggers.uBallHiRA_TRG", &TDCTriggers_uBallHira_TRG);
		t1_Data->SetBranchAddress("TDCTriggers.uBallNW_TRG", &TDCTriggers_uBallNW_TRG);
		
		//in one event:
		int EventNum = t1_Data->GetEntries();	//total num of event in 1 file
		
        
        int pid;
        double mass;
        
        for (int EventIndex=0;EventIndex<EventNum;EventIndex++)
		{
			if (EventIndex%100000==0 && EventIndex>0){cout<<EventIndex<< "events processed..."<<endl;}
			t1_Data->GetEntry(EventIndex);
            
			if (uBall_fmulti>uBall_fmulti_cut[0] && Hira_fmulti>=Hira_fmulti_cut[0] && Hira_fmulti<Hira_fmulti_cut[1] && uBall_fbhat>=b1 && uBall_fbhat<=b2 && TDCTriggers_uBallHira_TRG>-9990 )
            {	
                TotalPassedEvents++;
                for (int n=0;n<Hira_fmulti;n++)
                {
                    mass=100000;
                    pid=99;
                    
                    if (fAId[n]==1 && fZId[n]==1) {pid=0; mass = 938.272;}
                    else if (fAId[n]==2 && fZId[n]==1) {pid=1; mass = 1875.610; }
                    else if (fAId[n]==3 && fZId[n]==1) {pid=2; mass = 2808.9182;}
                    else if (fAId[n]==3 && fZId[n]==2) {pid=3; mass = 2808.3870;}
                    else if (fAId[n]==4 && fZId[n]==2) {pid=4; mass = 3727.374;}    
                    
                    double thetalab = HiraPos->GetTheta(fnumtel[n],fnumstripf[n],fnumstripb[n])*DegToRad();
                    double philab  = HiraPos->GetPhi(fnumtel[n],fnumstripf[n],fnumstripb[n])*DegToRad();
                    //lab to cms
                    double px = fMomentum[n]*TMath::Sin(thetalab)*TMath::Cos(philab);
                    double py = fMomentum[n]*TMath::Sin(thetalab)*TMath::Sin(philab);
                    double pz = fMomentum[n]*TMath::Cos(thetalab);				
                    double pznew = (1/TMath::Sqrt(1-betaCMS*betaCMS)) * (pz-betaCMS*TMath::Sqrt(fMomentum[n]*fMomentum[n] +mass*mass));
                
                    double p=TMath::Sqrt(px*px+py*py+pznew*pznew);
                    double pt= TMath::Sqrt(px*px+py*py);
                    double K=TMath::Sqrt(p*p+mass*mass)-mass;
                    double rap=0.5*TMath::Log((K+mass+pznew)/(K+mass-pznew));
                    double thetacms = TMath::ATan2(pt,pznew)*TMath::RadToDeg();
                    
                    double raplab=rap+0.5*TMath::Log((1+betaCMS)/(1-betaCMS));
                    double raplabnorm=(raplab/ybeam);
                    double ekinlab = TMath::Sqrt(pow(mass,2.)+pow(pz,2.)+pow(pt,2.))- mass;
                    
                    bool IsBad_Hira = Hira_BadMapper->IsBad_Hira(fnumtel[n]);
                    bool IsBad_CsI = Hira_BadMapper->IsBad_CsI(fnumtel[n],fnumstripf[n],fnumstripb[n]);
                    bool IsBad_StripX = Hira_BadMapper->IsBad_StripX(fnumtel[n],fnumstripf[n]);
                    bool IsBad_StripY = Hira_BadMapper->IsBad_StripY(fnumtel[n],fnumstripb[n]);

                    double EffGeo = Hira_GeoEfficiency->Get_GeoEff(thetalab*RadToDeg()); 
                    double ReactionEff = Hira_ReactionLost_Corrector->Get_ReactionLost_CorEff(fZId[n], fAId[n], fKinEnergy[n]);
                
                    bool BadMapCut = (IsBad_CsI ==0 && IsBad_Hira ==0&& IsBad_StripX==0 && IsBad_StripY==0 );
                    
                    bool KineCut = (ekinlab/fAId[n]>=Ekinlabcut[pid][0] && ekinlab/fAId[n]<=Ekinlabcut[pid][1]);
                    
                    bool AngleCut = (thetalab*TMath::RadToDeg() >= thetalabcut[0] && thetalab*TMath::RadToDeg() <= thetalabcut[1]);
                    
                    bool EffCut = (EffGeo>0 && ReactionEff>0);
                    
                    if (BadMapCut && KineCut && AngleCut && EffCut) 
                    {
                        h2_pta_yybeam_lab[0][pid]->Fill(raplabnorm, pt/fAId[n]);
                        h2_pta_yybeam_lab[1][pid]->Fill(raplabnorm, pt/fAId[n],1./ReactionEff);
                        h2_pta_yybeam_lab[2][pid]->Fill(raplabnorm, pt/fAId[n],1./ReactionEff/EffGeo);
                        
                        h2_ekin_theta_lab[0][pid]->Fill(thetalab*TMath::RadToDeg(), fKinEnergy[n]/fAId[n]);
                        h2_ekin_theta_lab[1][pid]->Fill(thetalab*TMath::RadToDeg(), fKinEnergy[n]/fAId[n],1./ReactionEff);
                        h2_ekin_theta_lab[2][pid]->Fill(thetalab*TMath::RadToDeg(), fKinEnergy[n]/fAId[n],1./ReactionEff/EffGeo);
                        
                        h2_ekin_theta_cms[0][pid]->Fill(thetacms, K/fAId[n]);
                        h2_ekin_theta_cms[1][pid]->Fill(thetacms, K/fAId[n],1./ReactionEff);
                        h2_ekin_theta_cms[2][pid]->Fill(thetacms, K/fAId[n],1./ReactionEff/EffGeo);
                        
                    }
                }//each particle				
            }
        }//each event
        t1_Data->Delete();
    }//each run file
    
    
    
    TFile* outfile  = new TFile(OutputFileName.c_str(),"RECREATE");
    
    for (int anal=0;anal<3;anal++)
    {
        for (int i=0;i<5;i++)
        {

            h2_pta_yybeam_lab[anal][i]->Scale(1./TotalPassedEvents);
            h2_ekin_theta_lab[anal][i]->Scale(1./TotalPassedEvents);
            h2_ekin_theta_cms[anal][i]->Scale(1./TotalPassedEvents);
            
            h2_pta_yybeam_lab[anal][i]->Write();
            h2_ekin_theta_lab[anal][i]->Write();
            h2_ekin_theta_cms[anal][i]->Write();    
        }
    }    

    outfile->Write();
    outfile->Write();
    
    
    
    
    
    
    
    
    
    
}
