#ifndef HiRA_Particle_h
#define HiRA_Particle_h 1

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "iostream"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "fstream"
#include "TCutG.h"
#include "TVector3.h"
#include "TLorentzVector.h"
using namespace TMath;


class HiRA_Particle : public TObject

{
	private:

	int pid;
	double mass;
	double px,py,pz;
	double E;//, Elab;// E is KE cMS;
	//TLorentzVector* LV;
	int fnumtel, fnumcsi, fnumstripb, fnumstripf;
	//bool IsBadMap;
	int iEvent;
	public:

	HiRA_Particle();
	HiRA_Particle(const HiRA_Particle &p);
	~HiRA_Particle();

	int GetPid();
	double GetMass();
	double GetPx();
	double GetPy();
	double GetPz();
	double GetE();
	//double GetElab();
	int GetTelescope();
	int GetCsI();
	int GetEStripb();
	int GetEStripf();
	//bool GetIsBadMap();
	//TLorentzVector* Get4Vector();
	int GetEventNum();
	
	

	void SetPid(int id);
	void SetMass(double m);
	void SetPx(double p1);
	void SetPy(double p2);
	void SetPz(double p3);
	void SetE(double e);
	//void SetElab(double elab);
	void SetTelescope(int tele);
	void SetCsI(int csi);
	void SetEStripb(int eb);
	void SetEStripf(int ef);
	//void SetIsBadMap(bool yes);
	//void Set4Vector(TLorentzVector* t);
	void SetEventNum(int e);
	ClassDef(HiRA_Particle,1)	
};


#endif








