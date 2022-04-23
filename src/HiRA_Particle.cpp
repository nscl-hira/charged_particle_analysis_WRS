#include "HiRA_Particle.hh"
ClassImp(HiRA_Particle);

HiRA_Particle::HiRA_Particle()
{
 ;  
}

HiRA_Particle::~HiRA_Particle()
{;}

HiRA_Particle::HiRA_Particle(const HiRA_Particle &p)
{
	pid = p.pid;
	mass = p.mass;
	px = p.px;
	py = p.py;
	pz = p.pz;
	E = p.E;
	//Elab = p.Elab;
	fnumtel = p.fnumtel;
	fnumcsi = p.fnumcsi;
	fnumstripb = p.fnumstripb;
	fnumstripf = p.fnumstripf;
	//IsBadMap = p.IsBadMap;
	//LV = p.LV;
	
}


int HiRA_Particle::GetPid() {return pid;}
double HiRA_Particle::GetMass() {return mass;}
double HiRA_Particle::GetPx() {return px;}
double HiRA_Particle::GetPy() {return py;}
double HiRA_Particle::GetPz() {return pz;}
double HiRA_Particle::GetE() {return E;}
//double HiRA_Particle::GetElab() {return Elab;}
int HiRA_Particle::GetTelescope() {return fnumtel;}
int HiRA_Particle::GetCsI() {return fnumcsi;}
int HiRA_Particle::GetEStripb() {return fnumstripb;}
int HiRA_Particle::GetEStripf() {return fnumstripf;}
//bool HiRA_Particle::GetIsBadMap() {return IsBadMap;}
//TLorentzVector* HiRA_Particle::Get4Vector() {return LV;}
int HiRA_Particle::GetEventNum() {return iEvent;}


void HiRA_Particle::SetPid(int id) {pid=id;}
void HiRA_Particle::SetMass(double  m) {mass=m;}
void HiRA_Particle::SetPx(double p1) {px=p1;}
void HiRA_Particle::SetPy(double p2) {py=p2;}
void HiRA_Particle::SetPz(double p3) {pz=p3;}
void HiRA_Particle::SetE(double e) {E=e;}
//void HiRA_Particle::SetElab(double elab) {Elab=elab;}
void HiRA_Particle::SetTelescope(int tele) {fnumtel=tele;}
void HiRA_Particle::SetCsI(int csi) {fnumcsi=csi;}
void HiRA_Particle::SetEStripb(int eb) {fnumstripb=eb;}
void HiRA_Particle::SetEStripf(int ef) {fnumstripf=ef;}
//void HiRA_Particle::SetIsBadMap(bool yes) {IsBadMap=yes;}
//void HiRA_Particle::Set4Vector(TLorentzVector* t) {LV=t;}
void HiRA_Particle::SetEventNum(int e) {iEvent = e;}




