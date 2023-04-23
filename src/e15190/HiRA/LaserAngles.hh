#ifndef LaserAngles_hh
#define LaserAngles_hh

#include <iostream>
#include <fstream>
#include <sstream>

#include "TString.h"

class LaserAngles
{
protected:
    static const unsigned int fNumTele = 12;
    static const unsigned int fNumCsI = 4;
    static const unsigned int fNumStripf = 32;
    static const unsigned int fNumStripb = 32;

public:
    LaserAngles(const std::string &path);
    ~LaserAngles(){};
    double ThetaLab[fNumTele][fNumStripf][fNumStripf]; //[HiraIndex][X Index][Y Index]
    double Phi[fNumTele][fNumStripf][fNumStripf];      //[HiraIndex][X Index][Y Index]

    double GetTheta(const int &iH, const int &iX, const int &iY) { return ThetaLab[iH][iX][iY]; }
    double GetPhi(const int &iH, const int &iX, const int &iY) { return Phi[iH][iX][iY]; }
};

#endif
