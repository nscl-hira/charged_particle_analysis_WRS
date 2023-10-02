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
    LaserAngles(const std::string &path)
    {
        std::ifstream infile(path.c_str());
        if (!infile)
        {
            std::cout << "file not found." << std::endl;
            std::exit(1);
        }
        infile.ignore(99, '\n');

        int HiraIndex, StripXIndex, StripYIndex;
        double theta, phi;

        auto strip = [](std::string &line) -> void
        {
            while (line[0] == ' ')
            {
                line = line.substr(1, line.length() - 1);
            }
            while (line[line.length() - 1] == ' ')
            {
                line = line.substr(0, line.length() - 1);
            }
        };

        for (std::string line; std::getline(infile, line);)
        {
            strip(line);
            std::istringstream ss(line);
            ss >> HiraIndex >> StripXIndex >> StripYIndex >> theta >> phi;
            this->ThetaLab[HiraIndex][StripXIndex][StripYIndex] = theta;
            this->Phi[HiraIndex][StripXIndex][StripYIndex] = phi;
        }
    }
    ~LaserAngles(){};
    double ThetaLab[fNumTele][fNumStripf][fNumStripf]; //[HiraIndex][X Index][Y Index]
    double Phi[fNumTele][fNumStripf][fNumStripf];      //[HiraIndex][X Index][Y Index]

    double GetTheta(const int &iH, const int &iX, const int &iY) { return ThetaLab[iH][iX][iY]; }
    double GetPhi(const int &iH, const int &iX, const int &iY) { return Phi[iH][iX][iY]; }
};

#endif
