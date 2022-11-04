
#ifndef reactionlost_h
#define reactionlost_h

#include "TMath.h"
#include "TF1.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

class reactionlost
{

protected:
    static const unsigned int ReactionLost_Cor_ParticleNum = 5; // currently, Sean only check the reaction lost correction up to 4He, the heavier particle will be taken as 4He.

public:
    reactionlost();
    ~reactionlost(){};
    void Init_Sean();
    double Get_ReactionLost_CorEff(const int &Z, const int &A, const double &Ekin_Total_Lab);

private:
    std::vector<TF1 *> f1_ReactionLost_CorEff;
    std::vector<std::string> ReactionLost_Cor_ParticleName;
    std::vector<int> ReactionLost_Cor_Z;
    std::vector<int> ReactionLost_Cor_A;
};

#endif
