#include "ReactionLost.hh"

ReactionLost::ReactionLost()
{
    Init_Sean();
}

void ReactionLost::Init_Sean()
{
    ReactionLost_Cor_Z = {1, 1, 1, 2, 2};
    ReactionLost_Cor_A = {1, 2, 3, 3, 4};

    std::vector<double> A = {-1.772780E-4, -3.124610E-4, -2.480730E-4, -9.824110E-5, -8.952670E-5};
    std::vector<double> B = {-1.120190E-5, -6.187160E-6, -4.943390E-6, -7.446480E-7, -5.984310E-7};

    ReactionLost_Cor_ParticleName = {"P", "D", "T", "3He", "4He"};
    f1_ReactionLost_CorEff.resize(ReactionLost_Cor_ParticleName.size());

    if (f1_ReactionLost_CorEff.size() > ReactionLost_Cor_ParticleNum)
    {
        throw std::invalid_argument("exceeding maximum number of particle in ReactionLost calculation.");
    }

    for (unsigned int iPID = 0; iPID < f1_ReactionLost_CorEff.size(); iPID++)
    {
        f1_ReactionLost_CorEff[iPID] = new TF1(("f1_ReactionLost_CorEff_" + ReactionLost_Cor_ParticleName[iPID]).c_str(), "TMath::Exp([0]*x*x+[1]*x)", 0, 500);
        f1_ReactionLost_CorEff[iPID]->SetParameter(0, B[iPID]);
        f1_ReactionLost_CorEff[iPID]->SetParameter(1, A[iPID]);
    }
}

double ReactionLost::Get_ReactionLost_CorEff(const int &Z, const int &A, const double &Ekin_Total_Lab)
{
    int Index = -1;
    for (int iPID = 0; iPID < f1_ReactionLost_CorEff.size(); iPID++)
    {
        if (Z == ReactionLost_Cor_Z[iPID] && A == ReactionLost_Cor_A[iPID])
        {
            Index = iPID;
            break;
        }
    }
    if (Index != -1)
    {
        return f1_ReactionLost_CorEff[Index]->Eval(Ekin_Total_Lab);
    }
    else
    {
        return f1_ReactionLost_CorEff[ReactionLost_Cor_ParticleNum - 1]->Eval(Ekin_Total_Lab);
    } // the others are set to same with the last particle( always is alpha ).
}