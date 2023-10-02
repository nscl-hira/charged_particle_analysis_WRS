#ifndef ReactionLost_hh
#define ReactionLost_hh

#include <map>
#include <array>
#include <string>
#include "TF1.h"
#include "TMath.h"

class ReactionLost
{
public:
    ReactionLost()
    {
        for (auto &[za, parameters] : reaction_loss_parameters)
        {
            std::string ptcl_name = accepted_particles[za];
            this->f1_ReactionLost_CorEff[ptcl_name] = new TF1(
                ("f1_ReactionLost_CorEff_" + ptcl_name).c_str(),
                "TMath::Exp([0] * x * x + [1] * x)",
                0, 800
                //
            );
            this->f1_ReactionLost_CorEff[ptcl_name]->SetParameter(0, parameters[1]);
            this->f1_ReactionLost_CorEff[ptcl_name]->SetParameter(1, parameters[0]);
        }
    }
    ~ReactionLost()
    {
        for (auto &[_, f1] : this->f1_ReactionLost_CorEff)
        {
            delete f1;
        }
    }

    double Get_ReactionLost_CorEff(const int &z, const int &a, const double &Ekin_Lab)
    {
        // currently, Sean only check the reaction lost correction up to 4He, the heavier particle will be taken as 4He.
        if (this->accepted_particles.count({z, a}) == 0)
        {
            return this->f1_ReactionLost_CorEff["4He"]->Eval(Ekin_Lab);
        }
        std::string name = this->accepted_particles[{z, a}];
        return this->f1_ReactionLost_CorEff[name]->Eval(Ekin_Lab);
    }

private:
    std::map<std::string, TF1 *> f1_ReactionLost_CorEff;
    std::map<std::array<int, 2>, std::array<double, 2>> reaction_loss_parameters = {
        {{1, 1}, {-1.772780E-4, -1.120190E-5}},
        {{1, 2}, {-3.124610E-4, -6.187160E-6}},
        {{1, 3}, {-2.480730E-4, -4.943390E-6}},
        {{2, 3}, {-9.824110E-5, -7.446480E-7}},
        {{2, 4}, {-8.952670E-5, -5.984310E-7}},
    };

    std::map<std::array<int, 2>, std::string> accepted_particles = {
        {{1, 1}, "p"},
        {{1, 2}, "d"},
        {{1, 3}, "t"},
        {{2, 3}, "3He"},
        {{2, 4}, "4He"},
    };
};

#endif
