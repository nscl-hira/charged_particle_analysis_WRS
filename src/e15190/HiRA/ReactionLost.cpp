#include "ReactionLost.hh"

E15190ReactionLost::E15190ReactionLost(const std::string &pth, const std::string &fcn_name)
{
    this->function_name = fcn_name;

    // get acceptable ame name
    this->_define_acceptable_particles();

    // set default values
    this->heavy_particle_substitute = this->default_heavy_particle_substitute;
    this->efficiency_range = this->default_efficiency_range;
    this->efficiency_function = this->default_efficiency_function;

    // read parameters
    this->_inititalize(pth);
}
void E15190ReactionLost::_inititalize(const std::string &pth)
{
    json reaction_loss_parameters;
    std::ifstream json_file(pth.c_str());
    json_file >> reaction_loss_parameters;

    std::string eff_fcn = (efficiency_function == "") ? this->default_efficiency_function : efficiency_function;

    // set up correction function for each acceptable particle
    for (auto &[name, parameters] : reaction_loss_parameters.items())
    {
        std::string fcn_name = this->function_name + "_" + name;

        // get ame chemical symbol from alias
        auto optional_symbol = AME::get_instance()->GetSymbol(name);
        std::string ame_name = (optional_symbol) ? optional_symbol.value() : "not found";

        // pass if the ame name is not in the acceptable list
        auto accepted = [this](const std::string &name)
        {
            auto iter = std::find(this->acceptable_particles.begin(), this->acceptable_particles.end(), name);
            return (iter != this->acceptable_particles.end());
        };

        if (!accepted(ame_name))
        {
            continue;
        }

        this->f1_ReactionLost_CorEff[ame_name] = new TF1(
            fcn_name.c_str(),
            eff_fcn.c_str(),
            this->efficiency_range[0],
            this->efficiency_range[1]
            //
        );
        this->f1_ReactionLost_CorEff[ame_name]->SetParameter(0, parameters["b"].get<double>());
        this->f1_ReactionLost_CorEff[ame_name]->SetParameter(1, parameters["a"].get<double>());
    }

    json_file.close();
}

double E15190ReactionLost::Get_ReactionLost_CorEff(const int &z, const int &a, const double &Ekin_Lab)
{
    auto optional_symbol = AME::get_instance()->GetSymbol(z, a);
    if (!optional_symbol)
    {
        return this->unphysical_effeciency;
    }

    std::string ame_name = optional_symbol.value();
    if (!this->_accepted(ame_name))
    {
        // use substitute particle if not accepted and change to ame-name
        ame_name = AME::get_instance()->GetSymbol(this->heavy_particle_substitute).value();
    }
    return this->f1_ReactionLost_CorEff[ame_name]->Eval(Ekin_Lab);
}

void E15190ReactionLost::_define_acceptable_particles()
{
    for (auto &name : this->default_acceptable_particles)
    {
        this->acceptable_particles.push_back(name);

        auto optional_symbol = AME::get_instance()->GetSymbol(name);
        if (!optional_symbol)
        {
            std::cout << "E15190ReactionLost::E15190ReactionLost: " << name << " is not found in AME database." << std::endl;
            continue;
        }

        this->acceptable_particles.push_back(optional_symbol.value());
    }
    return;
}

void E15190ReactionLost::set_heavy_particle_substitute(const std::string &name)
{
    this->heavy_particle_substitute = AME::get_instance()->GetSymbol(name).value_or(default_heavy_particle_substitute);
}

bool E15190ReactionLost::_accepted(const std::string &name)
{
    auto iter = std::find(this->acceptable_particles.begin(), this->acceptable_particles.end(), name);
    return (iter != this->acceptable_particles.end());
}