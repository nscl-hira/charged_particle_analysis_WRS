#ifndef ame_hh
#define ame_hh
#include <map>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

std::map<std::pair<int, int>, double> get_ame_mass_table()
{
    std::map<std::pair<int, int>, double> ame_mass_table;
    fs::path PROJECT_DIR = std::getenv("PROJECT_DIR");
    fs::path path = PROJECT_DIR / "database/ame/ame_mass.txt";
    std::ifstream infile(path.c_str());
    infile.ignore(99, '\n');

    std::string symbol;
    int Z, A;
    double mass;                       // MeV/c^2
    double binding_enregy_per_nucleon; // MeV
    while (infile >> symbol >> Z >> A >> mass >> binding_enregy_per_nucleon)
    {
        std::pair pair = std::make_pair(Z, A);
        ame_mass_table[pair] = mass;
    }

    return ame_mass_table;
}

#endif