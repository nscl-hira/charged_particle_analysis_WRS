#ifndef AME_hh
#define AME_hh

#include <array>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
namespace fs = std::filesystem;

struct AME
{
    std::map<std::array<int, 2>, double> ame_mass_table;
    void ReadAMETable(const std::string &filename = "");
    double GetMass(const int &Z, const int &A);
};

#endif