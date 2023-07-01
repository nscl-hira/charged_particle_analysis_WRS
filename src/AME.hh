#ifndef AME_hh
#define AME_hh

#include <map>
#include <string>
#include <optional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

#include <mutex>

class AME
{
public:
    AME(AME const &) = delete;
    AME &operator=(AME const &) = delete;
    ~AME() {}

    void ReadAMETable(const std::string &filename = "");
    void PrintTable(const int &head = 10) const;

    // We don't want getters to be thread-safe in our case.
    std::optional<double> GetMass(const std::string &symbol) const;
    std::optional<double> GetMass(const int &Z, const int &A) const;
    std::optional<double> _GetMassUnphysical(const int &Z, const int &A) const;

    std::optional<double> GetBindingEnergyPerNucleon(const std::string &symbol) const;

    std::optional<std::pair<int, int>> GetZA(const std::string &symbol) const;
    std::optional<int> GetA(const std::string &symbol) const;
    std::optional<int> GetZ(const std::string &symbol) const;

    std::optional<std::string> GetSymbol(const int &Z, const int &A) const;
    std::optional<std::string> GetSymbol(const std::string &alias) const;

    bool IsPhysical(const int &Z, const int &A) const;
    bool IsPhysical(const std::string &symbol) const;

    static AME *get_instance(const std::string &filename = "")
    {
        // Static local variable initialization is thread-safe
        // and will be initialized only once.
        static AME instance(filename);
        return &instance;
    }

private:
    explicit AME(const std::string &filename) : IsLoaded(false)
    {
        this->ReadAMETable(filename);
    }

    std::mutex m_mutex;
    bool IsLoaded;
    std::map<std::string, double> MassTable;
    std::map<std::string, double> BindingEnergyTable; // per nucleon
    std::map<std::string, std::pair<int, int>> ZATable;
    std::map<std::pair<int, int>, std::string> SymbolTable;

protected:
    // if Z and A are not physical, return A * NucleonMass as mass
    const double NucleonMass = 938.272;
    // default path relative to project directory
    const fs::path DefaultPath = "database/ame/ame_mass.txt";
    const std::map<std::string, std::string> alias = {
        {"n", "n1"},
        {"p", "h1"},
        {"d", "h2"},
        {"t", "h3"},
        {"3He", "he3"},
        {"4He", "he4"},
        {"alpha", "he4"},
    };
};

#endif