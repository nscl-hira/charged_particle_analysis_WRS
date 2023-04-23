#include "AME.hh"

void AME::ReadAMETable(const std::string &filename)
{
    std::string path = filename;
    fs::path project_dir = std::getenv("PROJECT_DIR");
    fs::path dat_path = project_dir / "database/ame/ame_mass.txt";

    if (path.empty())
    {
        path = dat_path;
    }
    std::ifstream infile(path.c_str());
    infile.ignore(99, '\n');
    int Z, A;
    double mass;
    while (infile >> Z)
    {
        infile >> A >> mass;
        this->ame_mass_table[{Z, A}] = mass;
    }
}

double AME::GetMass(const int &Z, const int &A)
{
    if (this->ame_mass_table.size() == 0)
    {
        this->ReadAMETable();
    }

    if (this->ame_mass_table.count({Z, A}) == 0)
    {
        if (Z > 0 && A >= Z)
        {
            // std::cout << "Warning: AME mass table does not contain Z = " << Z << ", A = " << A << ". Returning Mass = A * 938.272 MeV/c^2" << std::endl;
            return A * 938.272;
        }
        return 1e10;
    }
    return this->ame_mass_table[{Z, A}];
}