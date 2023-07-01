#include "AME.hh"

void AME::ReadAMETable(const std::string &filename)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (this->IsLoaded)
    {
        return;
    }
    fs::path path = filename;
    if (!fs::exists(path))
    {
        fs::path PROJECT_DIR = std::getenv("PROJECT_DIR");
        fs::path default_path = PROJECT_DIR / DefaultPath;
        path = default_path;
    }

    std::ifstream infile(path.c_str());
    infile.ignore(99, '\n');

    std::string symbol;
    int Z, A;
    double mass;                       // MeV/c^2
    double binding_enregy_per_nucleon; // MeV

    while (infile >> symbol >> Z >> A >> mass >> binding_enregy_per_nucleon)
    {
        std::pair pair = std::make_pair(Z, A);
        this->MassTable[symbol] = mass;
        this->ZATable[symbol] = pair;
        this->SymbolTable[pair] = symbol;
    }
    this->IsLoaded = true;
    infile.close();
    return;
}

std::optional<double> AME::GetMass(const std::string &symbol) const
{
    auto lower = [](const std::string &input) -> std::string
    {
        std::string output = input;
        std::transform(output.begin(), output.end(), output.begin(), [](unsigned char c)
                       { return std::tolower(c); });
        return output;
    };

    std::string lower_symbol = lower(symbol);

    if (this->MassTable.count(lower_symbol) == 0)
    {
        if (this->alias.count(symbol) == 0)
        {
            return std::nullopt;
        }
        std::string ame_symbol = this->alias.at(symbol);
        return this->MassTable.at(ame_symbol);
    }
    return this->MassTable.at(lower_symbol);
}

std::optional<std::string> AME::GetSymbol(const int &Z, const int &A) const
{
    if (this->SymbolTable.count({Z, A}) == 0)
    {
        return std::nullopt;
    }
    return this->SymbolTable.at({Z, A});
}

std::optional<std::string> AME::GetSymbol(const std::string &alias) const
{
    // if alias is a symbol, return itself
    if (this->ZATable.count(alias) == 1)
    {
        return alias;
    }
    // if alias is found in alias table, return the corresponding symbol
    return (this->alias.count(alias) == 0) ? std::nullopt : std::optional<std::string>(this->alias.at(alias));
}

std::optional<double> AME::GetMass(const int &Z, const int &A) const
{
    if (this->SymbolTable.count({Z, A}) == 0)
    {
        return this->_GetMassUnphysical(Z, A);
    }
    std::string symbol = this->SymbolTable.at({Z, A});
    return this->MassTable.at(symbol);
}

std::optional<double> AME::_GetMassUnphysical(const int &Z, const int &A) const
{
    if (Z > 0 && A >= Z)
    {
        return A * this->NucleonMass;
    }
    return std::nullopt;
}

std::optional<std::pair<int, int>> AME::GetZA(const std::string &symbol) const
{
    auto lower = [](const std::string &input) -> std::string
    {
        std::string output = input;
        std::transform(output.begin(), output.end(), output.begin(), [](unsigned char c)
                       { return std::tolower(c); });
        return output;
    };

    std::string lower_symbol = lower(symbol);

    if (this->ZATable.count(lower_symbol) == 0)
    {
        if (this->alias.count(symbol) == 0)
        {
            return std::nullopt;
        }
        return this->ZATable.at(alias.at(symbol));
    }
    return this->ZATable.at(lower_symbol);
}

std::optional<int> AME::GetZ(const std::string &symbol) const
{
    auto za = this->GetZA(symbol);
    return (za) ? std::optional<int>(za->first) : std::nullopt;
}

std::optional<int> AME::GetA(const std::string &symbol) const
{
    auto za = this->GetZA(symbol);
    return (za) ? std::optional<int>(za->second) : std::nullopt;
}

bool AME::IsPhysical(const int &Z, const int &A) const
{
    return this->GetMass(Z, A) != std::nullopt;
}
bool AME::IsPhysical(const std::string &symbol) const
{
    return this->GetMass(symbol) != std::nullopt;
}