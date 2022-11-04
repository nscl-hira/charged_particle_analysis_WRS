#include "badmap.hh"

badmap::badmap(const std::string &path, const std::string &version)
{
    this->BadMapVersion = version;
    this->BadMapPath = path;

    for (unsigned int iHira = 0; iHira < this->fNumTele; iHira++)
    {

        this->BadHiraMap[iHira] = 0;

        for (unsigned int StripIndex = 0; StripIndex < this->fNumStripf; StripIndex++)
        {
            this->BadStripMap[iHira][StripIndex][0] = 0;
            this->BadStripMap[iHira][StripIndex][1] = 0;
        }
        for (unsigned int iCsI = 0; iCsI < fNumCsI; iCsI++)
        {
            this->BadCsIMap[iHira][iCsI] = 0;
        }
    }

    Read_BadMap_Strip(Form("%s/BadMap_%s/Hira_BadMap_Strip_%s.dat", this->BadMapPath.c_str(), this->BadMapVersion.c_str(), this->BadMapVersion.c_str()));

    Read_BadMap_CsI(Form("%s/BadMap_%s/Hira_BadMap_CsI_%s.dat", this->BadMapPath.c_str(), this->BadMapVersion.c_str(), this->BadMapVersion.c_str()));

    Read_BadMap_Tele(Form("%s/BadMap_%s/Hira_BadMap_Tele_%s.dat", this->BadMapPath.c_str(), this->BadMapVersion.c_str(), this->BadMapVersion.c_str()));

    //     std::cout << "Initializing Badmap..." << std::endl;
}

badmap::~badmap() {}

void badmap::Read_BadMap_Strip(const std::string &path)
{

    std::ifstream infile(path.c_str());
    if (!infile)
    {
        std::cout << "file not found." << std::endl;
        return;
    }
    infile.ignore(99, '\n');
    while (infile)
    {
        int FrontBack = 0;
        int HiraIndex = 0;
        int StripIndex = 0;
        bool IsBad = 1;
        infile >> FrontBack >> HiraIndex >> StripIndex >> IsBad;
        BadStripMap[HiraIndex][StripIndex][FrontBack] = IsBad;

        infile.ignore(99, '\n');
        char first;
        infile.get(first);
        infile.putback(first);
    }
}

void badmap::Read_BadMap_CsI(const std::string &path)
{

    std::ifstream infile(path.c_str());
    if (!infile)
    {
        std::cout << "file not found." << std::endl;
        return;
    }

    infile.ignore(99, '\n');

    while (infile)
    {
        int HiraIndex = 0;
        int CsIIndex = 0;
        bool IsBad = 1;
        infile >> HiraIndex >> CsIIndex >> IsBad;
        BadCsIMap[HiraIndex][CsIIndex] = IsBad;

        infile.ignore(99, '\n');
        char first;
        infile.get(first);
        infile.putback(first);
    }
}

void badmap::Read_BadMap_Tele(const std::string &path)
{
    std::ifstream infile(path.c_str());
    if (!infile)
    {
        std::cout << "file not found." << std::endl;
        return;
    }
    infile.ignore(99, '\n');

    int HiraIndex = 0;
    bool IsBad = 0;
    while (infile >> HiraIndex)
    {
        infile >> IsBad;
        BadHiraMap[HiraIndex] = IsBad;
    }
}

bool badmap::IsBad_StripX(const int &HiraIndex, const int &StripX_Index)
{
    return (BadStripMap[HiraIndex][StripX_Index][0] == 1) || (StripX_Index == 0 || StripX_Index == 15 || StripX_Index == 16 || StripX_Index == 31);
    return 0;
}

bool badmap::IsBad_StripY(const int &HiraIndex, const int &StripY_Index)
{
    return (this->BadStripMap[HiraIndex][StripY_Index][1] == 1) || (StripY_Index == 0 || StripY_Index == 15 || StripY_Index == 16 || StripY_Index == 31);
}

bool badmap::IsBad_CsI(const int &HiraIndex, const int &CsI_Index)
{
    return BadCsIMap[HiraIndex][CsI_Index] == 1;
}

bool badmap::IsBad_CsI(const int &HiraIndex, const int &X_StripIndex, const int &Y_StripIndex)
{
    int CsIIndex = 0;
    if (X_StripIndex <= 15 && Y_StripIndex >= 16)
    {
        CsIIndex = 0;
    }
    else if (X_StripIndex <= 15 && Y_StripIndex <= 15)
    {
        CsIIndex = 1;
    }
    else if (X_StripIndex >= 16 && Y_StripIndex <= 15)
    {
        CsIIndex = 2;
    }
    else if (X_StripIndex >= 16 && Y_StripIndex >= 16)
    {
        CsIIndex = 3;
    }

    return IsBad_CsI(HiraIndex, CsIIndex);
}

bool badmap::IsBad_Hira(const int &HiraIndex)
{
    return this->BadHiraMap[HiraIndex] == 1;
}