#ifndef badmap_hh
#define badmap_hh

#include "TString.h"
#include <iostream>
#include <fstream>

class badmap
{
protected:
    static const unsigned int fNumTele = 12;
    static const unsigned int fNumCsI = 4;
    static const unsigned int fNumStripf = 32;
    static const unsigned int fNumStripb = 32;

public:
    badmap(const std::string &path, const std::string &version);
    ~badmap();
    std::string GetBadMapVersion() { return this->BadMapVersion; }

    bool BadHiraMap[fNumTele];
    bool BadCsIMap[fNumTele][fNumCsI];
    bool BadStripMap[fNumTele][fNumStripf][2]; //[HiraIndex][StripIndex][Front,Back]

    void Read_BadMap_Strip(const std::string &path);
    void Read_BadMap_CsI(const std::string &path);
    void Read_BadMap_Tele(const std::string &path);

    bool IsBad_StripX(const int &HiraIndex, const int &StripX_Index);
    bool IsBad_StripY(const int &HiraIndex, const int &StripY_Index);
    bool IsBad_CsI(const int &HiraIndex, const int &CsI_Index);
    bool IsBad_CsI(const int &HiraIndex, const int &StripX_Index, const int &StripY_Index);
    bool IsBad_Hira(const int &HiraIndex);

private:
    std::string BadMapVersion;
    std::string BadMapPath;
};

#endif