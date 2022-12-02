#ifndef particle_h
#define particle_h

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "TMath.h"

struct hira_particle
{
    // general properties
    int nid, zid;
    double px, py, pz;

    double ekinlab, ekincms, plab, pcms;
    double etrans, pt, pzcms, rapidity_lab, rapidity_cms;
    double pseudo_rapidity_lab, pseudo_rapidity_cms;
    double thetalab, thetacms, phi;
    int aid;

    int fnumtel, fnumcsi, fnumstripf, fnumstripb;
    double geoeff = 1.;
    double reaction_eff = 1.;

    std::string name = "";
    std::string coordinate = "hira";
    void set_uball_coordinate()
    {
        this->coordinate = "uball";
    }
    void set_geometrical_efficiency(const double &eff) { this->geoeff = eff; }
    void set_reaction_efficiency(const double &eff) { this->reaction_eff = eff; }
    void set_detector_id(const int &tel, const int &csi, const int &ef, const int &eb)
    {
        this->fnumtel = tel;
        this->fnumcsi = csi;
        this->fnumstripf = ef;
        this->fnumstripb = eb;
    }

    void init(const double &betacms);
};

#endif
