#ifndef particle_h
#define particle_h

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "TMath.h"

struct particle
{
    // general properties
    int nid, zid;
    double px, py, pz;

    double ekinlab, ekincms, plab, pcms;
    double etrans, pt, pzcms, rapidity_lab, rapidity_cms;
    double pseudo_rapidity_lab, pseudo_rapidity_cms;
    double thetalab, thetacms, phi;
    int aid;
    std::string name = "";
    std::string coordinate = "hira";
    void set_uball_coordinate()
    {
        this->coordinate = "uball";
    }

    void init(const double &betacms);
};

#endif
