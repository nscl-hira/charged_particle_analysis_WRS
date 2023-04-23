#include "LaserAngles.hh"

LaserAngles::LaserAngles(const std::string &path)
{
    std::ifstream infile(path.c_str());
    if (!infile)
    {
        std::cout << "file not found." << std::endl;
        std::exit(1);
    }
    infile.ignore(99, '\n');

    int HiraIndex, StripXIndex, StripYIndex;
    double theta, phi;

    auto strip = [](std::string &line) -> void
    {
        while (line[0] == ' ')
        {
            line = line.substr(1, line.length() - 1);
        }
        while (line[line.length() - 1] == ' ')
        {
            line = line.substr(0, line.length() - 1);
        }
    };

    for (std::string line; std::getline(infile, line);)
    {
        strip(line);
        std::istringstream ss(line);
        ss >> HiraIndex >> StripXIndex >> StripYIndex >> theta >> phi;
        this->ThetaLab[HiraIndex][StripXIndex][StripYIndex] = theta;
        this->Phi[HiraIndex][StripXIndex][StripYIndex] = phi;
    }
}