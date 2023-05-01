#include "TFile.h"
#include "TChain.h"

#include "AME.hh"
#include "Particle.hh"
#include "Physics.hh"
#include "e15190/RunInfo.hh"
#include "e15190/HiRA.hh"
#include "e15190/Reader.hh"

fs::path PROJECT_DIR = std::getenv("PROJECT_DIR");
constexpr int UBALL_MULTI_MIN = 5;
constexpr double TDC_TRIGGER_THRESHOLD = -9990;

class ArgumentParser
{
public:
    std::string reaction = "";

    fs::path dir_data = "/data/HiRA_Cali";
    fs::path dir_geoeff = PROJECT_DIR / "database/e15190/GeoEff";
    fs::path path_pixel_angles = PROJECT_DIR / "database/e15190/Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat";
    fs::path path_bimp = PROJECT_DIR / "database/e15190/microball/bimp_mapping/Ca48Ni64E140.dat";
    fs::path path_runinfo = PROJECT_DIR / "database/e15190/RunInfo.dat";

    std::string output_file = "./test.root";
    std::array<int, 2> cut_on_multiplicity = {5, 128};
    std::array<double, 2> cut_on_impact_parameter = {0., 0.4};

    ArgumentParser(int argc, char *argv[])
    {
        options = {
            {"help", no_argument, 0, 'h'},
            {"reaction", required_argument, 0, 'r'},
            {"dir_data", required_argument, 0, 'i'},
            {"output", required_argument, 0, 'o'},
            {"cut_on_multiplicity", required_argument, 0, 'c'},
            {"cut_on_impact_parameter", required_argument, 0, 'b'},
            {0, 0, 0, 0},
        };

        int option_index = 0;
        int opt;
        while ((opt = getopt_long(argc, argv, "hir:o:c:b:", options.data(), &option_index)) != -1)
        {
            switch (opt)
            {
            case 'r':
            {
                this->reaction = optarg;
                break;
            }

            case 'i':
            {
                this->dir_data = optarg;
                break;
            }

            case 'o':
            {
                this->output_file = optarg;
                break;
            }
            case 'c':
            {
                std::string cut_on_multiplicity = optarg;
                std::cout << cut_on_multiplicity << std::endl;

                std::istringstream iss(cut_on_multiplicity);
                iss >> this->cut_on_multiplicity[0] >> this->cut_on_multiplicity[1];
                break;
            }

            case 'b':
            {
                std::string cut_on_impact_parameter = optarg;
                std::istringstream iss(cut_on_impact_parameter);
                iss >> this->cut_on_impact_parameter[0] >> this->cut_on_impact_parameter[1];
                break;
            }

            case '?':
            {
                std::cout << "Got unknown option." << std::endl;
                break;
            }

            case 'h':
            {
                this->help();
                exit(1);
            }

            default:
            {
                std::cout << "Got unknown parse returns: " << optarg << std::endl;
            }
            }
        }

        if (this->cut_on_multiplicity[0] < UBALL_MULTI_MIN)
        {
            std::cout << "Warning: cut_on_multiplicity[0] < " << UBALL_MULTI_MIN << " , set to UBALL_MULTI_MIN." << std::endl;
            this->cut_on_multiplicity[0] = UBALL_MULTI_MIN;
        }

        std::cout << "reaction: " << this->reaction << std::endl;
        std::cout << "output_file: " << this->output_file << std::endl;
        std::cout << "cut_on_multiplicity: " << this->cut_on_multiplicity[0] << " " << this->cut_on_multiplicity[1] << std::endl;
        std::cout << "cut_on_impact_parameter: " << this->cut_on_impact_parameter[0] << " " << this->cut_on_impact_parameter[1] << std::endl;
    }
    void help()
    {
        const char *msg = R"(
            -r      reaction tag, e.g. `Ca48Ni64E140`
            -i      a list of input ROOT files, separated by space.
            -o      ROOT file output path.
        )";
        std::cout << msg << std::endl;
    }

protected:
    std::vector<option> options;
};