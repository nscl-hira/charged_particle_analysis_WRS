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
    long nevents;
    bool debug = false;
    std::string reaction = "";
    std::string output_file = "./test.root";
    std::array<int, 2> cut_on_multiplicity = {5, 128};
    std::array<double, 2> cut_on_impact_parameter = {0., 0.4};

    fs::path dir_data = "/data/HiRA_Cali";
    fs::path dir_geoeff = PROJECT_DIR / "database/e15190/hira/GeoEff";
    fs::path path_pixel_angles = PROJECT_DIR / "database/e15190/hira/Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat";
    fs::path path_bimp = PROJECT_DIR / "database/e15190/microball/bimp_mapping";
    fs::path path_runinfo = PROJECT_DIR / "database/e15190/RunInfo.dat";
    fs::path path_coverage = PROJECT_DIR / "database/e15190/hira/kinematic_cut.json";
    fs::path path_reaction_lost = PROJECT_DIR / "database/e15190/hira/reaction_loss.json";

    ArgumentParser(int argc, char *argv[])
    {
        options = {
            {"help", no_argument, 0, 'h'},
            {"nevents", required_argument, 0, 'n'},
            {"reaction", required_argument, 0, 'r'},
            {"output", required_argument, 0, 'o'},
            {"cut_on_multiplicity", required_argument, 0, 'c'},
            {"cut_on_impact_parameter", required_argument, 0, 'b'},
            {"debug", no_argument, 0, 'd'},
            {"dir_data", required_argument, 0, 0},
            {"dir_geoeff", required_argument, 0, 0},
            {"path_pixel_angles", required_argument, 0, 0},
            {"path_bimp", required_argument, 0, 0},
            {"path_runinfo", required_argument, 0, 0},
            {"path_coverage", required_argument, 0, 0},
            {"path_reaction_lost", required_argument, 0, 0},
            {0, 0, 0, 0},
        };

        int option_index = 0;
        int opt;
        while ((opt = getopt_long(argc, argv, "hn:r:o:c:b:", options.data(), &option_index)) != -1)
        {
            switch (opt)
            {
            case 'n':
            {
                this->nevents = std::stol(optarg);
                break;
            }
            case 'r':
            {
                this->reaction = optarg;
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

            case 'd':
            {
                this->debug = true;
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
            case 0:
            {
                std::string field = options[option_index].name;
                if (field == "dir_data")
                    this->dir_data = optarg;
                else if (field == "dir_geoeff")
                    this->dir_geoeff = optarg;
                else if (field == "path_pixel_angles")
                    this->path_pixel_angles = optarg;
                else if (field == "path_bimp")
                    this->path_bimp = optarg;
                else if (field == "path_runinfo")
                    this->path_runinfo = optarg;
                else if (field == "path_coverage")
                    this->path_coverage = optarg;
                else if (field == "path_reaction_lost")
                    this->path_reaction_lost = optarg;
                else
                    std::cout << "Got unknown option: " << field << std::endl;
                break;
            }
            default:
            {
                std::cout << "Got unknown parse returns: " << optarg << std::endl;
            }
            }
        }

        // handle invalid cuts on events
        if (this->cut_on_multiplicity[0] < UBALL_MULTI_MIN)
        {
            std::cout << "Warning: cut_on_multiplicity[0] < " << UBALL_MULTI_MIN << " , set to UBALL_MULTI_MIN." << std::endl;
            this->cut_on_multiplicity[0] = UBALL_MULTI_MIN;
        }
        if (this->cut_on_impact_parameter[0] < 0.)
        {
            std::cout << "Warning: cut_on_impact_parameter[0] < 0, set to 0." << std::endl;
            this->cut_on_impact_parameter[0] = 0.;
        }

        // handle paths
        this->path_bimp = this->path_bimp / (this->reaction + ".dat");
    }

    void help()
    {
        const char *msg = R"(
            -n      number of events to analyze.
            -r      reaction tag, e.g. `Ca48Ni64E140`.
            -o      ROOT file output path.
            -c      cut on multiplicity, e.g. `5 128`.
            -b      cut on impact parameter, e.g. `0. 0.4`.
            
            --dir_data              directory of the data files.         
            --dir_geoeff            directory of the GeometricEfficiency files.
            --path_pixel_angles     path to the file mapping detector id to angles.
            --path_bimp             directory containing the impact parameter mapping files for different reaction systems.
            --path_runinfo          path to the file containing information of each experiment Run.
            --path_coverage         path to the json file containing the kinetmatic cut of HiRA10.
            --path_reaction_lost    path to the json file containing the reaction lost coefficents for different particles in HiRA10.
        )";
        std::cout << msg << std::endl;
    }

    void report()
    {
        std::cout << "--------------------------------------------------------------------------" << std::endl;
        std::cout << "reaction: " << this->reaction << std::endl;
        std::cout << "cut_on_multiplicity: " << this->cut_on_multiplicity[0] << " " << this->cut_on_multiplicity[1] << std::endl;
        std::cout << "cut_on_impact_parameter: " << this->cut_on_impact_parameter[0] << " " << this->cut_on_impact_parameter[1] << std::endl;

        // print the paths to the input files
        std::cout << '\n';
        std::cout << "dir_data: " << this->dir_data << std::endl;
        std::cout << "dir_geoeff: " << this->dir_geoeff << std::endl;
        std::cout << "path_pixel_angles: " << this->path_pixel_angles << std::endl;
        std::cout << "path_bimp: " << this->path_bimp << std::endl;
        std::cout << "path_runinfo: " << this->path_runinfo << std::endl;
        std::cout << "output_file: " << this->output_file << std::endl;
        std::cout << "--------------------------------------------------------------------------" << std::endl;
    }

protected:
    std::vector<option> options;
};