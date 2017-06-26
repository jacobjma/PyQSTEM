#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include "data_containers.h"
#include "stemlib.h"

namespace shapes {

    class QSTEM {
    public:
        MULS muls;
        WavePtr wave;
        int trans_array_state;
        int wave_state;
        int atoms_state;
        int box_state;

        std::string mode;
        QSTEM(std::string mode);
        ~QSTEM();

        void set_atoms(int natom, int atomKinds, const std::vector< std::vector<double> > & pos,
                        const std::vector<double> & occ,const std::vector<double> & q,
                        const std::vector<int> & Znum);
        void set_positions(int natom, const std::vector< std::vector<double> > & pos);
        void set_box(const std::vector<double> &, int nonPeriod, int nonPeriodZ, float cellDiv);

        void get_resolution(float* resolutionX,float* resolutionY);
        void get_potential_samples(int* potNx, int* potNy, int* slices);
        void get_probe_samples(int* nx, int* ny);
        void get_scan_range(float* scanXStart, float* scanXStop, int* scanXN, float* scanYStart, float* scanYStop, int* scanYN);
        void get_potential_extent(float* potOffsetX,float* potOffsetY,float* potSizeX,float* potSizeY);
        void get_energy(float* energy);

        void allocate_potential();
        void build_potential(int slices);
        void set_potential(const std::vector< std::vector< std::vector< std::vector<double> > > > & potential,
                const std::vector<int> & size, const std::vector<double> & extent);
        std::vector <std::vector <std::vector <std::vector <double> > > >get_potential_or_transfunc(float * resolutionX,float* resolutionY,float* sliceThickness);

        void calculate_transfunc();
        void set_scan_range(float scanXStart,float scanXStop,int scanXN,float scanYStart,float scanYStop,int scanYN);
        void free_wave();
        void build_wave(int type,float v0,int nx,int ny,float resolutionX,float resolutionY);
        void build_probe(float v0,float alpha,int nx,int ny,float resolutionX,float resolutionY,std::unordered_map<std::string,float> abberations);
        void set_wave(std::vector <std::vector <std::vector <double> > >wave,float v0,int nx,int ny,float resolutionX,float resolutionY);
        std::vector <std::vector <std::vector <double> > > get_wave(float * resolutionX,float* resolutionY,float *v0);
        //void create_detector(std::string name,float rInside,float rOutside,float shiftX,float shiftY,int detectorNum);

        void create_detectors(std::vector<std::string> names,std::vector<std::vector<double>> properties, int detectorNum);
        void run(int displayProgInterval);
        std::vector <std::vector <double> > read_detector(int detNum);

    };
}
