#include <vector>
#include <string>
#include <functional>
#include <map>
#include "EOS.H"
#include "riem.H"
#include "operators.H"
#include "boundary.H"

class solver{

    public:

        enum Cons{FRAC = 0, RHO1, RHO2, XMOM, YMOM, ENE};
        enum Prim{UX = 3, UY, PRES};
        enum dir{YDIR=0, XDIR};

        bool PRIM;


        void init(std::vector< std::vector< std::array<double,6> > >);


        void setDt();
        void laxFriedrichs();
        std::array<double,6> fBurgersLF(std::array<double,6>)const;
        std::array<double,6> finvectLF(std::array<double,6>)const;
        std::array<double,6> fEuler(std::array<double,6>, EOS*,bool);
        void Richt();
        void FORCE();
        void SLIC(bool);
        void MUSCL(bool);
        void HLLCGodunov(bool);

        bool isGodunov;

        std::function<void(bool)> fluxMethod;

        std::array<double,6> HLLC(std::array<double,6>,std::array<double,6>, EOS*, int, int, bool);
        std::array<double,6> laxFhalf(const std::array<double,6>, const std::array<double,6>, EOS*, bool)const;
        std::array<double,6> LF(const std::array<double,6>, const std::array<double,6>, EOS*, bool)const;
        void pointsUpdate(bool);
        void RK2(bool);
        void run();

        // slope limiting
        void calcHalfSlopes(bool,bool);
        void calcr(bool);
        std::array<double,6> static calcSlope(double omega, std::array<double,6> slopeleft, std::array<double,6> sloperight);
        void calcUBars(bool,bool);
        void updUBars(bool,bool);
        std::function<std::array<double,6>(std::array<double,6>)> slopeLim;
        
        std::array<double,6> minbee(std::array<double,6>);
        std::array<double,6> superbee(std::array<double,6>);
        std::array<double,6> vanLeer(std::array<double,6>);

        //

        void transmissiveBC();
        std::function<void(fluid&)> setBCs;

        std::function<std::array<double,6>(std::array<double,6>,EOS*,bool)> flux;

        void setWriteInterval(double);
        bool checkWrite;
        void writeData() const;

        int ghosts();
        std::array<double,2> get_dxdy()const;

        boundary& get_boundary(){
            return bound;
        };


        double dt;

        std::string dirName;
        std::string fileName;

        std::map<std::string,int> varMap;

        // constructor

        solver(double,double,double,double,double,double,int,int,int,double);

        // utility functions

        static void print_arr(std::vector< std::array<double,6> >, int);
        static void print_vect(std::array<double,6>);

        // Source terms

        void sourceUpdate();
        void allaireSource(bool,double);
        std::vector< std::array<double,6> > sourceTerm(const std::vector< std::array<double,6> > uVals, int);
        int alpha;

        // EOS handling

        void setEOS(std::array<EOS*,2> eosInstances){
            eos = eosInstances;
        };

        std::array<EOS*,2> eos = {nullptr,nullptr};


 

    protected:

        fluid f;
        boundary bound;

        double timeMulti;
        int nCellsX, nCellsY, nG;
    
        double dx;
        double dy;

        double cSound;

        int splitFlip{0};

        std::vector< std::array<double,6> > sourceResult;

        double x0,x1,y0,y1;

        double startTime;
        double endTime;
        double time;
        double slopeWeight;

        double cour;

        double writeInterval;

        double dtReducer = 1;
};

