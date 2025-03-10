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

        enum Cons{FRAC = 0, RHO1, RHO2, XMOM, ENE};
        enum Prim{UX = 3, PRES};

        bool PRIM;


        void init(std::vector< std::vector< std::array<double,5> > >);


        void setDt();
        void laxFriedrichs();
        std::array<double,5> fBurgersLF(std::array<double,5>)const;
        std::array<double,5> finvectLF(std::array<double,5>)const;
        std::array<double,5> fEuler(std::array<double,5>, EOS* eos);
        void Richt();
        void FORCE();
        void SLIC();
        void MUSCL();
        void HLLCGodunov();

        std::function<void()> fluxMethod;

        std::array<double,5> HLLC(std::array<double,5>,std::array<double,5>, EOS*, int, int);
        std::array<double,5> laxFhalf(const std::array<double,5>, const std::array<double,5>, EOS* eos)const;
        std::array<double,5> LF(const std::array<double,5>, const std::array<double,5>, EOS* eos)const;
        void pointsUpdate();
        void run();

        // slope limiting
        void calcHalfSlopes(bool);
        void calcr();
        std::array<double,5> static calcSlope(double omega, std::array<double,5> slopeleft, std::array<double,5> sloperight);
        void calcUBars(bool);
        void updUBars(bool);
        std::function<std::array<double,5>(std::array<double,5>)> slopeLim;
        
        std::array<double,5> minbee(std::array<double,5>);
        std::array<double,5> superbee(std::array<double,5>);
        std::array<double,5> vanLeer(std::array<double,5>);

        //

        void transmissiveBC();
        std::function<void(fluid&)> setBCs;

        std::function<std::array<double,5>(std::array<double,5>,EOS*)> flux;

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

        static void print_arr(std::vector< std::array<double,5> >, int);
        static void print_vect(std::array<double,5>);

        // Source terms

        void sourceUpdate();
        void allaireSource();
        std::vector< std::array<double,5> > sourceTerm(const std::vector< std::array<double,5> > uVals, int);
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

        

        std::vector< std::array<double,5> > sourceResult;

        double x0,x1,y0,y1;

        double startTime;
        double endTime;
        double time;
        double slopeWeight;

        double cour;

        double writeInterval;
};

