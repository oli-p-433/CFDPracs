#include <vector>
#include <string>
#include <functional>
#include <map>
#include "EOS.H"
#include "riem.H"
#include "operators.H"

class solver{

    public:

        enum Cons{FRAC = 0, RHO1, RHO2, XMOM, ENE};
        enum Prim{UX = 3, PRES};

        bool PRIM;

        std::vector< std::array<double,5> > uPlus1;

        // level set
        std::vector< double > phi;
        std::vector< double > phiPlus1;



        void init(std::vector< std::array<double,5> >);


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

        std::array<double,5> HLLC(std::array<double,5>,std::array<double,5>, EOS*, int);
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

        std::function<std::array<double,5>(std::array<double,5>,EOS*)> flux;

        void setWriteInterval(double);
        bool checkWrite;
        void writeData() const;

        int ghosts();
        double get_dx()const;


        double dt;

        std::string dirName;
        std::string fileName;

        std::map<std::string,int> varMap;

        // constructor

        solver(double,double,double,double,int,int,double);

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

        double timeMulti;

    
        double dx;
        int nCells;
        int nGhost;
        double cSound;

        std::vector< std::array<double,5> > u; //indexed as u[i][j] (ith array, jth row)
        std::vector< std::array<double,5> > uPrev; //indexed as u[i][j] (ith array, jth row)
        std::vector< std::array<double,5> > fluxes;
        std::vector< double > sStars;


        std::vector< std::array<double,5> > uPrim;
        std::vector< std::array<double,5> > halfSlopes;
        std::vector< std::array<double,5> > r;
        std::vector< std::array<double,5> > uBarL;
        std::vector< std::array<double,5> > uBarR;
        std::vector< std::array<double,5> > uBarLupd;
        std::vector< std::array<double,5> > uBarRupd;

        std::vector< std::array<double,5> > sourceResult;

        std::vector< std::array<double,5> > interfaceStates;  


        double x0;
        double x1;
        double startTime;
        double endTime;
        double time;
        double slopeWeight;

        double cour;

        double writeInterval;
};


/*
std::array<double, 3> operator-(const std::array<double, 3>&, const std::array<double, 3>&);

std::array<double, 3> operator+(const std::array<double, 3>&, const std::array<double, 3>&);

std::array<double, 3> operator*(const double,const std::array<double, 3>&);

std::array<double, 3> operator*(const std::array<double, 3>&,const std::array<double, 3>&);

std::array<double, 3> operator/(const std::array<double, 3>&, const double);

std::array<double, 3> elementDivide(const std::array<double, 3>&, const std::array<double, 3>&);
*/

