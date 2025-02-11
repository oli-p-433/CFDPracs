#include <vector>
#include <string>
#include <functional>
#include <map>
#include "EOS.H"
#include "riem.H"
#include "operators.H"

class solver{
    public:
        std::vector< std::array<double,3> > u1Plus1;
        std::vector< std::array<double,3> > u2Plus1;

        // level set
        std::vector< double > phi;
        std::vector< double > phiPlus1;

        // GFM
        std::vector<int> interface_cells;
        void findBoundary();
        void interfaceRiem();
        void setGhostFluids();
        void reinitPhi();

        void init(std::vector< std::array<double,3> >,std::vector< std::array<double,3> >);
        void phiInit(std::vector<double>);


        void setDt();
        void laxFriedrichs();
        std::array<double,3> fBurgersLF(std::array<double,3>)const;
        std::array<double,3> finvectLF(std::array<double,3>)const;
        std::array<double,3> fEuler(std::array<double,3>, EOS* eos);
        void Richt();
        void FORCE();
        void SLIC();
        void godunov();
        std::array<double,3> laxFhalf(const std::array<double,3>, const std::array<double,3>, EOS* eos)const;
        std::array<double,3> LF(const std::array<double,3>, const std::array<double,3>, EOS* eos)const;
        void pointsUpdate();
        void phiUpdate();
        void run();

        // slope limiting
        void calcHalfSlopes();
        void calcr();
        std::array<double,3> static calcSlope(double omega, std::array<double,3> slopeleft, std::array<double,3> sloperight);
        void calcUBars();
        void updUBars();
        std::function<std::array<double,3>(std::array<double,3>)> slopeLim;
        
        std::array<double,3> minbee(std::array<double,3>);
        std::array<double,3> superbee(std::array<double,3>);
        std::array<double,3> vanLeer(std::array<double,3>);

        //

        void transmissiveBC();

        std::function<std::array<double,3>(std::array<double,3>,EOS*)> flux;

        void setWriteInterval(double);
        bool checkWrite;
        void writeData(std::string) const;

        int ghosts();
        double get_dx()const;


        double dt;

        std::string dirName;
        std::string fileName;

        std::map<std::string,int> varMap;

        // constructor

        solver(double,double,double,double,int,int,double);

        // utility functions

        static void print_arr(std::vector< std::array<double,3> >, int);
        static void print_vect(std::array<double,3>);

        // Source terms

        void sourceUpdate();
        std::vector< std::array<double,3> > sourceTerm(const std::vector< std::array<double,3> > uVals, int);
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

        std::vector< std::array<double,3> > u1; //indexed as u[i][j] (ith array, jth row)
        std::vector< std::array<double,3> > fluxes1;

        // second material
        std::vector< std::array<double,3> > u2; //indexed as u[i][j] (ith array, jth row)
        std::vector< std::array<double,3> > fluxes2;

        std::vector< std::array<double,3> > halfSlopes1;
        std::vector< std::array<double,3> > r1;
        std::vector< std::array<double,3> > u1BarL;
        std::vector< std::array<double,3> > u1BarR;
        std::vector< std::array<double,3> > u1BarLupd;
        std::vector< std::array<double,3> > u1BarRupd;
        // material 2
        std::vector< std::array<double,3> > halfSlopes2;
        std::vector< std::array<double,3> > r2;
        std::vector< std::array<double,3> > u2BarL;
        std::vector< std::array<double,3> > u2BarR;
        std::vector< std::array<double,3> > u2BarLupd;
        std::vector< std::array<double,3> > u2BarRupd;

        std::vector< std::array<double,3> > sourceResult;

        std::vector< std::array<double,3> > interfaceStates;  


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

