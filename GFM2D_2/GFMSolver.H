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

        enum Cons{RHO=0, XMOM, YMOM, ENE};
        enum Prim{UX = 1, UY, PRES};
        enum dir{YDIR=0, XDIR};

        bool PRIM;

        double timeMulti;



        void init(std::vector< std::vector< std::array<double,4> > >, fluid&);
        std::array<double,4> set_vals(int, double, double, double, double);
        fluid& get_fluid(bool);
        void phiInit(std::vector< std::vector < double > >);


        //GFM
        void findBoundary();
        void calcInterfaceNormals(fluid&);
        void calcInterface(fluid&);
        void interpInterfaceStates(fluid&);

        std::array<double,4> bilinear4(double,double,fluid&,std::array<double,2>,int);
        std::array<double,4> bilinear5(double,double,fluid&,std::array<double,2>,int);
        void interfaceRiem(fluid&);
        void resolveVelocities(fluid&);
        void calcStarStates(fluid&);
        void calcPhiGrad();
        std::vector < std::vector < std::array<double,2 > > > phiNormals;
        void calcnDotPhiNormals();
        void setInterface();
        void setGhostFluids();
        std::array<double,4> solveQ(bool,bool,std::array<double,4>,std::array<double,4>,double);
        void solveSweepPoint(int,int);
        void GFFastSweeping();
        void GFFastSweeping2();
        void fastSweepExtrapolation(std::vector<std::vector<std::array<double, 4>>>& grid,
            std::vector<std::array<int, 2>>& interfaceCells,std::vector<std::vector<double>>& phi,
            fluid& f, int numSweeps);
        void eikonalDt();
        void dtExtrap();
        double extrapDt;

        std::array<double,4> exactRiemann(std::array<double,4>,std::array<double,4>, EOS*, bool);
        

        void setDt();
        void laxFriedrichs();
        std::array<double,4> fBurgersLF(std::array<double,4>)const;
        std::array<double,4> finvectLF(std::array<double,4>)const;
        std::array<double,4> fEuler(std::array<double,4>, EOS*,bool);
        void Richt();
        void FORCE();
        void SLIC(bool,fluid&,EOS*);
        void MUSCL(bool,fluid&,EOS*);
        void HLLCGodunov(bool,fluid&,EOS*);

        bool isGodunov;

        std::function<void(bool,fluid&,EOS*)> fluxMethod;
        std::function<std::array<double,4>(std::array<double,4>,std::array<double,4>, EOS*)> fluxRiemannSolver;

        std::array<double,4> HLLC(std::array<double,4>,std::array<double,4>, EOS*, int, int, bool);
        std::array<double,4> laxFhalf(const std::array<double,4>, const std::array<double,4>, EOS*, bool)const;
        std::array<double,4> LF(const std::array<double,4>, const std::array<double,4>, EOS*, bool)const;


        void pointsUpdate(bool,fluid&);
        void RK2(bool, fluid&,EOS*);
        void runGFM();
        void run();

        // Level set

        void initialisePhi();
        void reinitPhiFastSweeping(int);
        void reinitPhiIterative();
        void solveEikonalPoint(int,int);

        void updateLevelSetPoint(std::vector<std::vector<double>> &phi, int i, int j);
        void reinitialiseLevelSet(std::vector<std::vector<double>> &phi, int maxIter);

        void neighbourAvg(fluid& f, int, int);
        std::array<double,4> neighbourAvgInterp(fluid& f, int, int);
        void fixFreshlyCleared();

        void phiBC();
        void phiUpdate(int);
        void phiMUSCLupdate(int);

        // slope limiting
        void calcHalfSlopes(bool,bool,fluid&);
        void calcr(bool,fluid&);
        std::array<double,4> static calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight);
        void calcUBars(bool,bool,fluid&,EOS*);
        void updUBars(bool,bool,fluid&,EOS*);
        std::function<std::array<double,4>(std::array<double,4>)> slopeLim;
        
        std::array<double,4> minbee(std::array<double,4>);
        std::array<double,4> superbee(std::array<double,4>);
        std::array<double,4> vanLeer(std::array<double,4>);

        //

        void transmissiveBC(fluid&);
        std::function<void(fluid&)> setBCs;

        std::function<std::array<double,4>(std::array<double,4>,EOS*,bool)> flux;

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

        static void print_arr(std::vector< std::array<double,4> >, int);
        static void print_vect(std::array<double,4>);
        void printInterfaceArrays(std::vector<std::vector<double>>, fluid&, const std::string&);

        // Source terms

        void sourceUpdate();
        void allaireSource(bool,double);
        std::vector< std::array<double,4> > sourceTerm(const std::vector< std::array<double,4> > uVals, int);
        int alpha;

        // EOS handling

        void setEOS(std::array<EOS*,2> eosInstances){
            eos = eosInstances;
        };

        std::array<EOS*,2> eos = {nullptr,nullptr};


 

    protected:

        fluid fluid1, fluid2;
        boundary bound;

        int nCellsX, nCellsY, nG;
    
        double dx;
        double dy;

        double cSound;

        int splitFlip{0};

        std::vector< std::array<double,4> > sourceResult;

        double x0,x1,y0,y1;

        double startTime;
        double endTime;
        double time;
        double slopeWeight;

        double cour;

        double writeInterval;

        double dtReducer = 1;

        // level set
        std::vector< std::vector < double > > phi;
        std::vector< std::vector < double > > phiPlus1;
        std::vector< std::vector < double > > phiOld;
        std::vector<std::array<int, 2>> freshlyCleared;

        std::vector< std::vector < double > > phiGrads;

        double pInf1, pInf2;
};

