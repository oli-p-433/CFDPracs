#include <vector>
#include <string>
#include <functional>
#include <map>
#include "EOS.H"
#include "riem.H"
#include "boundary.H"
#include "operators.H"

class solver{
    public:

        enum Cons{RHO = 0, XMOM, YMOM, ENE};
        enum Cons2{BASE0 = 0, UX, UY, PRES};
        enum Cons3{YDIR=0,XDIR=1};

        bool direction; //

        std::vector<std::string> variables;

        void init(std::vector< std::vector< std::array<double,4> > >, fluid& fluid);\
        void phiInit(std::vector< std::vector < double > >);

        // GFM
        void findBoundary();
        void calcInterfaceNormals(fluid&);
        void calcInterface(fluid&);
        void interpInterfaceStates(fluid&);
        std::array<double,4> bilinear(double,double,fluid&,std::array<double,2>);
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
        void fastSweepExtrapolation(std::vector<std::vector<std::array<double, 4>>>& grid,
            std::vector<std::array<int, 2>>& interfaceCells,std::vector<std::vector<double>>& phi,
            fluid& f, int numSweeps);
        void eikonalDt();
        void dtExtrap();
        double extrapDt;


        void setDt();
        void laxFriedrichs(fluid&,EOS*);
        void richt(fluid&,EOS*);
        void FORCE(fluid&,EOS*);
        std::array<double,4> fBurgersLF(std::array<double,4>)const;
        std::array<double,4> finvectLF(std::array<double,4>)const;
        std::array<double,4> fEuler(std::array<double,4>,EOS*);

        void SLIC(fluid&, EOS*);
        void MUSCL(fluid&, EOS*);
        void godunov(fluid&, EOS*);

        std::function<void(fluid&, EOS*)> fluxMethod;
        std::function<std::array<double,4>(std::array<double,4>,std::array<double,4>, EOS*)> fluxRiemannSolver;
        
        std::array<double,4> HLL(std::array<double,4>,std::array<double,4>, EOS*);
        std::array<double,4> HLLC(std::array<double,4>,std::array<double,4>, EOS*);
        std::array<std::array<double,4>,2> HLLCStates(std::array<double,4>,std::array<double,4>, EOS*);
        

        std::array<double,4> laxFhalf(const std::array<double,4>, const std::array<double,4>, EOS*)const;
        std::array<double,4> LF(const std::array<double,4>, const std::array<double,4>, EOS*)const;
        void pointsUpdate(fluid&);
        void phiUpdate(int);

        void initialisePhi();
        void reinitPhiFastSweeping(int);
        void reinitPhiIterative();
        void solveEikonalPoint(int,int);

        void neighbourAvg(fluid& f, int, int);
        void fixFreshlyCleared();


        void run();
        double calcSoundSpeed(std::array<double,4>);

        std::array<double,4> exactRiemann(std::array<double,4>,std::array<double,4>, EOS*);

        // slope limiting
        void calcHalfSlopes(fluid&);
        void calcr(fluid&);
        std::array<double,4> static calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight);
        void calcUBars(fluid&);
        void updUBars(fluid&, EOS*);
        std::function<std::array<double,4>(std::array<double,4>)> slopeLim;
        
        std::array<double,4> minbee(std::array<double,4>);
        std::array<double,4> vanLeer(std::array<double,4>);
        std::array<double,4> vanAlbada(std::array<double,4>);

        //

        void transmissiveBC(fluid&);
        void reflectiveBC(fluid&);
        void cylTransmissiveBC(fluid&);
        void phiBC();
        void scalarTransBC(std::vector< std::vector < double >>&);
        template <size_t N>
        void vectorTransBC(std::vector< std::vector < std::array<double,N> > >&);

        std::function<void(fluid&)> setBCs;

        std::function<std::array<double,4>(std::array<double,4>,EOS*)> flux;

        void setWriteInterval(double);
        bool checkWrite;
        void writeData(std::string) const;

        // getters
        int ghosts();
        double get_dx()const;
        double get_dy()const;
        fluid& get_fluid(bool);
        boundary& get_boundary();


        double dt;

        std::string dirName;
        std::string fileName;

        std::map<std::string,int> varMap;
        // constructor

        solver(double,double,double,double,double,double,int,int,int,double,double);

        // utility functions

        static void printPrec(double,int);
        static void print_arr(std::vector < std::vector< std::array<double,4> > >, int);
        void print_var(int);
        static void print_vect(std::vector < std::array<double,4> >);
        static void print_state(std::array<double,4>);
        static std::array<double,4> set_vals(double,double,double,double);
        static void resize2D(int,int,std::vector < std::vector < double >>&);
        static void resize2Db(int,int,std::vector < std::vector < std::array<double,2> >>&);
        static void resize2Dc(int,int,std::vector < std::vector < std::array<double,4> >>&);
        static void printScalarField(const std::vector < std::vector < double > >&);
        void printInterfaceArrays(std::vector<std::vector<double>>, fluid&, const std::string&);

        // Source term modelling

        void sourceUpdate();
        void sourceTerm(const std::vector < std::vector< std::array<double,4> > > uVals);
 
        // EOS handling

        void setEOS(std::array<EOS*,2> eosInstances){
            eos = eosInstances;
        };

        std::array<EOS*,2> eos = {nullptr,nullptr};

    private:

        double dx;
        double dy;
        int nCellsX;
        int nCellsY;

        int nGhost;
        double cSound;

        fluid fluid1, fluid2;

        boundary bound;

        // level set
        std::vector< std::vector < double > > phi;
        std::vector< std::vector < double > > phiPlus1;
        std::vector< std::vector < double > > phiOld;
        std::vector<std::array<int, 2>> freshlyCleared;

        std::vector< std::vector < double > > phiGrads;



        double x0;
        double x1;
        double y0;
        double y1;
        double startTime;
        double endTime;
        double time;
        double slopeWeight;
        int splitFlip = 0;

        double dtReducer;
        int maxIter;

        double cour;

        double writeInterval;

        std::vector<double> errs;

        const double gamma;
};

/*
template<size_t N>
std::array<double, N> operator-(const std::array<double, N>&, const std::array<double, N>&);

template<size_t N>
std::array<double, N> operator+(const std::array<double, N>&, const std::array<double, N>&);

template<size_t N>
std::array<double, N> operator*(const double,const std::array<double, N>&);

template<size_t N>
std::array<double, N> operator*(const std::array<double, N>&,const std::array<double, N>&);

template<size_t N>
std::array<double, N> operator/(const std::array<double, N>&, const double);

template<size_t N>
std::array<double, N> elementDivide(const std::array<double, N>&, const std::array<double, N>&);
*/

inline void resize2D(int szY, int szX, std::vector < std::vector < std::array<double,4> > >& arr){
    arr.resize(szY);
    for (int j = 0; j < szY; ++j){
        arr[j].resize(szX);
    }
}



