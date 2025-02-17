#include <vector>
#include <string>
#include <functional>
#include <map>
#include "EOS.H"
#include "riem.H"
#include "fluid.H"

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
        std::vector< std::array< int,2 > > interfaceCells;
        void findBoundary();
        std::vector < std::array<double,2> > interfaceNormals;
        void calcPhiGrads();
        std::vector < std::array < std::array<double,2>,2 > > interfacePositions;
        void calcInterface();
        std::vector <  std::array< std::array<double,4>,2 > > interfaceStates;
        void interpInterfaceStates();
        std::array<double,4> bilinear(double,double,fluid&);
        void interfaceRiem();
        std::vector < std::array<double,4> > riemInterfaceStates;
        void resolveVelocities();
        std::vector < std::array < std::array<double,2>,2> > resolvedVelocities;
        void calcStarStates();
        std::vector < std::array < std::array<double,4>,2> > starStates;
        void calcPhiNormals();
        std::vector < std::vector < std::array<double,2 > > > phiNormals;
        void calcnDotPhiNormals();
        std::vector < std::vector < std::array < double,4 > > > uExtrap;
        std::vector < std::vector < std::array < double,4 > > > nDotGradPhi;
        void setInterface();
        void setGhostFluids();
        void eikonalDt();
        double extrapDt;


        void setDt();
        void laxFriedrichs();
        std::array<double,4> fBurgersLF(std::array<double,4>)const;
        std::array<double,4> finvectLF(std::array<double,4>)const;
        std::array<double,4> fEuler(std::array<double,4>,EOS*);

        void SLIC(fluid&, EOS*);
        void MUSCL(fluid&, EOS*);
        std::array<double,4> laxFhalf(const std::array<double,4>, const std::array<double,4>, EOS*)const;
        std::array<double,4> LF(const std::array<double,4>, const std::array<double,4>, EOS*)const;
        void pointsUpdate(fluid&);
        void phiUpdate();
        void run();
        double calcSoundSpeed(std::array<double,4>);

        std::array<double,4> riemannSolver(std::array<double,4>,std::array<double,4>);

        // slope limiting
        void calcHalfSlopes(fluid&);
        void calcr(fluid&);
        std::array<double,4> static calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight);
        void calcUBars(fluid&);
        void updUBars(fluid&, EOS*);
        std::function<std::array<double,4>(std::array<double,4>)> slopeLim;
        
        std::array<double,4> minbee(std::array<double,4>);
        std::array<double,4> vanLeer(std::array<double,4>);

        //

        void transmissiveBC(fluid&);
        void cylTransmissiveBC(fluid&);
        void phiBC();

        std::function<std::array<double,4>(std::array<double,4>,EOS*)> flux;

        void setWriteInterval(double);
        bool checkWrite;
        void writeData(std::string) const;

        // getters
        int ghosts();
        double get_dx()const;
        double get_dy()const;
        fluid& get_fluid(bool);


        double dt;

        std::string dirName;
        std::string fileName;

        std::map<std::string,int> varMap;
        // constructor

        solver(double,double,double,double,double,double,int,int,double,double);

        // utility functions

        static void print_arr(std::vector < std::vector< std::array<double,4> > >, int);
        static void print_vect(std::vector < std::array<double,4> >);
        static void print_state(std::array<double,4>);
        static std::array<double,4> set_vals(double,double,double,double);
        static void resize2D(int,int,std::vector < std::vector < double >>&);
        static void resize2Db(int,int,std::vector < std::vector < std::array<double,2> >>&);
        static void resize2Dc(int,int,std::vector < std::vector < std::array<double,4> >>&);

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
        int nCells;
        int nGhost;
        double cSound;

        fluid fluid1, fluid2; // Add these lines

        // level set
        std::vector< std::vector < double > > phi;
        std::vector< std::vector < double > > phiPlus1;



        double x0;
        double x1;
        double y0;
        double y1;
        double startTime;
        double endTime;
        double time;
        double slopeWeight;

        double cour;

        double writeInterval;

        std::vector<double> errs;

        const double gamma;
};

std::array<double, 4> operator-(const std::array<double, 4>&, const std::array<double, 4>&);

std::array<double, 4> operator+(const std::array<double, 4>&, const std::array<double, 4>&);

std::array<double, 4> operator*(const double,const std::array<double, 4>&);

std::array<double, 4> operator*(const std::array<double, 4>&,const std::array<double, 4>&);

std::array<double, 4> operator/(const std::array<double, 4>&, const double);

std::array<double, 4> elementDivide(const std::array<double, 4>&, const std::array<double, 4>&);


inline void resize2D(int szY, int szX, std::vector < std::vector < std::array<double,4> > >& arr){
    arr.resize(szY);
    for (int j = 0; j < szY; ++j){
        arr[j].resize(szX);
    }
}



