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

        std::vector< std::vector< std::array<double,4> > > uPlus1;

        void init(std::vector< std::vector< std::array<double,4> > >);

        void setDt();
        void laxFriedrichs();
        std::array<double,4> fBurgersLF(std::array<double,4>)const;
        std::array<double,4> finvectLF(std::array<double,4>)const;
        std::array<double,4> fEuler(std::array<double,4>,EOS*);
        void Richt();
        void FORCE();
        void SLIC();
        void MUSCL();
        std::array<double,4> laxFhalf(const std::array<double,4>, const std::array<double,4>)const;
        std::array<double,4> LF(const std::array<double,4>, const std::array<double,4>)const;
        void pointsUpdate();
        void run();
        double calcSoundSpeed(std::array<double,4>);

        std::array<double,4> riemannSolver(std::array<double,4>,std::array<double,4>);

        // slope limiting
        void calcHalfSlopes();
        void calcr();
        std::array<double,4> static calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight);
        void calcUBars();
        void updUBars();
        std::function<std::array<double,4>(std::array<double,4>)> slopeLim;
        
        std::array<double,4> minbee(std::array<double,4>);
        std::array<double,4> vanLeer(std::array<double,4>);

        //

        void transmissiveBC();
        void cylTransmissiveBC();

        std::function<std::array<double,4>(std::array<double,4>,EOS*)> flux;

        void setWriteInterval(double);
        bool checkWrite;
        void writeData(std::string) const;

        int ghosts();
        double get_dx()const;
        double get_dy()const;


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

        std::vector< std::vector< std::array<double,4> > > u; //indexed as u[i][j] (ith row, jth column)
        std::vector< std::vector< std::array<double,4> > > fluxesX;
        std::vector< std::vector< std::array<double,4> > > fluxesY;

        std::vector< std::vector< std::array<double,4> > > halfSlopesX;
        std::vector< std::vector< std::array<double,4> > > halfSlopesY;
        std::vector< std::vector< std::array<double,4> > > rX;
        std::vector< std::vector< std::array<double,4> > > rY;
        std::vector< std::vector< std::array<double,4> > > uBarLX;
        std::vector< std::vector< std::array<double,4> > > uBarLY;
        std::vector< std::vector< std::array<double,4> > > uBarRX;
        std::vector< std::vector< std::array<double,4> > > uBarRY;
        std::vector< std::vector< std::array<double,4> > > uBarLupdX;
        std::vector< std::vector< std::array<double,4> > > uBarLupdY;
        std::vector< std::vector< std::array<double,4> > > uBarRupdX;
        std::vector< std::vector< std::array<double,4> > > uBarRupdY;

        std::vector< std::vector< std::array<double,4> > > sourceResult;


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



