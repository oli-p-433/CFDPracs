#include <vector>
#include <string>
#include <functional>
#include <map>
#include "fluid.H"

class solverGFM{

    public:

        void init(const std::vector< std::array<double,3> >& init, fluid* f); // sets the initial values for the fluid
        void phiInit(std::vector<double>); // sets initial values for phi

        void setDt();

        std::array<double,3> fEuler(const EOS* eos, std::array<double,3>); // euler flux calculation
        std::function<std::array<double,3>(const EOS* eos, std::array<double,3>)> flux; // placeholder function for changing flux
        
        std::array<double,3> laxFhalf(const EOS* eos, const std::array<double,3>, const std::array<double,3>)const;
        std::array<double,3> LF(const EOS* eos, const std::array<double,3>, const std::array<double,3>)const;
 
        // slope limited scheme
        void SLIC();
        void calcHalfSlopes(fluid*);
        void calcr(fluid*);
        std::array<double,3> static calcSlope(double omega, std::array<double,3> slopeleft, std::array<double,3> sloperight);
        void calcUBars(fluid*);
        void updUBars(fluid*);
        std::function<std::array<double,3>(std::array<double,3>)> slopeLim;
        
        std::array<double,3> minbee(std::array<double,3>);
        std::array<double,3> superbee(std::array<double,3>);
        std::array<double,3> vanLeer(std::array<double,3>);

        void pointsUpdate();
        void levelSetUpdate();
        void run();

        void transmissiveBC(); // sets transmissive BCs

        void setWriteInterval(double);
        bool checkWrite;
        void writeData(std::string) const;

        int ghosts();
        double get_dx()const;


        double dt;

        std::string dirName;
        std::string fileName;

        std::map<std::string,int> varMap;

        // fluid handling
        void setFluids(std::vector<fluid*> fluidsVect){
            fluids = fluidsVect;
        };

        // constructor
        solverGFM(double,double,double,double,int,int,double);


    protected:
        double x0;
        double x1;
        double startTime;
        double endTime;
        double time;

        double cour;

        double writeInterval;
        double timeMulti;

        double dx;
        int nCells;
        int nGhost;

        std::vector<fluid*> fluids = {nullptr};
        std::vector<std::array<double,3>> temp_u;

        std::vector< double > phi;
        std::vector< double > phiPlus1;




};

std::array<double, 3> operator-(const std::array<double, 3>&, const std::array<double, 3>&);

std::array<double, 3> operator+(const std::array<double, 3>&, const std::array<double, 3>&);

std::array<double, 3> operator*(const double,const std::array<double, 3>&);

std::array<double, 3> operator*(const std::array<double, 3>&,const std::array<double, 3>&);

std::array<double, 3> operator/(const std::array<double, 3>&, const double);

std::array<double, 3> elementDivide(const std::array<double, 3>&, const std::array<double, 3>&);