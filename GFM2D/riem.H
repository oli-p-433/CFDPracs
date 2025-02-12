#include <vector>
#include <string>
#include <functional>
#include <map>

class riemann{
    public:
        bool direction;
        // exact Riemann solver
        double calcSoundSpeed(std::array<double,3>,bool);

        std::array<double,4> exctRiemann();
        double fRiemann(double pStar);
        double fShock(double, std::array<double,3>, bool);
        double fRare(double, std::array<double,3>, bool);
        double fPrm(double, std::array<double,3>, bool);
        std::array<double,4> solveRiem(double, std::array<double,3>, bool);
        std::array<double,3> rareFan(std::array<double,3>, bool, double);

        std::array<std::string,5> waveSignature(double);
        std::array<double,5> wavePositions(std::array<std::string,5>, double, std::array<double,4>, std::array<double,4>);
        std::array<double,5> wavePos;
        
        void exactRiemannSolution(double, double);

        double gamma(bool);
        double pInf(bool);

        //

        std::string dirName;

        std::map<std::string,int> varMap;

        void writeData(std::string, std::vector<std::array<double,3>>) const;

        // converters

        std::array<double,3> primToConsv(std::array<double,3>,bool)const;
        std::array<double,3> consvToPrim(std::array<double,3>,bool)const;



        // Constructor

        riemann(double,double,std::array<double,4>,std::array<double,4>,bool,double,double,double,double,int,double pInfLeft=0, double pInfRight = 0);


    private:

        double gammaL;
        double gammaR;

        std::array<double,3> state1;
        std::array<double,3> state2;

        std::array<double,4> state1orig;
        std::array<double,4> state2orig;

        double x0;
        double x1;
        double x_disc;
        double endTime;

        double dx;
        int nCells;

        double pInfLeft;
        double pInfRight;

        // riemann problem solution
        std::vector< std::array<double,3> > uRiemExact;



};
