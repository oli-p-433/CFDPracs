#include <vector>
#include <string>
#include <functional>
#include <map>
#include "EOS.H"

class fluid{
    public:

        std::vector< std::array<double,3> > u; //indexed as u[i][j] (ith array, jth row)
        std::vector< std::array<double,3> > fluxes;

        std::vector< std::array<double,3> > halfSlopes;
        std::vector< std::array<double,3> > r;
        std::vector< std::array<double,3> > uBarL;
        std::vector< std::array<double,3> > uBarR;
        std::vector< std::array<double,3> > uBarLupd;
        std::vector< std::array<double,3> > uBarRupd;

        static void print_arr(std::vector< std::array<double,3> >, int);
        static void print_vect(std::array<double,3>);

        // EOS handling:
        void setEOS(EOS* eosInstance){
        eos = eosInstance;
        };

        EOS* eos = nullptr;

        fluid(int,int); // ncells, nghost

        // destructor
        virtual ~fluid() = default;

        // getters
        const double& get_slopeWeight() const { return slopeWeight; };

    protected:
        double cSound;

        double slopeWeight;

};