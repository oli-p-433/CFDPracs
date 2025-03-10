#include <iostream>
#include <vector>
#include <array>

class fluid{
    public:
        std::vector< std::vector < std::array<double,5> > > u, uPlus1, uPrev, uPrim, fluxes;
        std::vector< std::vector < double > > sStars;

        std::vector< std::vector< std::array<double,5> > > halfSlopes;
        std::vector< std::vector< std::array<double,5> > > r;
        std::vector< std::vector< std::array<double,5> > > uBarL;
        std::vector< std::vector< std::array<double,5> > > uBarR;
        std::vector< std::vector< std::array<double,5> > > uBarLupd;
        std::vector< std::vector< std::array<double,5> > > uBarRupd;

        fluid(int,int,int);


    private:
        int nCellsX, nCellsY, nGhost;


};