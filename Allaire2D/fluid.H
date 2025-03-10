#include <iostream>
#include <vector>
#include <array>

class fluid{
    public:
        std::vector< std::vector < std::array<double,6> > > u, uPlus1, uPrev, uPrim, fluxesX, fluxesY;
        std::vector< std::vector < double > > sStarsX, sStarsY;

        std::vector< std::vector< std::array<double,6> > > halfSlopesX;
        std::vector< std::vector< std::array<double,6> > > rX;
        std::vector< std::vector< std::array<double,6> > > uBarLX;
        std::vector< std::vector< std::array<double,6> > > uBarRX;
        std::vector< std::vector< std::array<double,6> > > uBarLupdX;
        std::vector< std::vector< std::array<double,6> > > uBarRupdX;

        fluid(int,int,int);


    private:
        int nCellsX, nCellsY, nGhost;


};