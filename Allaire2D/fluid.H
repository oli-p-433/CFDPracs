#include <iostream>
#include <vector>
#include <array>

class fluid{
    public:
        std::vector< std::vector < std::array<double,6> > > u, uPlus1, uPrev, uPrim, fluxesX, fluxesY;
        std::vector< std::vector < double > > sStarsX, sStarsY;

        std::vector< std::vector< std::array<double,6> > > halfSlopesX, halfSlopesY;
        std::vector< std::vector< std::array<double,6> > > rX, rY;
        std::vector< std::vector< std::array<double,6> > > uBarLX, uBarLY;
        std::vector< std::vector< std::array<double,6> > > uBarRX, uBarRY;
        std::vector< std::vector< std::array<double,6> > > uBarLupdX, uBarLupdY;
        std::vector< std::vector< std::array<double,6> > > uBarRupdX, uBarRupdY;

        fluid(int,int,int);


    private:
        int nCellsX, nCellsY, nGhost;


};