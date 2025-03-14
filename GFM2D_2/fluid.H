#include <iostream>
#include <vector>
#include <array>

class fluid{
    public:
        std::vector< std::vector < std::array<double,4> > > u, uPlus1, uPrev, uPrim, fluxesX, fluxesY;
        std::vector< std::vector < double > > sStarsX, sStarsY;

        std::vector< std::vector< std::array<double,4> > > halfSlopesX, halfSlopesY;
        std::vector< std::vector< std::array<double,4> > > rX, rY;
        std::vector< std::vector< std::array<double,4> > > uBarLX, uBarLY;
        std::vector< std::vector< std::array<double,4> > > uBarRX, uBarRY;
        std::vector< std::vector< std::array<double,4> > > uBarLupdX, uBarLupdY;
        std::vector< std::vector< std::array<double,4> > > uBarRupdX, uBarRupdY;

        std::vector< std::array< int,2 > > interfaceCells;
        std::vector < std::array<double,2> > interfaceNormals;
        std::vector < std::array < std::array<double,2>,2 > > interfacePositions;
        std::vector <  std::array< std::array<double,4>,2 > > interfaceStates;
        std::vector < std::array<double,4> > riemInterfaceStates;
        std::vector < std::array < std::array<double,2>,2> > resolvedVelocities;
        std::vector < std::array<double,4> > starStates;
        std::vector < std::vector < std::array < double,4 > > > uExtrap;
        std::vector < double > interfacePhis;

        fluid(int,int,int);


    private:
        int nCellsX, nCellsY, nGhost;


};