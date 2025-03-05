#include <vector>
#include <array>
#include <iostream>

class fluid {
    public:
        fluid(int n_CellsX, int nCellsY, int n_Ghosts);
        std::vector<std::vector<std::array<double, 4>>> u, uPlus1, fluxesX, fluxesY, halfSlopesX, halfSlopesY, rX, rY, uBarLX, uBarLY, uBarRX, uBarRY, uBarLupdX, uBarLupdY, uBarRupdX, uBarRupdY, sourceResult, nDotGradPhi;
        
        static void resizef2D(size_t rows, size_t cols, std::vector<std::vector<std::array<double, 4>>>& vec);
        std::vector< std::array< int,2 > > interfaceCells;
        std::vector < std::array<double,2> > interfaceNormals;
        std::vector < std::array < std::array<double,2>,2 > > interfacePositions;
        std::vector <  std::array< std::array<double,4>,2 > > interfaceStates;
        std::vector < std::array<double,4> > riemInterfaceStates;
        std::vector < std::array < std::array<double,2>,2> > resolvedVelocities;
        std::vector < std::array<double,4> > starStates;
        std::vector < std::vector < std::array < double,4 > > > uExtrap;
        std::vector < double > interfacePhis;

    
    private:
        int nCellsX, nCellsY, nGhost;
    

};
