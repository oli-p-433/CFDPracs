#include <vector>
#include <array>
#include <iostream>

class fluid {
    public:
        fluid(int n_Cells, int n_Ghosts);
        std::vector<std::vector<std::array<double, 4>>> u, uPlus1, fluxesX, fluxesY, halfSlopesX, halfSlopesY, rX, rY, uBarLX, uBarLY, uBarRX, uBarRY, uBarLupdX, uBarLupdY, uBarRupdX, uBarRupdY, sourceResult, nDotGradPhi;
        static void resizef2D(size_t rows, size_t cols, std::vector<std::vector<std::array<double, 4>>>& vec);

    
    private:
        int nCells, nGhost;
    

};
