#include "fluid.H" 

fluid::fluid(int n_Cells, int n_Ghosts)
    : nCells(n_Cells), nGhost(n_Ghosts) {

    resizef2D(nCells + 2 * nGhost, nCells + 2 * nGhost, u);
    resizef2D(nCells + 2 * nGhost, nCells + 2 * nGhost, uPlus1);
    resizef2D(nCells + 2, nCells + 1, fluxesX);
    resizef2D(nCells + 1, nCells + 2, fluxesY);
    resizef2D(nCells + 2, nCells + 3, halfSlopesX);
    resizef2D(nCells + 3, nCells + 2, halfSlopesY);
    resizef2D(nCells + 2, nCells + 2, rX);
    resizef2D(nCells + 2, nCells + 2, rY);
    resizef2D(nCells + 2, nCells + 2, uBarLX);
    resizef2D(nCells + 2, nCells + 2, uBarLY);
    resizef2D(nCells + 2, nCells + 2, uBarRX);
    resizef2D(nCells + 2, nCells + 2, uBarRY);
    resizef2D(nCells + 2, nCells + 2, uBarLupdX);
    resizef2D(nCells + 2, nCells + 2, uBarLupdY);
    resizef2D(nCells + 2, nCells + 2, uBarRupdX);
    resizef2D(nCells + 2, nCells + 2, uBarRupdY);
    resizef2D(nCells, nCells, sourceResult);
}

void fluid::resizef2D(size_t rows, size_t cols, std::vector<std::vector<std::array<double, 4>>>& vec) {
    vec.resize(rows);
    for (auto& row : vec) {
        row.resize(cols);
    }
}