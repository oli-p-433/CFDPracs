#include "fluid.H" 

fluid::fluid(int n_CellsX, int n_CellsY, int n_Ghosts)
    : nCellsX(n_CellsX), nCellsY(n_CellsY), nGhost(n_Ghosts) {

    resizef2D(nCellsY + 2 * nGhost, nCellsX + 2 * nGhost, u);
    resizef2D(nCellsY + 2 * nGhost, nCellsX + 2 * nGhost, uPlus1);
    resizef2D(nCellsY + 2 * nGhost, nCellsX + 2 * nGhost, nDotGradPhi);
    resizef2D(nCellsY, nCellsX + 1, fluxesX);
    resizef2D(nCellsY + 1, nCellsX, fluxesY);
    resizef2D(nCellsY + 2, nCellsX + 3, halfSlopesX);
    resizef2D(nCellsY + 3, nCellsX + 2, halfSlopesY);
    resizef2D(nCellsY + 2, nCellsX + 2, rX);
    resizef2D(nCellsY + 2, nCellsX + 2, rY);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarLX);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarLY);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarRX);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarRY);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarLupdX);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarLupdY);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarRupdX);
    resizef2D(nCellsY + 2, nCellsX + 2, uBarRupdY);
    resizef2D(nCellsY, nCellsX, sourceResult);


}

void fluid::resizef2D(size_t rows, size_t cols, std::vector<std::vector<std::array<double, 4>>>& vec) {
    vec.resize(rows);
    for (auto& row : vec) {
        row.resize(cols);
    }
}