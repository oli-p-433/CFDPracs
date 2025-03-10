#include "fluid.H"
#include "operators.H"

fluid::fluid(int nX, int nY, int nG)
    : nCellsX(nX), nCellsY(nY), nGhost(nG){
        resizeVector(u,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(uPlus1,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(uPrev,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(uPrim,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(sStars,nCellsY,nCellsX+1);
        resizeVector(fluxes,nCellsY,nCellsX+1);
        resizeVector(halfSlopes,nCellsY,nCellsX+3);
        resizeVector(r,nCellsY,nCellsX+2);
        resizeVector(uBarL,nCellsY,nCellsX+2);
        resizeVector(uBarR,nCellsY,nCellsX+2);
        resizeVector(uBarLupd,nCellsY,nCellsX+2);
        resizeVector(uBarRupd,nCellsY,nCellsX+2);

}