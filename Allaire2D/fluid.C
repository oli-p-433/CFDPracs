#include "fluid.H"
#include "operators.H"

fluid::fluid(int nX, int nY, int nG)
    : nCellsX(nX), nCellsY(nY), nGhost(nG){
        resizeVector(u,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(uPlus1,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(uPrev,nCellsY+2*nGhost,nCellsX+2*nGhost);
        resizeVector(uPrim,nCellsY+2*nGhost,nCellsX+2*nGhost);

        resizeVector(sStarsX,nCellsY,nCellsX+1);
        resizeVector(fluxesX,nCellsY,nCellsX+1);
        resizeVector(halfSlopesX,nCellsY,nCellsX+3);
        resizeVector(rX,nCellsY,nCellsX+2);
        resizeVector(uBarLX,nCellsY,nCellsX+2);
        resizeVector(uBarRX,nCellsY,nCellsX+2);
        resizeVector(uBarLupdX,nCellsY,nCellsX+2);
        resizeVector(uBarRupdX,nCellsY,nCellsX+2);

        resizeVector(sStarsY,nCellsY+1,nCellsX);
        resizeVector(fluxesY,nCellsY+1,nCellsX);
        resizeVector(halfSlopesY,nCellsY+3,nCellsX);
        resizeVector(rY,nCellsY+2,nCellsX);
        resizeVector(uBarLY,nCellsY+2,nCellsX);
        resizeVector(uBarRY,nCellsY+2,nCellsX);
        resizeVector(uBarLupdY,nCellsY+2,nCellsX);
        resizeVector(uBarRupdY,nCellsY+2,nCellsX);
        

}