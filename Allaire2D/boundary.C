#include "boundary.H"

// alpha rho1 rho2 vx vy p
// Reflective left boundary: update left ghost columns (j = 0 .. nGhost-1)
// for interior rows (i = nGhost .. nCellsY+nGhost-1).
void boundary::reflectiveLeftBC(fluid &fluid) {
    for (int i = nGhost; i < nCellsY + nGhost; ++i) {
        for (int j = 0; j < nGhost; ++j) {
            // Copy interior state from first interior column (j = nGhost)
            fluid.u[i][j] = fluid.u[i][nGhost];
            fluid.uPlus1[i][j] = fluid.u[i][nGhost];
            // Flip the x momentum (index 1)
            fluid.u[i][j][3] = -fluid.u[i][nGhost][3];
            fluid.uPlus1[i][j][3] = -fluid.u[i][nGhost][3];
        }
    }
}

// Reflective right boundary: update right ghost columns 
// (j = nCellsX+nGhost .. nCellsX+2*nGhost-1) for interior rows.
void boundary::reflectiveRightBC(fluid &fluid) {
    for (int i = nGhost; i < nCellsY + nGhost; ++i) {
        for (int j = nCellsX + nGhost; j < nCellsX + 2 * nGhost; ++j) {
            // Copy interior state from last interior column (j = nCellsX+nGhost-1)
            fluid.u[i][j] = fluid.u[i][nCellsX + nGhost - 1];
            fluid.uPlus1[i][j] = fluid.u[i][nCellsX + nGhost - 1];
            // Flip the x momentum (index 1)
            fluid.u[i][j][3] = -fluid.u[i][nCellsX + nGhost - 1][3];
            fluid.uPlus1[i][j][3] = -fluid.u[i][nCellsX + nGhost - 1][3];
        }
    }
}

// Reflective bottom boundary: update bottom ghost rows 
// (i = 0 .. nGhost-1) for interior columns.
void boundary::reflectiveBottomBC(fluid &fluid) {
    for (int i = 0; i < nGhost; ++i) {
        for (int j = nGhost; j < nCellsX + nGhost; ++j) {
            // Copy interior state from first interior row (i = nGhost)
            fluid.u[i][j] = fluid.u[nGhost][j];
            fluid.uPlus1[i][j] = fluid.u[nGhost][j];
            // Flip the y momentum (index 2)
            fluid.u[i][j][4] = -fluid.u[nGhost][j][4];
            fluid.uPlus1[i][j][4] = -fluid.u[nGhost][j][4];
        }
    }
}

// Reflective top boundary: update top ghost rows 
// (i = nCellsY+nGhost .. nCellsY+2*nGhost-1) for interior columns.
void boundary::reflectiveTopBC(fluid &fluid) {
    for (int i = nCellsY + nGhost; i < nCellsY + 2 * nGhost; ++i) {
        for (int j = nGhost; j < nCellsX + nGhost; ++j) {
            // Copy interior state from last interior row (i = nCellsY+nGhost-1)
            fluid.u[i][j] = fluid.u[nCellsY + nGhost - 1][j];
            fluid.uPlus1[i][j] = fluid.u[nCellsY + nGhost - 1][j];
            // Flip the y momentum (index 2)
            fluid.u[i][j][4] = -fluid.u[nCellsY + nGhost - 1][j][4];
            fluid.uPlus1[i][j][4] = -fluid.u[nCellsY + nGhost - 1][j][4];
        }
    }
}

// Transmissive left boundary: update left ghost columns (j = 0 .. nGhost-1)
// for interior rows (i = nGhost .. nCellsY+nGhost-1).
void boundary::transmissiveLeftBC(fluid &fluid) {
    for (int i = nGhost; i < nCellsY + nGhost; ++i) {
        for (int j = 0; j < nGhost; ++j) {
            // Copy the interior state from the first interior column (j = nGhost)
            fluid.u[i][j] = fluid.u[i][nGhost];
            fluid.uPlus1[i][j] = fluid.u[i][nGhost];
        }
    }
}

// Transmissive right boundary: update right ghost columns 
// (j = nCellsX+nGhost .. nCellsX+2*nGhost-1) for interior rows.
void boundary::transmissiveRightBC(fluid &fluid) {
    for (int i = nGhost; i < nCellsY + nGhost; ++i) {
        for (int j = nCellsX + nGhost; j < nCellsX + 2 * nGhost; ++j) {
            // Copy the interior state from the last interior column (j = nCellsX+nGhost-1)
            fluid.u[i][j] = fluid.u[i][nCellsX + nGhost - 1];
            fluid.uPlus1[i][j] = fluid.u[i][nCellsX + nGhost - 1];
        }
    }
}

// Transmissive bottom boundary: update bottom ghost rows 
// (i = 0 .. nGhost-1) for interior columns.
void boundary::transmissiveBottomBC(fluid &fluid) {
    for (int i = 0; i < nGhost; ++i) {
        for (int j = nGhost; j < nCellsX + nGhost; ++j) {
            // Copy the interior state from the first interior row (i = nGhost)
            fluid.u[i][j] = fluid.u[nGhost][j];
            fluid.uPlus1[i][j] = fluid.u[nGhost][j];
        }
    }
}

// Transmissive top boundary: update top ghost rows 
// (i = nCellsY+nGhost .. nCellsY+2*nGhost-1) for interior columns.
void boundary::transmissiveTopBC(fluid &fluid) {
    for (int i = nCellsY + nGhost; i < nCellsY + 2 * nGhost; ++i) {
        for (int j = nGhost; j < nCellsX + nGhost; ++j) {
            // Copy the interior state from the last interior row (i = nCellsY+nGhost-1)
            fluid.u[i][j] = fluid.u[nCellsY + nGhost - 1][j];
            fluid.uPlus1[i][j] = fluid.u[nCellsY + nGhost - 1][j];
        }
    }
}

// Update the bottom-left corner ghost cells (i = 0..nGhost-1, j = 0..nGhost-1)
// The behavior depends on whether the left and bottom boundaries are reflective.
// 'leftReflective' and 'bottomReflective' indicate the type for that wall.
void boundary::updateBottomLeftCorner(fluid &fluid, bool leftReflective, bool bottomReflective) {
    // Use the interior cell at (nGhost, nGhost) as the base.
    // (You might choose a different cell if that better reflects your grid layout.)
    for (int i = 0; i < nGhost; ++i) {
        for (int j = 0; j < nGhost; ++j) {
            fluid.u[i][j] = fluid.u[nGhost][nGhost];
            fluid.uPlus1[i][j] = fluid.u[nGhost][nGhost];
            
            // If the left wall is reflective, flip the x momentum (index 1)
            if (leftReflective) {
                fluid.u[i][j][3] = -fluid.u[nGhost][nGhost][3];
                fluid.uPlus1[i][j][3] = -fluid.u[nGhost][nGhost][3];
            }
            // If the bottom wall is reflective, flip the y momentum (index 2)
            if (bottomReflective) {
                fluid.u[i][j][4] = -fluid.u[nGhost][nGhost][4];
                fluid.uPlus1[i][j][4] = -fluid.u[nGhost][nGhost][4];
            }
        }
    }
}

void boundary::updateTopRightCorner(fluid &fluid, bool rightReflective, bool topReflective) {
    for (int i = nCellsY + nGhost; i < nCellsY + 2 * nGhost; ++i) {
        for (int j = nCellsX + nGhost; j < nCellsX + 2 * nGhost; ++j) {
            fluid.u[i][j] = fluid.u[nCellsY + nGhost - 1][nCellsX + nGhost - 1];
            fluid.uPlus1[i][j] = fluid.u[nCellsY + nGhost - 1][nCellsX + nGhost - 1];
            
            if (rightReflective) {
                fluid.u[i][j][3] = -fluid.u[nCellsY + nGhost - 1][nCellsX + nGhost - 1][3];
                fluid.uPlus1[i][j][3] = -fluid.u[nCellsY + nGhost - 1][nCellsX + nGhost - 1][3];
            }
            if (topReflective) {
                fluid.u[i][j][4] = -fluid.u[nCellsY + nGhost - 1][nCellsX + nGhost - 1][4];
                fluid.uPlus1[i][j][4] = -fluid.u[nCellsY + nGhost - 1][nCellsX + nGhost - 1][4];
            }
        }
    }
}

// Reflective bottom-right corner: ghost cells with
// i = 0 .. nGhost-1 and j = nCellsX+nGhost .. nCellsX+2*nGhost-1.
// 'rightReflective' and 'bottomReflective' indicate whether the right
// and bottom boundaries are reflective.
void boundary::updateBottomRightCorner(fluid &fluid, bool rightReflective, bool bottomReflective) {
    for (int i = 0; i < nGhost; ++i) {
        for (int j = nCellsX + nGhost; j < nCellsX + 2 * nGhost; ++j) {
            // Use the interior cell at (nGhost, nCellsX+nGhost-1) as the base.
            fluid.u[i][j] = fluid.u[nGhost][nCellsX + nGhost - 1];
            fluid.uPlus1[i][j] = fluid.u[nGhost][nCellsX + nGhost - 1];
            // If the right wall is reflective, flip the x momentum (index 1)
            if (rightReflective) {
                fluid.u[i][j][3] = -fluid.u[nGhost][nCellsX + nGhost - 1][3];
                fluid.uPlus1[i][j][3] = -fluid.u[nGhost][nCellsX + nGhost - 1][3];
            }
            // If the bottom wall is reflective, flip the y momentum (index 2)
            if (bottomReflective) {
                fluid.u[i][j][4] = -fluid.u[nGhost][nCellsX + nGhost - 1][4];
                fluid.uPlus1[i][j][4] = -fluid.u[nGhost][nCellsX + nGhost - 1][4];
            }
        }
    }
}

// Reflective top-left corner: ghost cells with
// i = nCellsY+nGhost .. nCellsY+2*nGhost-1 and j = 0 .. nGhost-1.
// 'leftReflective' and 'topReflective' indicate whether the left
// and top boundaries are reflective.
void boundary::updateTopLeftCorner(fluid &fluid, bool leftReflective, bool topReflective) {
    for (int i = nCellsY + nGhost; i < nCellsY + 2 * nGhost; ++i) {
        for (int j = 0; j < nGhost; ++j) {
            // Use the interior cell at (nCellsY+nGhost-1, nGhost) as the base.
            fluid.u[i][j] = fluid.u[nCellsY + nGhost - 1][nGhost];
            fluid.uPlus1[i][j] = fluid.u[nCellsY + nGhost - 1][nGhost];
            // If the left wall is reflective, flip the x momentum (index 1)
            if (leftReflective) {
                fluid.u[i][j][3] = -fluid.u[nCellsY + nGhost - 1][nGhost][3];
                fluid.uPlus1[i][j][3] = -fluid.u[nCellsY + nGhost - 1][nGhost][3];
            }
            // If the top wall is reflective, flip the y momentum (index 2)
            if (topReflective) {
                fluid.u[i][j][4] = -fluid.u[nCellsY + nGhost - 1][nGhost][4];
                fluid.uPlus1[i][j][4] = -fluid.u[nCellsY + nGhost - 1][nGhost][4];
            }
        }
    }
}


boundary::boundary(int nx, int ny, int ng) : nCellsX(nx), nCellsY(ny), nGhost(ng) {
}

