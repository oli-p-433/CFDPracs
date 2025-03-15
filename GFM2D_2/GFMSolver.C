#include "GFMSolver.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <omp.h>

void solver::run(){
    alpha = 0;
    std::cout << "pre-bc" << std::endl;

    setBCs(fluid1);
    setBCs(fluid2);
    phiBC();
    std::cout << "set BCs" << std::endl;

    splitFlip = 0;

    do{

                //time = 0+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        
        phiBC();
        //std::cout << "finding boundary" << std::endl;
        this->findBoundary();
        //std::cout << "number of interface cells = " << interfaceCells.size() << std::endl;
        //std::cout << "boundary found" << std::endl;

        //printInterfaceArray(phi, interfaceCells, "interface.txt");

        this->calcInterfaceNormals(fluid1);
        this->calcInterfaceNormals(fluid2);
        //std::cout << "phi grads calculated" << std::endl;
        this->calcInterface(fluid1);
        this->calcInterface(fluid2);
        //std::cout << "interface calculated" << std::endl;

        //time = 0;
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        
        this->interpInterfaceStates(fluid1);
        this->interpInterfaceStates(fluid2);
        //std::cout << "interface states interpolated" << std::endl;

        this->resolveVelocities(fluid1);
        this->resolveVelocities(fluid2);
        //std::cout << "resolved velocities" << std::endl;
        this->interfaceRiem(fluid1);
        this->interfaceRiem(fluid2);
        //std::cout << "interface riem" << std::endl;
        this->calcStarStates(fluid1);
        this->calcStarStates(fluid2);
        //std::cout << "star states calcd" << std::endl;

        this->GFFastSweeping();
        //fastSweepExtrapolation(fluid1.u,fluid1.interfaceCells,phi,fluid1,1000);
        //fastSweepExtrapolation(fluid2.u,fluid2.interfaceCells,phi,fluid2,1000);
        //this->setGhostFluids();
        //std::cout << "set ghost fluids" << std::endl; 
        
        this->setDt();
        double origDt = dt;
        
        //realTime = time;

        std::cout << "time is " << time << std::endl;

        //std::cout << "setting phiBC" << std::endl;
        phiBC();
        //std::cout << "updating phi" << std::endl;

        //time = 3+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        
        phiUpdate(splitFlip);
        //std::cout << "phi updated" << std::endl;

        //time = 4+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        fixFreshlyCleared();
        //std::cout << "freshly cleared fixed" << std::endl;

        //time = 5+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        //time = 0.0;
        
        if (splitFlip % 1 == 0){
           phiBC();
           reinitialiseLevelSet(phi,2);
           //std::cout << "phi reinitialised" << std::endl;
        }
        
        splitFlip++;
        if (splitFlip%2 == 0){
            //this->fluxMethod(XDIR,fluid1,eos[0]); this->RK2(XDIR,fluid1,eos[0]);
            //this->fluxMethod(XDIR,fluid2,eos[1]); this->RK2(XDIR,fluid2,eos[1]);
            //this->fluxMethod(YDIR,fluid1,eos[0]); this->RK2(YDIR,fluid1,eos[0]);
            //this->fluxMethod(YDIR,fluid2,eos[1]); this->RK2(YDIR,fluid2,eos[1]);

            this->fluxMethod(XDIR,fluid1,eos[0]); this->RK2(XDIR,fluid1,eos[0]);
            this->fluxMethod(XDIR,fluid2,eos[1]); this->RK2(XDIR,fluid2,eos[1]);
            this->fluxMethod(YDIR,fluid1,eos[0]); this->RK2(YDIR,fluid1,eos[0]);
            this->fluxMethod(YDIR,fluid2,eos[1]); this->RK2(YDIR,fluid2,eos[1]);

        } else {
            //this->fluxMethod(YDIR,fluid1,eos[0]); this->RK2(YDIR,fluid1,eos[0]);
            //this->fluxMethod(YDIR,fluid2,eos[1]); this->RK2(YDIR,fluid2,eos[1]);
            //this->fluxMethod(XDIR,fluid1,eos[0]); this->RK2(XDIR,fluid1,eos[0]);
            //this->fluxMethod(XDIR,fluid2,eos[1]); this->RK2(XDIR,fluid2,eos[1]);

            this->fluxMethod(YDIR,fluid1,eos[0]); this->RK2(YDIR,fluid1,eos[0]);
            this->fluxMethod(YDIR,fluid2,eos[1]); this->RK2(YDIR,fluid2,eos[1]);
            this->fluxMethod(XDIR,fluid1,eos[0]); this->RK2(XDIR,fluid1,eos[0]);
            this->fluxMethod(XDIR,fluid2,eos[1]); this->RK2(XDIR,fluid2,eos[1]);

        }

        setBCs(fluid1);
        setBCs(fluid2);

        // empty interface cell buffer
        this->fluid1.interfaceCells={};
        this->fluid2.interfaceCells={};
        this->freshlyCleared={};

        if (checkWrite==1){
            std::cout << "------------------ writing " << time << " ------------------"<< std::endl;
            writeData();
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData();
};

void solver::findBoundary(){ // checks for change in sign of level set. emplaces indices of cell with phi <= 0.
    for (size_t i = nG; i < phi.size()-nG; ++i) {  // Loop through all cells (except last)
        for (size_t j = nG; j < phi[0].size()-nG; ++j){
            //if (phi[i][j] == 0.0){
            //    std::array<int,2> cell = {i,j};
            //    fluid2.interfaceCells.push_back(cell); // 1
            //}
            if ((std::signbit(phi[i][j]) != std::signbit(phi[i+1][j])) || 
                (std::signbit(phi[i][j]) != std::signbit(phi[i][j+1]))){  // Sign change detected
                //std::cout << phi[i] << " " << phi[i+1] << std::endl;
                if (phi[i][j] == 0.0){
                    std::array<int,2> cell = {i,j};
                    fluid2.interfaceCells.push_back(cell); // 1
                }
                if (std::signbit(phi[i][j]) == 1 && phi[i][j] != 0.0){ // phi_ij negative
                    std::array<int,2> cell = {i,j};
                    fluid2.interfaceCells.push_back(cell);
                    if (std::signbit(phi[i+1][j]) == 0){
                        std::array<int,2> cell = {i+1,j};
                        fluid1.interfaceCells.push_back(cell);
                    }
                    if (std::signbit(phi[i][j+1]) == 0){
                        std::array<int,2> cell = {i,j+1};
                        fluid1.interfaceCells.push_back(cell);
                    }
                } else if (std::signbit(phi[i][j]) == 0 && phi[i][j] != 0.0){ // phi_ij +ve
                    std::array<int,2> cell = {i,j};
                    fluid1.interfaceCells.push_back(cell);
                    if (std::signbit(phi[i+1][j]) == 1){
                        std::array<int,2> cell = {i+1,j};
                        fluid2.interfaceCells.push_back(cell); //1
                    }
                    if (std::signbit(phi[i][j+1]) == 1){
                        std::array<int,2> cell = {i,j+1};
                        fluid2.interfaceCells.push_back(cell); //1
                    }
                }
                //std::cout << "interface at cell " << i << std::endl; 
            } 
        }
    };
    int ny = phi.size()-nG-1;
    int nx = phi.size()-nG-1;
    if (phi[ny-1][nx] * phi[ny][nx] < 0){ // corner cases
        if (std::signbit(phi[ny][nx])==0){
            std::array<int,2> cell1 = {ny-1,nx};
            std::array<int,2> cell2 = {ny,nx};
            fluid2.interfaceCells.push_back(cell1); //1
            fluid1.interfaceCells.push_back(cell2);
        }
        //std::array<int,2> cell1 = (phi[n][n]>0) ? std::array<int,2>{n-1,n} : std::array<int,2>{n,n};
        //fluid1.interfaceCells.push_back(cell1);
    }
    if (phi[ny][nx-1] * phi[ny][nx] < 0){
        if (std::signbit(phi[ny][nx])==0){
            std::array<int,2> cell1 = {ny,nx-1};
            std::array<int,2> cell2 = {ny,nx};
            fluid2.interfaceCells.push_back(cell1); //1
            fluid1.interfaceCells.push_back(cell2);
        }
        //std::array<int,2> cell = (phi[n][n]>0) ? std::array<int,2>{n,n-1} : std::array<int,2>{n,n};
        //fluid1.interfaceCells.push_back(cell);
    }
    

    for (fluid& f : {std::ref(fluid1), std::ref(fluid2)}){
        f.interfaceStates.resize(f.interfaceCells.size());
        f.interfaceNormals.resize(f.interfaceCells.size());
        f.interfacePositions.resize(f.interfaceCells.size());
        f.resolvedVelocities.resize(f.interfaceCells.size());
        f.riemInterfaceStates.resize(f.interfaceCells.size());
        f.starStates.resize(f.interfaceCells.size());
        
    }
    // optional - print found boundary
    /*
    std::cout << std::endl; std::cout << "interface cells" << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        std::cout << interfaceCells[i][0] << " ";
    }
    std::cout << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        std::cout << interfaceCells[i][1] << " ";
    }
        */

}

void solver::printInterfaceArrays(std::vector<std::vector<double>> field, fluid& f, const std::string& filename) {
    // === Original arrays: ones, interface normals, and p ===
    std::vector<std::vector<double>> onesArray = field;
    std::vector<std::vector<double>> interfaceXArray = field;
    std::vector<std::vector<double>> interfaceYArray = field;
    std::vector<std::vector<double>> pArray = field;
    
    for (auto& row : onesArray)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceXArray)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceYArray)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : pArray)
        std::fill(row.begin(), row.end(), 0.0);
    
    // === Interface states arrays (8 arrays) ===
    // f.interfaceStates[i] is a std::array<std::array<double,4>,2>
    // We'll create two sets of 4 arrays.
    std::vector<std::vector<double>> interfaceState0_0 = field;
    std::vector<std::vector<double>> interfaceState0_1 = field;
    std::vector<std::vector<double>> interfaceState0_2 = field;
    std::vector<std::vector<double>> interfaceState0_3 = field;
    
    std::vector<std::vector<double>> interfaceState1_0 = field;
    std::vector<std::vector<double>> interfaceState1_1 = field;
    std::vector<std::vector<double>> interfaceState1_2 = field;
    std::vector<std::vector<double>> interfaceState1_3 = field;
    
    for (auto& row : interfaceState0_0)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState0_1)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState0_2)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState0_3)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState1_0)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState1_1)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState1_2)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : interfaceState1_3)
        std::fill(row.begin(), row.end(), 0.0);
    
    // === Star states arrays (4 arrays) ===
    // f.starStates[i] is a std::array<double,4>
    std::vector<std::vector<double>> starState0 = field;
    std::vector<std::vector<double>> starState1 = field;
    std::vector<std::vector<double>> starState2 = field;
    std::vector<std::vector<double>> starState3 = field;
    
    for (auto& row : starState0)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : starState1)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : starState2)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : starState3)
        std::fill(row.begin(), row.end(), 0.0);
    
    // === Fill arrays at interface cell locations ===
    // f.interfaceCells[i] gives the row and col for the i-th interface cell.
    for (int i = 0; i < f.interfaceCells.size(); ++i) {
        int r = f.interfaceCells[i][0];
        int c = f.interfaceCells[i][1];
        
        // Original arrays:
        onesArray[r][c] = 1.0;
        interfaceXArray[r][c] = f.interfaceNormals[i][0];
        interfaceYArray[r][c] = f.interfaceNormals[i][1];
        pArray[r][c] = f.u[r][c][3];
        
        // Interface states: for each cell, there are two sets (index 0 and 1) of 4 variables.
        interfaceState0_0[r][c] = f.interfaceStates[i][0][0];
        interfaceState0_1[r][c] = f.interfaceStates[i][0][1];
        interfaceState0_2[r][c] = f.interfaceStates[i][0][2];
        interfaceState0_3[r][c] = f.interfaceStates[i][0][3];
        
        interfaceState1_0[r][c] = f.interfaceStates[i][1][0];
        interfaceState1_1[r][c] = f.interfaceStates[i][1][1];
        interfaceState1_2[r][c] = f.interfaceStates[i][1][2];
        interfaceState1_3[r][c] = f.interfaceStates[i][1][3];
        
        // Star states: 4 variables per cell.
        starState0[r][c] = f.starStates[i][0];
        starState1[r][c] = f.starStates[i][1];
        starState2[r][c] = f.starStates[i][2];
        starState3[r][c] = f.starStates[i][3];
    }
    
    // === Write each array to a file ===
    auto writeArrayToFile = [](const std::vector<std::vector<double>>& arr, const std::string& fname) {
        std::ofstream outFile(fname);
        if (outFile.is_open()) {
            for (const auto& row : arr) {
                for (const auto& value : row)
                    outFile << std::fixed << std::setprecision(2) << value << " ";
                outFile << std::endl;
            }
            outFile.close();
        } else {
            std::cerr << "Unable to open file: " << fname << std::endl;
        }
    };
    
    // Original files:
    writeArrayToFile(onesArray, filename + "_ones.txt");
    writeArrayToFile(interfaceXArray, filename + "_x.txt");
    writeArrayToFile(interfaceYArray, filename + "_y.txt");
    writeArrayToFile(pArray, filename + "_p.txt");
    
    // Interface state files (8 files)
    writeArrayToFile(interfaceState0_0, filename + "_interfaceState0_0.txt");
    writeArrayToFile(interfaceState0_1, filename + "_interfaceState0_1.txt");
    writeArrayToFile(interfaceState0_2, filename + "_interfaceState0_2.txt");
    writeArrayToFile(interfaceState0_3, filename + "_interfaceState0_3.txt");
    
    writeArrayToFile(interfaceState1_0, filename + "_interfaceState1_0.txt");
    writeArrayToFile(interfaceState1_1, filename + "_interfaceState1_1.txt");
    writeArrayToFile(interfaceState1_2, filename + "_interfaceState1_2.txt");
    writeArrayToFile(interfaceState1_3, filename + "_interfaceState1_3.txt");
    
    // Star state files (4 files)
    writeArrayToFile(starState0, filename + "_starState0.txt");
    writeArrayToFile(starState1, filename + "_starState1.txt");
    writeArrayToFile(starState2, filename + "_starState2.txt");
    writeArrayToFile(starState3, filename + "_starState3.txt");
}


void solver::calcInterfaceNormals(fluid& f){
    //std::cout << std::endl;
    //printScalarField(phi);
    
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        int x = f.interfaceCells[i][1];
        int y = f.interfaceCells[i][0];
        double phi_y;
        double phi_x;
        if (x == nG){
            phi_x = (phi[y][x+1]-phi[y][x])/dx;
        } else if (x == nCellsX+nG-1){
            phi_x = (phi[y][x]-phi[y][x-1])/dx;
        } else {
            phi_x = (phi[y][x+1]-phi[y][x-1])/(2*dx);
        }
        if (y == nG){
            phi_y = (phi[y+1][x]-phi[y][x])/dy;
        } else if (y == nCellsY+nG-1){
            phi_y = (phi[y][x]-phi[y-1][x])/dy;
        } else {
            phi_y = (phi[y+1][x]-phi[y-1][x])/(2*dy);
        }

        std::array<double,2> grad = {phi_x/(sqrt(phi_x*phi_x+phi_y*phi_y)),phi_y/(sqrt(phi_x*phi_x+phi_y*phi_y))};
        f.interfaceNormals[i] = grad;
    }

    /*
    std::cout << std::endl; std::cout << "interface normals" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfaceNormals[i][0],2);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfaceNormals[i][1],2);
    }
    */
    

    
}

void solver::calcInterface(fluid& f){
    int fl1 = (&f == &fluid1) ? -1 : 1;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        double x_i = x0 + (f.interfaceCells[i][1]-nG+0.5)*dx; // -2 for ghost cells, +0.5 for cell centre
        double y_i = y0 + (f.interfaceCells[i][0]-nG+0.5)*dy;
        //std::cout << "boundary cell position " << x_i << " " << y_i << std::endl;

        double x = x_i - phi[f.interfaceCells[i][0]][f.interfaceCells[i][1]]*f.interfaceNormals[i][0];
        double y = y_i - phi[f.interfaceCells[i][0]][f.interfaceCells[i][1]]*f.interfaceNormals[i][1];
        //std::cout << "interface position " << x << " " << y << std::endl;

        // lPos always in fluid1, rPos always in fluid2
        std::array<double,2> lPos = {x-1.5*dx*f.interfaceNormals[i][0],y-1.5*dy*f.interfaceNormals[i][1]};
        std::array<double,2> rPos = {x+1.5*dx*f.interfaceNormals[i][0],y+1.5*dy*f.interfaceNormals[i][1]};
        //std::cout << "left interface position " << lPos[0] << " " << lPos[1] << std::endl;
        //std::cout << "right interface position " << rPos[0] << " " << rPos[1] << std::endl;
        f.interfacePositions[i] = {lPos,rPos};
    }

    // optional - print positions
    /*
    std::cout << std::endl; std::cout << "left interface positions" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfacePositions[i][0][0],2);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfacePositions[i][0][1],2);
    }
    std::cout << std::endl; std::cout << "right interface positions" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfacePositions[i][1][0],2);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfacePositions[i][1][1],2);
    }
    */
    

}

void solver::interpInterfaceStates(fluid& f){ // Calculates interpolated states (primitive) at a set distance normal to the interface
    //int fl1 = (&f == &fluid1) ? -1 : 1;
    for (size_t i = 0; i < f.interfacePositions.size(); ++i){
        //std::cout << "which fluid: " << fl1 << std::endl;
        //std::cout << "boundary cell" << f.interfaceCells[i][0] << " " << f.interfaceCells[i][1] << std::endl;
        //std::cout << "phi value" << phi[f.interfaceCells[i][0]][f.interfaceCells[i][1]] << std::endl;
        //std::cout << "interface pos"; printState(f.interfacePositions[i][0]); printState(f.interfacePositions[i][1]);
        std::array<double,4> lState = bilinear4(f.interfacePositions[i][0][0],f.interfacePositions[i][0][1],fluid1, f.interfaceNormals[i],i);
        std::array<double,4> rState = bilinear4(f.interfacePositions[i][1][0],f.interfacePositions[i][1][1],fluid2, f.interfaceNormals[i],i);
        //print_state(lState); print_state(rState);
        
        f.interfaceStates[i] = {eos[0]->consvToPrim(lState),eos[1]->consvToPrim(rState)};
    }


    // optional - print states
    /*
    std::cout << std::endl; std::cout << "L/R UY" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfaceStates[i][0][UY],3);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfaceStates[i][1][UY],3);
    }
    std::cout << std::endl; std::cout << "L/R rho" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfaceStates[i][0][RHO],3);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.interfaceStates[i][1][RHO],3);
    }
    std::cout << std::endl;
    */
    
}

std::array<double,4> solver::bilinear(double x, double y, fluid& fluid, std::array<double,2> norm){
    int rows = nCellsY;
    int cols = nCellsX;

    // Compute indices of grid cell containing the point
    int j = ((x - x0) / dx) + nG - 0.5;
    int i = ((y - y0) / dy) + nG - 0.5;

    //print_state(fluid.u[i][j]);

    /*
    
    if (phi[i][j] > 0){
        if (norm[0] >= 0 && norm[1] >= 0){
            jR = j+1; iR = i+1;
        } else if (norm[0] >= 0){
            jR = j+1; iR = i-1;
        } else if (norm[0] < 0){
            jR = j-1; iR = i-1;
        } else {
            jR = j-1; iR = i+1;
        }
    } else {
        if (norm[0] >= 0 && norm[1] >= 0){
            jR = j-1; iR = i-1;
        } else if (norm[0] >= 0){
            jR = j-1; iR = i+1;
        } else if (norm[0] < 0){
            jR = j+1; iR = i+1;
        } else {
            jR = j+1; iR = i-1;
        }
    }
        */

    int jR,iR;
    if (phi[i][j] > 0) {
        jR = j + (norm[0] >= 0 ? 1 : -1);
        iR = i + (norm[1] >= 0 ? 1 : -1);
    } else {
        jR = j - (norm[0] >= 0 ? 1 : -1);
        iR = i - (norm[1] >= 0 ? 1 : -1);
    }

    if (i == 0 || i == (2*nG+nCellsY) || j==0 || j == (2*nG+nCellsX)){
        return fluid.u[i][j];
    }

    /*
    if (std::signbit(phi[i][j]) != std::signbit(phi[iR][j]) || std::signbit(phi[i][jR]) != std::signbit(phi[i][j])
        || std::signbit(phi[iR][jR]) != std::signbit(phi[iR][j]) || std::signbit(phi[iR][jR]) != std::signbit(phi[i][jR])){
        return fluid.u[i][j];
    } else {
        return fluid.u[i][j];
    }
        */

    if ((iR >= 0 && iR < rows+2*nG) == 0){
        std::cout << "y x ij = " << y << " " << x << " " << i << " " << j << std::endl;
        throw std::runtime_error("i interpolation out of bounds");
    }
    if ((jR >= 0 && jR < cols+2*nG) == 0){
        std::cout << "y x ij = " << y << " " << x << " " << i << " " << j << std::endl;
        throw std::runtime_error("j interpolation out of bounds");
    }


    //std::cout << "i = " << i << " j = " << j << "iR = " << iR << "jR = " << jR << std::endl;

    //if (i < 0 || i > rows - 1 || j < 0 || j > cols - 1) {
        //std::cout << i << " " << j << std::endl;
        //return {NAN,NAN,NAN,NAN};
        //throw std::runtime_error("interpolation out of bounds");
    //};

    // Physical coordinates of the four surrounding points
    /*
    double xL = x0 + (j - nG + 0.5) * dx;
    double xR = x0 + (j + 1 - nG + 0.5) * dx;
    double yL = y0 + (i - nG + 0.5) * dy;
    double yR = y0 + (i + 1 - nG + 0.5) * dy;
    */
    double xL = x0 + (j - nG + 0.5) * dx;
    double xR = x0 + (jR - nG + 0.5) * dx;
    double yL = y0 + (i - nG + 0.5) * dy;
    double yR = y0 + (iR - nG + 0.5) * dy;
    
    std::array<double,4> f00 = fluid.u[i][j];
    std::array<double,4> f10 = fluid.u[iR][j];
    std::array<double,4> f01 = fluid.u[i][jR];
    std::array<double,4> f11 = fluid.u[iR][jR];

    assert(f00[RHO] > 0 && f01[RHO] > 0 && f10[RHO] > 0 && f11[RHO] > 0);
    assert(f00[PRES] > 0 && f01[PRES] > 0 && f10[PRES] > 0 && f11[PRES] > 0);

    
    //std::cout << "interpolation points" << std::endl;
    //std::cout << i << " " << j << " " << f00[RHO] << " " << f00[UX]<< " "  << f00[UY]<< " "  << f00[PRES]<< " "  << std::endl;
    //std::cout << i << " " << jR << " " << f01[RHO] << " " << f01[UX]<< " "  << f01[UY]<< " "  << f01[PRES]<< " "  << std::endl;
    //std::cout << iR << " " << j << " " << f10[RHO] << " " << f10[UX]<< " "  << f10[UY]<< " "  << f10[PRES]<< " "  << std::endl;
    //std::cout << iR << " " << jR << " " << f11[RHO] << " " << f11[UX]<< " "  << f11[UY]<< " "  << f11[PRES]<< " "  << std::endl;
    
    
    

    //std::cout << "neighbour values" << f00[0] << " " << f10[0] << " " << f01[0] << " " << f11[0] << std::endl;

    // Interpolation along x at y0 and y1
    std::array<double,4> f_x_y0 = (xR - x) / (xR - xL) * f00 + (x - xL) / (xR - xL) * f01;
    std::array<double,4> f_x_y1 = (xR - x) / (xR - xL) * f10 + (x - xL) / (xR - xL) * f11;
    
    // Interp along y
    std::array<double,4> res = (yR - y) / (yR - yL) * f_x_y0 + (y - yL) / (yR - yL) * f_x_y1;
    //std::cout << "result:" << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << std::endl;
    if (res[RHO] <= 0){
        //std::cout << i << " " << iR << " " << j << " " << jR << std::endl;
        //throw std::runtime_error("Negative density extrapolated");
        double minRho = 1e10;
        for (double i : {f00[RHO],f01[RHO],f10[RHO],f11[RHO]}){
            minRho = std::min(minRho,i);
        }
        res[RHO] = minRho;
    }
    
    if (res[PRES] <= 0){
        //std::cout << i << " " << iR << " " << j << " " << jR << std::endl;
        //throw std::runtime_error("Negative density extrapolated");
        double minPres = 1e10;
        for (double i : {f00[PRES],f01[PRES],f10[PRES],f11[PRES]}){
            minPres = std::min(minPres,i);
        }
        res[PRES] = minPres;
    }
    
    return res;
}

std::array<double,4> solver::bilinear4(double x, double y, fluid& fluid, std::array<double,2> norm, int index) {
    int rows = nCellsY;
    int cols = nCellsX;

    // Expected sign for the fluid: -1 for fluid1 and 1 for fluid2.
    int fl1 = (&fluid == &fluid1) ? -1 : 1;

    // Compute indices of the grid cell containing the point (centered by subtracting 0.5)
    double jdub = ((x - x0) / dx) + nG - 0.5;
    double idub = ((y - y0) / dy) + nG - 0.5;

    int j = static_cast<int>(std::round(jdub));
    int i = static_cast<int>(std::round(idub));

    // --- selection of interpolation points in the normal direction ---
    int jR, iR;
    if (fl1 == 1) {
        // For fluid2 (real when phi > 0): use positive offset if normal component is >= 0
        jR = j + (norm[0] >= 0 ? 1 : -1);
        iR = i + (norm[1] >= 0 ? 1 : -1);
    } else {
        // For fluid1 (real when phi <= 0): reverse the offset compared to fluid2
        jR = j - (norm[0] >= 0 ? 1 : -1);
        iR = i - (norm[1] >= 0 ? 1 : -1);
    }

    // Check boundaries for the base cell; if near boundaries, simply return its value if valid.
    if (i <= 0 || i >= (2*nG + nCellsY - 1) ||
        j <= 0 || j >= (2*nG + nCellsX - 1)) {
        // Treat the base cell as any other cell.
        if (std::signbit(phi[i][j]) != std::signbit(fl1))
            throw std::runtime_error("Base cell outside real fluid (boundary case)");
        return fluid.u[i][j];
    }
    // Check boundaries for the neighbor cell.
    if (iR < 0 || iR >= (2*nG + nCellsY) ||
        jR < 0 || jR >= (2*nG + nCellsX)) {
        std::cout << "y x ij = " << y << " " << x << " " << i << " " << j << std::endl;
        throw std::runtime_error("Interpolation neighbor out of bounds");
    }

    // Physical coordinates of cell centers.
    double xL      = x0 + (j - nG + 0.5) * dx;
    double xR_coord = x0 + (jR - nG + 0.5) * dx;
    double yL      = y0 + (i - nG + 0.5) * dy;
    double yR_coord = y0 + (iR - nG + 0.5) * dy;

    // Retrieve the fluid state at the four grid points.
    std::array<double,4> f00 = fluid.u[i][j];    // Base cell.
    std::array<double,4> f10 = fluid.u[iR][j];     // Cell vertically (or diagonally) offset.
    std::array<double,4> f01 = fluid.u[i][jR];     // Cell horizontally offset.
    std::array<double,4> f11 = fluid.u[iR][jR];     // Diagonal neighbor.

    // Sanity-check: density and pressure must be positive in all cells.
    assert(f00[RHO] > 0 && f01[RHO] > 0 && f10[RHO] > 0 && f11[RHO] > 0);
    assert(f00[PRES] > 0 && f01[PRES] > 0 && f10[PRES] > 0 && f11[PRES] > 0);

    // --- (b) Mark valid cells: each cell is valid only if its sign matches fl1.
    bool valid00 = (std::signbit(phi[i][j])  == std::signbit(fl1));
    bool valid10 = (std::signbit(phi[iR][j]) == std::signbit(fl1));
    bool valid01 = (std::signbit(phi[i][jR]) == std::signbit(fl1));
    bool valid11 = (std::signbit(phi[iR][jR])== std::signbit(fl1));

    // Compute standard bilinear interpolation weights.
    double denom_x = (xR_coord - xL);
    double denom_y = (yR_coord - yL);
    if (denom_x == 0 || denom_y == 0) {
        if (!valid00)
            throw std::runtime_error("Degenerate cell spacing and invalid base cell");
        return fluid.u[i][j];
    }
    double w00 = ((xR_coord - x) / denom_x) * ((yR_coord - y) / denom_y);
    double w01 = ((x - xL) / denom_x) * ((yR_coord - y) / denom_y);
    double w10 = ((xR_coord - x) / denom_x) * ((y - yL) / denom_y);
    double w11 = ((x - xL) / denom_x) * ((y - yL) / denom_y);

    // Zero out the weights of cells that are not in the correct fluid.
    if (!valid00) w00 = 0;
    if (!valid01) w01 = 0;
    if (!valid10) w10 = 0;
    if (!valid11) w11 = 0;

    double wsum = w00 + w01 + w10 + w11;
    if (wsum <= 0) {
        //std::cout << fl1 << " " << " " << idub << " " << jdub << " " << i << " " << j << " " << iR << " " << jR << " " << norm[0] << " " << norm[1] << std::endl;
        //throw std::runtime_error("No valid neighbors in the correct fluid for bilinear interpolation");
        return fluid.u[fluid.interfaceCells[index][0]][fluid.interfaceCells[index][1]];
    }
    // Normalize the weights and compute the weighted sum.
    w00 /= wsum;
    w01 /= wsum;
    w10 /= wsum;
    w11 /= wsum;
    std::array<double,4> res = w00 * f00 + w01 * f01 + w10 * f10 + w11 * f11;

    // Safeguards: ensure density and pressure remain positive.
    if (res[RHO] <= 0) {
        double minRho = std::min({f00[RHO], f01[RHO], f10[RHO], f11[RHO]});
        res[RHO] = minRho;
    }
    if (res[PRES] <= 0) {
        double minPres = std::min({f00[PRES], f01[PRES], f10[PRES], f11[PRES]});
        res[PRES] = minPres;
    }
    
    return res;
}

std::array<double,4> solver::bilinear2(double x, double y, fluid& fluid, std::array<double,2> norm) {
    int rows = nCellsY;
    int cols = nCellsX;

    int fl1 = (&fluid == &fluid1) ? -1 : 1;

    // Compute indices of grid cell containing the point (centered by subtracting 0.5)
    int j = ((x - x0) / dx) + nG - 0.5;
    int i = ((y - y0) / dy) + nG - 0.5;

    if (std::signbit(phi[i][j]) != std::signbit(fl1)){
        throw std::runtime_error("point outside real fluid");
    }


    // --- (a) Symmetric neighbor selection ---
    int jR, iR;
    if (phi[i][j] > 0) {
        // For fluid2 (real when phi > 0): use positive offset if normal component is >= 0
        jR = j + (norm[0] >= 0 ? 1 : -1);
        iR = i + (norm[1] >= 0 ? 1 : -1);
    } else {
        // For fluid1 (real when phi <= 0): reverse the offset compared to fluid2
        jR = j - (norm[0] >= 0 ? 1 : -1);
        iR = i - (norm[1] >= 0 ? 1 : -1);
    }

    // Check boundaries for the base cell; if near boundaries, simply return its value.
    if (i <= 0 || i >= (2*nG + nCellsY - 1) ||
        j <= 0 || j >= (2*nG + nCellsX - 1)) {
        return fluid.u[i][j];
    }
    // Check boundaries for the neighbor cell.
    if (iR < 0 || iR >= (2*nG + nCellsY) ||
        jR < 0 || jR >= (2*nG + nCellsX)) {
        std::cout << "y x ij = " << y << " " << x << " " << i << " " << j << std::endl;
        throw std::runtime_error("Interpolation neighbor out of bounds");
    }

    // Physical coordinates of cell centers.
    double xL = x0 + (j - nG + 0.5) * dx;
    double xR_coord = x0 + (jR - nG + 0.5) * dx;
    double yL = y0 + (i - nG + 0.5) * dy;
    double yR_coord = y0 + (iR - nG + 0.5) * dy;

    // Retrieve the fluid state at the four grid points.
    std::array<double,4> f00 = fluid.u[i][j];   // Base cell.
    std::array<double,4> f10 = fluid.u[iR][j];    // Cell vertically (or diagonally) offset.
    std::array<double,4> f01 = fluid.u[i][jR];    // Cell horizontally offset.
    std::array<double,4> f11 = fluid.u[iR][jR];     // Diagonal neighbor.

    // Sanity-check: density and pressure must be positive in all cells.
    assert(f00[RHO] > 0 && f01[RHO] > 0 && f10[RHO] > 0 && f11[RHO] > 0);
    assert(f00[PRES] > 0 && f01[PRES] > 0 && f10[PRES] > 0 && f11[PRES] > 0);

    // --- (b) Check that all cells are in the real fluid ---
    // Determine the "real" fluid based on the base cell.
    bool baseIsFluid2 = (phi[i][j] > 0);  // true means fluid2 is real; false means fluid1 is real.
    // For each cell, check that its phi agrees with the base cell.
    bool valid00 = true;  // Base cell is assumed valid.
    bool valid10 = ((phi[iR][j] > 0) == baseIsFluid2);
    bool valid01 = ((phi[i][jR] > 0) == baseIsFluid2);
    bool valid11 = ((phi[iR][jR] > 0) == baseIsFluid2);

    // Compute standard bilinear interpolation weights.
    double denom_x = (xR_coord - xL);
    double denom_y = (yR_coord - yL);
    if (denom_x == 0 || denom_y == 0) {
        // Degenerate cell spacing; return the base cell's value.
        return fluid.u[i][j];
    }
    double w00 = ((xR_coord - x) / denom_x) * ((yR_coord - y) / denom_y);
    double w01 = ((x - xL) / denom_x) * ((yR_coord - y) / denom_y);
    double w10 = ((xR_coord - x) / denom_x) * ((y - yL) / denom_y);
    double w11 = ((x - xL) / denom_x) * ((y - yL) / denom_y);

    // Set weights to zero for cells not in the real fluid.
    if (!valid00) w00 = 0;
    if (!valid01) w01 = 0;
    if (!valid10) w10 = 0;
    if (!valid11) w11 = 0;

    double wsum = w00 + w01 + w10 + w11;
    std::array<double,4> res;
    if (wsum > 0) {
        // Normalize the weights and compute the weighted sum.
        w00 /= wsum;
        w01 /= wsum;
        w10 /= wsum;
        w11 /= wsum;
        res = w00 * f00 + w01 * f01 + w10 * f10 + w11 * f11;
    } else {
        // If no valid neighbors remain, fall back to the base cell.
        res = fluid.u[i][j];
    }

    // Safeguards: ensure density and pressure remain positive.
    if (res[RHO] <= 0) {
        double minRho = std::min({f00[RHO], f01[RHO], f10[RHO], f11[RHO]});
        res[RHO] = minRho;
    }
    if (res[PRES] <= 0) {
        double minPres = std::min({f00[PRES], f01[PRES], f10[PRES], f11[PRES]});
        res[PRES] = minPres;
    }
    
    return res;
}

std::array<double,4> solver::bilinear3(double x, double y, fluid& fluid, std::array<double,2> norm) {
    int rows = nCellsY;
    int cols = nCellsX;

    // Compute indices of grid cell containing the point (centered by subtracting 0.5)
    int j = ((x - x0) / dx) + nG - 0.5;
    int i = ((y - y0) / dy) + nG - 0.5;
    
    return neighbourAvgInterp(fluid,i,j);
}

void solver::resolveVelocities(fluid& f){
    /*
    resolves normal and tangential components of velocities of the interpolated states
    */
    for (size_t i=0; i<f.interfaceStates.size(); ++i){
        double vNormL = f.interfaceStates[i][0][UX]*f.interfaceNormals[i][0]+f.interfaceStates[i][0][UY]*f.interfaceNormals[i][1];
        double vNormR = f.interfaceStates[i][1][UX]*f.interfaceNormals[i][0]+f.interfaceStates[i][1][UY]*f.interfaceNormals[i][1];

        //double vTanL = sqrt((f.interfaceStates[i][0][UX]*f.interfaceStates[i][0][UX]+f.interfaceStates[i][0][UY]*f.interfaceStates[i][0][UY])-vNormL*vNormL);
        //double vTanR = sqrt((f.interfaceStates[i][1][UX]*f.interfaceStates[i][1][UX]+f.interfaceStates[i][1][UY]*f.interfaceStates[i][1][UY])-vNormR*vNormR);
        double vTanL = f.interfaceStates[i][0][UX]*f.interfaceNormals[i][1]-f.interfaceStates[i][0][UY]*f.interfaceNormals[i][0];
        double vTanR = f.interfaceStates[i][1][UX]*f.interfaceNormals[i][1]-f.interfaceStates[i][1][UY]*f.interfaceNormals[i][0];
        
        std::array<double,2> vL = {vNormL,vTanL};
        std::array<double,2> vR = {vNormR,vTanR};
        f.resolvedVelocities[i] = {vL,vR};
    }

    /*
    std::cout << std::endl; std::cout << "vL" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.resolvedVelocities[i][0][0],3);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.resolvedVelocities[i][0][1],3);
    }
    std::cout << std::endl; std::cout << "vR" << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.resolvedVelocities[i][1][0],3);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        printPrec(f.resolvedVelocities[i][1][1],3);
    }
    */
    
}

std::array<double,4> solver::exactRiemann(std::array<double,4> left,std::array<double,4> right, EOS* eos, bool direction){
    double d = (direction == XDIR) ? dx : dy;
    riemann solution(eos->get_gamma(),eos->get_gamma(),eos->consvToPrim(left),eos->consvToPrim(right),direction,-0.5*d,0.5*d,0,dt,2,0,0);
    std::array<double,4> result = solution.exctRiemann();
    if (result[RHO] <= 0){
        std::array<double,4> l1 = eos->consvToPrim(left);
        std::array<double,4> r1 = eos->consvToPrim(right);
        printState(l1); printState(r1);
        throw std::runtime_error("flux riemann solution rho < 0");
    }
    //if (sqrt(result[UX]*result[UX]) < 1e-4){result[UX] = 0.0;}
    //if (sqrt(result[UY]*result[UY]) < 1e-4){result[UY] = 0.0;}
    return flux(result,eos,direction);
}

void solver::interfaceRiem(fluid& f){
    bool left = (&f == &fluid1) ? false : true;

    for (size_t i=0; i<f.interfaceStates.size(); ++i){
        std::array<double,4> lState = {f.interfaceStates[i][0][0],f.resolvedVelocities[i][0][0],f.resolvedVelocities[i][0][1],f.interfaceStates[i][0][3]};  
        std::array<double,4> rState = {f.interfaceStates[i][1][0],f.resolvedVelocities[i][1][0],f.resolvedVelocities[i][1][1],f.interfaceStates[i][1][3]};
        
        //std::cout << "riem " << f.interfaceCells[i][0] << " " << f.interfaceCells[i][1] << " bool left = " << left << " normal:" << f.interfaceNormals[i][0] <<" " << f.interfaceNormals[i][1]<< std::endl;
        
        //print_state(f.interfaceStates[i][0]); print_state(f.interfaceStates[i][1]);
        //print_state(lState); print_state(rState);

        riemann solution(eos[0]->get_gamma(),eos[1]->get_gamma(),lState,rState,1,-1.5*(dx+dy)/2.0,1.5*(dx+dy)/2.0,0,dt,2,0,0); // need p_inf setup
        std::array<double,4> result = solution.interfaceRiemann(left);
        //std::array<std::array<double,4>,2> result = HLLCStates(lState,rState,eos[0]);
        //print_state(result[0]); print_state(result[1]);
        //if (result[RHO] <= 0){
        //    print_state(lState); print_state(rState);
        //    throw std::runtime_error("interface riemann solution rho < 0");
        //}
        f.riemInterfaceStates[i] = result;

        //f.riemInterfaceStates[i] = left ? result[0] : result[1];
        
    }
}

void solver::calcStarStates(fluid& f){
    for (size_t i=0; i<f.interfaceStates.size(); ++i){
            
            double uStarx = f.riemInterfaceStates[i][UX]*f.interfaceNormals[i][0]+f.riemInterfaceStates[i][UY]*f.interfaceNormals[i][1];
            double uStary = f.riemInterfaceStates[i][UX]*f.interfaceNormals[i][1]-f.riemInterfaceStates[i][UY]*f.interfaceNormals[i][0];

            std::array<double,4> uStar = {f.riemInterfaceStates[i][0],uStarx,uStary,f.riemInterfaceStates[i][3]};
            f.starStates[i] = uStar;

            //std::cout << "star states left: " << uStarL[0] << " " << uStarL[1] << " " << uStarL[2] << " " << uStarL[3] << std::endl;
    }
}

void solver::calcPhiGrad(){
    phiBC();

    #pragma omp parallel for collapse (2)
    for (int y = nG; y < nCellsY+nG; ++y){
        for (int x = nG; x < nCellsX+nG; ++x){
            double phi_y;
            double phi_x;
            if (x == nG){
                phi_x = (phi[y][x+1]-phi[y][x])/dx;
            } else if (x == nCellsX+nG-1){
                phi_x = (phi[y][x]-phi[y][x-1])/dx;
            } else {
                phi_x = (phi[y][x+1]-phi[y][x-1])/(2*dx);
            }
            if (y == nG){
                phi_y = (phi[y+1][x]-phi[y][x])/dy;
            } else if (y == nCellsY+nG-1){
                phi_y = (phi[y][x]-phi[y-1][x])/dy;
            } else {
                phi_y = (phi[y+1][x]-phi[y-1][x])/(2*dy);
            }


            std::array<double,2> norm = {0,0};
            if(((phi_x == 0) && (phi_y == 0)) == 1){
                norm = {0,0};
            } else {
                norm[0] = phi_x/(sqrt(phi_x*phi_x+phi_y*phi_y)); norm[1] = phi_y/(sqrt(phi_x*phi_x+phi_y*phi_y));
            }
            //std::cout << norm[0] << " " << norm[1] << std::endl;
            phiNormals[y][x] = norm;
            phiGrads[y][x] = sqrt(phi_x*phi_x+phi_y*phi_y);
        }
    }
    //vectorTransBC(phiNormals);
    //scalarTransBC(phiGrads);
    //printScalarField(phiGrads);
}

void solver::setInterface(){ // ind = 0 for fluid1, 1 for fluid2

    // setting uExtrap to the ghost fluids
    //std::cout << "no of star states = " << starStates.size() << std::endl;
    for (size_t i = 0; i<fluid1.starStates.size(); ++i){
        fluid1.u[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]] = eos[0]->primToConsv(fluid1.starStates[i]); // interface states set to u*R
    }
    for (size_t i = 0; i<fluid2.starStates.size(); ++i){
        fluid2.u[fluid2.interfaceCells[i][0]][fluid2.interfaceCells[i][1]] = eos[1]->primToConsv(fluid2.starStates[i]); // interface states set to u*R
    }
    setBCs(fluid1);
    setBCs(fluid2);
    /*
    for (size_t i = 0; i<fluid1.starStates.size(); ++i){
        fluid2.u[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]] = eos[1]->primToConsv(fluid2.starStates[i][1]); // interface states set to u*R
        if (phi[fluid1.interfaceCells[i][0]+1][fluid1.interfaceCells[i][1]] > 0){ // surrounding states set to u*L
            fluid1.u[fluid1.interfaceCells[i][0]+1][fluid1.interfaceCells[i][1]] = eos[0]->primToConsv(fluid1.starStates[i][0]);
        }
        if (phi[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]+1] > 0){
            fluid1.u[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]+1] = eos[0]->primToConsv(fluid1.starStates[i][0]);
        }
        if (phi[fluid1.interfaceCells[i][0]-1][fluid1.interfaceCells[i][1]] > 0){
            fluid1.u[fluid1.interfaceCells[i][0]-1][fluid1.interfaceCells[i][1]] = eos[0]->primToConsv(fluid1.starStates[i][0]);
        }
        if (phi[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]-1] > 0){
            fluid1.u[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]-1] = eos[0]->primToConsv(fluid1.starStates[i][0]);
        }
    }
    */

}

int sign(double value) {
    return (value <= 0) ? -1 : 1;
}

std::array<double,4> solver::solveQ(bool a, bool b, std::array<double,4> Qx,std::array<double,4> Qy, double alpha){
    std::array<double,4> Q;
    int signA = 2 * a - 1;
    int signB = 2 * b - 1;
    //std::cout << "A,B: " << signA << " " << signB << std::endl;
    //print_state(Qx); print_state(Qy);
    std::array<double,4> numer = ((signA*Qx)/alpha) + signB*Qy;
    //print_state(numer);
    double denom = signA/alpha + signB;
    Q = numer/denom;
    //print_state(Q);
    return Q;
}

void solver::solveSweepPoint(int i, int j){
    std::array<double,4> Qx, Qy;
    if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
        && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
            
            double alpha = (dx/dy)*phiNormals[i][j][1]/(phiNormals[i][j][0]);
            if (std::isnan(alpha) == 1){
                alpha = dx/dy;
            } else if (std::isinf(alpha) == 1){
                alpha = 1e15;
            } else if (alpha == 0){
                alpha = 1e-15;
            }
            //std::cout << alpha << std::endl;

            if (phi[i][j] <= 0){
                Qx = (phiNormals[i][j][0] > 0) ? fluid2.u[i][j+1] : fluid2.u[i][j-1];
                Qy = (phiNormals[i][j][1] > 0) ? fluid2.u[i+1][j] : fluid2.u[i-1][j];
                fluid2.u[i][j] = solveQ((phiNormals[i][j][0] > 0),(phiNormals[i][j][1] > 0),Qx,Qy,alpha);
            } else {
                Qx = (phiNormals[i][j][0] <= 0) ? fluid1.u[i][j+1] : fluid1.u[i][j-1];
                Qy = (phiNormals[i][j][1] <= 0) ? fluid1.u[i+1][j] : fluid1.u[i-1][j];
                fluid1.u[i][j] = solveQ((phiNormals[i][j][0] <= 0),(phiNormals[i][j][1] <= 0),Qx,Qy,alpha);
            }
        }
}

void solver::GFFastSweeping2(){
    double maxIter = (splitFlip == 0) ? 20 : 1;
    calcPhiGrad();
    //std::cout << "calced phi grad" << std::endl;
    setInterface();
    //std::cout << "set interface" << std::endl;

    for (int iter = 0; iter < maxIter; ++iter){
        // Sweep 1: i ascending, j ascending
        for (int i = nG; i < fluid1.u.size()-nG; ++i){
            for (int j = nG; j < fluid1.u[0].size()-nG; ++j){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 2: j ascending, i ascending
        for (int j = nG; j < fluid1.u[0].size()-nG; ++j){
            for (int i = nG; i < fluid1.u.size()-nG; ++i){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 3: i ascending, j descending
        for (int i = nG; i < fluid1.u.size()-nG; ++i){
            for (int j = fluid1.u[0].size()-1-nG; j >= nG; --j){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 4: j ascending, i descending
        for (int j = nG; j < fluid1.u[0].size()-nG; ++j){
            for (int i = fluid1.u.size()-1-nG; i >= nG; --i){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 5: i descending, j ascending
        for (int i = fluid1.u.size()-1-nG; i >= nG; --i) {
            for (int j = nG; j < fluid1.u[0].size()-nG; ++j) {
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 6: i descending, j descending
        for (int i = fluid1.u.size()-1-nG; i >= nG; --i) {
            for (int j = fluid1.u[0].size()-1-nG; j >= nG; --j) {
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 7: j descending, i ascending
        for (int j = fluid1.u[0].size()-1-nG; j >= nG; --j) {
            for (int i = nG; i < fluid1.u.size()-nG; ++i) {
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);

        // Sweep 8: j descending, i descending
        for (int j = fluid1.u[0].size()-1-nG; j >= nG; --j) {
            for (int i = fluid1.u.size()-1-nG; i >= nG; --i) {
                solveSweepPoint(i,j);
            }
        }
        // Final boundary condition update after all sweeps
        setBCs(fluid1);
        transmissiveBC(fluid2);
    }
}

void solver::GFFastSweeping(){
    double maxIter = (splitFlip == 0) ? 20 : 1;
    calcPhiGrad();
    //std::cout << "calced phi grad" << std::endl;
    setInterface();
    //std::cout << "set interface" << std::endl;
    //printInterfaceArrays(phi, fluid1, "interface1");
    //printInterfaceArrays(phi, fluid2, "interface2");
    //std::cout << "phi normals calculated" << std::endl;
    for (int iter = 0; iter < maxIter; ++iter){
        for (int i = nG; i < fluid1.u.size()-nG; ++i){
            for (int j = nG; j < fluid1.u[0].size()-nG; ++j){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);
        //transmissiveBC(fluid2);
        for (int j = nG; j < fluid1.u[0].size()-nG; ++j){
            for (int i = nG; i < fluid1.u.size()-nG; ++i){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);
        for (int i = nG; i < fluid1.u.size()-nG; ++i){
            for (int j = fluid1.u[0].size()-1-nG; j >= nG; --j){
                solveSweepPoint(i,j);
            }
        }
        setBCs(fluid1); setBCs(fluid2);
        for (int j = nG; j < fluid1.u[0].size()-nG; ++j){
            for (int i = fluid1.u.size()-1-nG; i >= nG; --i){
                solveSweepPoint(i,j);
            }
        }

        setBCs(fluid1);
        transmissiveBC(fluid2);
        
    }
}

void solver::fastSweepExtrapolation(std::vector<std::vector<std::array<double, 4>>>& grid,
                            std::vector<std::array<int, 2>>& interfaceCells,std::vector<std::vector<double>>& phi,
                            fluid& f, int numSweeps = 100) {

    bool fl1 = (&f == &fluid1) ? true : false;
    // Determine grid dimensions.
    int NX = grid.size();
    if (NX == 0) return;
    int NY = grid[0].size();

    setInterface();

    // Create a 2D boolean mask to mark known (interface) cells.
    std::vector<std::vector<bool>> known(NX, std::vector<bool>(NY, false));
    for (const auto &cell : interfaceCells) {
        int i = static_cast<int>(cell[0]);
        int j = static_cast<int>(cell[1]);
        if (i >= 0 && i < NX && j >= 0 && j < NY) {
            known[i][j] = true;
        }
    }

    // Lambda function to update a cell's value from its neighbors,
    // but only if the cell is not a known (interface) cell.
    auto updateCell = [&](int i, int j) {
        if (known[i][j] || (phi[i][j]<=0)==fl1) return;  // Skip known cells.

        std::array<double, 4> sum = {0.0, 0.0, 0.0, 0.0};
        int count = 0;

        // Up neighbor
        if (i > 0) {
            for (int k = 0; k < 4; k++) {
                sum[k] += grid[i - 1][j][k];
            }
            count++;
        }
        // Left neighbor
        if (j > 0) {
            for (int k = 0; k < 4; k++) {
                sum[k] += grid[i][j - 1][k];
            }
            count++;
        }
        // Down neighbor
        if (i < NX - 1) {
            for (int k = 0; k < 4; k++) {
                sum[k] += grid[i + 1][j][k];
            }
            count++;
        }
        // Right neighbor
        if (j < NY - 1) {
            for (int k = 0; k < 4; k++) {
                sum[k] += grid[i][j + 1][k];
            }
            count++;
        }

        // Update the cell if any neighbors exist.
        if (count > 0) {
            for (int k = 0; k < 4; k++) {
                grid[i][j][k] = sum[k] / count;
            }
        }
    };

    // Perform multiple sweeps in four different orders to propagate known values.
    for (int sweep = 0; sweep < numSweeps; sweep++) {
        // Order 1: i increasing, j increasing.
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                updateCell(i, j);
            }
        }
        // Order 2: i decreasing, j decreasing.
        for (int i = NX - 1; i >= 0; i--) {
            for (int j = NY - 1; j >= 0; j--) {
                updateCell(i, j);
            }
        }
        // Order 3: i increasing, j decreasing.
        for (int i = 0; i < NX; i++) {
            for (int j = NY - 1; j >= 0; j--) {
                updateCell(i, j);
            }
        }
        // Order 4: i decreasing, j increasing.
        for (int i = NX - 1; i >= 0; i--) {
            for (int j = 0; j < NY; j++) {
                updateCell(i, j);
            }
        }
    }
}

void solver::SLIC(bool direction,fluid& f, EOS* eos){
    if (PRIM == 1){
        //std::cout << "starting primitive SLIC" << std::endl;
        for (size_t i = 0; i < f.uPrim.size(); ++i){
            for (size_t j = 0; j < f.uPrim[0].size(); ++j){
                f.uPrim[i][j] = eos->consvToPrim(f.u[i][j]);
            }
        }

    } else {
        //std::cout << "starting conservative SLIC" << std::endl;
    }

    //std::cout << "uPrim calculated, calc halfslopes" << std::endl;
    calcHalfSlopes(PRIM,direction,f);

    //std::cout << "calcr" << std::endl;
    calcr(direction,f);

    //std::cout << "calcUBars" << std::endl;
    calcUBars(PRIM,direction,f,eos);

    //std::cout << "updUBars" << std::endl;
    updUBars(PRIM,direction,f,eos);

    if (direction == XDIR){
        for (size_t i = 0; i < f.fluxesX.size(); ++i) {
            for (size_t j = 0; j < f.fluxesX[0].size(); ++j){
                f.fluxesX[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(f.uBarRX[i][j],f.uBarLX[i][j+1],eos,XDIR)),eos,XDIR)+LF(f.uBarRX[i][j],f.uBarLX[i][j+1],eos,XDIR)); // fluxesX stored in conserved variable form
            }
        }
    } else {
        for (size_t j = 0; j < f.fluxesY[0].size(); ++j) {
            for (size_t i = 0; i < f.fluxesY.size(); ++i){
                f.fluxesY[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(f.uBarRY[i][j],f.uBarLY[i+1][j],eos,YDIR)),eos,YDIR)+LF(f.uBarRY[i][j],f.uBarLY[i+1][j],eos,YDIR)); // fluxesX stored in conserved variable form
            }
        }
    }

    //std::cout << std::endl << "-----------------------fluxesX above-------------------------" << std::endl;
}

void solver::MUSCL(bool direction,fluid& f, EOS* eos){
    
    if (PRIM == 1){
        //std::cout << "starting Primitive MUSCL" << std::endl;
        for (size_t i = 0; i < f.uPrim.size(); ++i){
            for (size_t j = 0; j < f.uPrim[0].size(); ++j){
                f.uPrim[i][j] = eos->consvToPrim(f.u[i][j]);
            }
        }
    } else {
        //std::cout << "starting conservative MUSCL" << std::endl;
    }

    calcHalfSlopes(PRIM,direction,f);

    calcr(direction,f);

    calcUBars(PRIM,direction,f,eos);

    updUBars(PRIM,direction,f,eos);

    if (direction == XDIR){
        #pragma omp parallel for collapse(2)
        for (size_t i=0; i < f.fluxesX.size(); ++i){
            for (size_t j = 0; j < f.fluxesX[0].size(); ++j){
                f.fluxesX[i][j] = HLLC(f.uBarRX[i][j],f.uBarLX[i][j+1],eos,i,j,XDIR);
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (size_t j=0; j < f.fluxesY[0].size(); ++j){
            for (size_t i = 0; i < f.fluxesY.size(); ++i){
                f.fluxesY[i][j] = HLLC(f.uBarRY[i][j],f.uBarLY[i+1][j],eos,i,j,YDIR);
            }
        }
    }
}

void solver::HLLCGodunov(bool direction,fluid& f, EOS* eos){
    if (direction == XDIR){
        for (size_t i=0; i < f.fluxesX.size(); ++i){
            for (size_t j = 0; j < f.fluxesX[0].size(); ++j){
                f.fluxesX[i][j] = HLLC(f.u[i+nG][j+nG-1],f.u[i+nG][j+nG],eos,i,j,XDIR);
            }
        }
    } else {
        for (size_t j=0; j < f.fluxesY[0].size(); ++j){
            for (size_t i = 0; i < f.fluxesY.size(); ++i){
                f.fluxesY[i][j] = HLLC(f.u[i+nG-1][j+nG],f.u[i+nG][j+nG],eos,i,j,YDIR);
            }
        }
    }
}

std::array<double,4> solver::HLLC(std::array<double,4> left,std::array<double,4> right, EOS* eos, int i, int j, bool direction){ // takes conserved variables
    double SL = 0.0, SR = 0.0;
    std::array<double,4> LPrim = eos->consvToPrim(left);
    std::array<double,4> RPrim = eos->consvToPrim(right);

    //printState(left); printState(right);
    //printState(LPrim); printState(RPrim);

    // Pressure - based wavespeed estimation
    /*
    double pEst = 0.5*(LPrim[PRES]+RPrim[PRES])-0.5*(RPrim[UX]-LPrim[UX])*0.25*(LPrim[RHO]+RPrim[RHO])*(eos->calcSoundSpeed(LPrim)+eos->calcSoundSpeed(RPrim));
    double qL,qR;
    
    qL = (pEst < LPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/LPrim[PRES])-1));
    qR = (pEst < RPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/RPrim[PRES])-1));

    SL = LPrim[UX] - eos->calcSoundSpeed(LPrim)*qL;
    SR = RPrim[UX] + eos->calcSoundSpeed(RPrim)*qR;
    */
    // easy wavespeed
    SL = (direction == XDIR) ? std::min(LPrim[UX] - eos->calcSoundSpeed(LPrim), RPrim[UX] - eos->calcSoundSpeed(RPrim)) : std::min(LPrim[UY] - eos->calcSoundSpeed(LPrim), RPrim[UY] - eos->calcSoundSpeed(RPrim));
    SR = (direction == XDIR) ? std::max(LPrim[UX] + eos->calcSoundSpeed(LPrim), RPrim[UX] + eos->calcSoundSpeed(RPrim)) : std::max(LPrim[UY] + eos->calcSoundSpeed(LPrim), RPrim[UY] + eos->calcSoundSpeed(RPrim));
    //std::cout << SL << " " << SR << std::endl;
    if ((std::isnan(SL) == 1) || (std::isnan(SR) == 1)){
        std::cout << "SL SR " << SL << " " << SR << std::endl;
    }

    double rhoL = left[RHO];
    double rhoR = right[RHO];

    // Calculate the numerator of sStar
    double numerator = (direction == XDIR) ? RPrim[PRES] - LPrim[PRES] + rhoL * LPrim[UX] * (SL - LPrim[UX]) - rhoR * RPrim[UX] * (SR - RPrim[UX])
                                        : RPrim[PRES] - LPrim[PRES] + rhoL * LPrim[UY] * (SL - LPrim[UY]) - rhoR * RPrim[UY] * (SR - RPrim[UY]);
    
    // Calculate the denominator of sStar
    double denominator = (direction == XDIR) ? rhoL * (SL - LPrim[UX]) - rhoR * (SR - RPrim[UX])
                                            : rhoL * (SL - LPrim[UY]) - rhoR * (SR - RPrim[UY]);

    
    // Calculate sStar
    double sStar = 0;
    if (numerator == 0 && denominator == 0){
        sStar = 0;
    } else if (denominator == 0){
        sStar = 1e15;
    } else {
        sStar = numerator / denominator;
    }
    
    
    if (std::isnan(sStar) == 1){
        std::cout << "numerator " <<  numerator << std::endl;
        std::cout << "denominator " << denominator << std::endl; 
        throw std::runtime_error("sStar is nan");
    }
    
    //
    double prefL = (direction == XDIR) ? (SL-LPrim[UX])/(SL-sStar) : (SL-LPrim[UY])/(SL-sStar);
    double prefR = (direction == XDIR) ? (SR-RPrim[UX])/(SR-sStar) : (SR-RPrim[UY])/(SR-sStar);

    double eL,eR;
    eL = (direction == XDIR) ? (left[ENE]/rhoL)+(sStar-LPrim[UX])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UX])))
                            : (left[ENE]/rhoL)+(sStar-LPrim[UY])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UY])));
    eR = (direction == XDIR) ? (right[ENE]/rhoR)+(sStar-RPrim[UX])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UX])))
                            : (right[ENE]/rhoR)+(sStar-RPrim[UY])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UY])));


    std::array<double,4> uHLLCL,uHLLCR;
    uHLLCL = (direction == XDIR) ? std::array<double,4>{left[RHO]*prefL,rhoL*prefL*sStar, rhoL*prefL*LPrim[UY], rhoL*prefL*eL}
                                : std::array<double,4>{left[RHO]*prefL, rhoL*prefL*LPrim[UX], rhoL*prefL*sStar, rhoL*prefL*eL}; // ?
    
    uHLLCR = (direction == XDIR) ? std::array<double,4>{right[RHO]*prefR, rhoR*prefR*sStar, rhoR*prefR*RPrim[UY], rhoR*prefR*eR}
                                : std::array<double,4>{right[RHO]*prefR, rhoR*prefR*RPrim[UX], rhoR*prefR*sStar, rhoR*prefR*eR}; // ?


    if (0 <= SL){
        return flux(LPrim,eos,direction);
    } else if (SL < 0 && 0 <= sStar){
        return flux(LPrim,eos,direction)+SL*(uHLLCL-left);
    } else if (sStar < 0 && 0 <= SR){
        return flux(RPrim,eos,direction)+SR*(uHLLCR-right);
    } else if (SR < 0){
        return flux(RPrim,eos,direction);
    } else {
        throw std::runtime_error("wavespeed condition invalid");
    }

}



void solver::pointsUpdate(bool direction, fluid& f){
    int fl1 = (&f == &fluid1) ? -1 : 1;
    if (direction == XDIR){
        #pragma omp parallel for collapse(2)
        for (size_t i=nG; i < f.u.size()-nG; ++i){
            for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
                //std::cout << (dt / dx) * (fluxesX[i-nG+1] - fluxesX[i-nG])[0] << std::endl;
                // if (std::signbit(phi[i][j]) == std::signbit(fl1)){
                //     f.uPlus1[i][j] = f.u[i][j] - (dt / dx) * (f.fluxesX[i-nG][j-nG+1] - f.fluxesX[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
                // } else {
                //     f.uPlus1[i][j] = f.u[i][j];
                // }
                f.uPlus1[i][j] = f.u[i][j] - (dt / dx) * (f.fluxesX[i-nG][j-nG+1] - f.fluxesX[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update

            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (size_t j=nG; j < f.u[0].size()-nG; ++j){
            for (size_t i = nG; i < f.u.size() - nG; ++i) {
                //std::cout << (dt / dx) * (fluxesX[i-nG+1] - fluxesX[i-nG])[0] << std::endl;
                // if (std::signbit(phi[i][j]) == std::signbit(fl1)){
                //     f.uPlus1[i][j] = f.u[i][j] - (dt / dy) * (f.fluxesY[i-nG+1][j-nG] - f.fluxesY[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
                // } else {
                //     f.uPlus1[i][j] = f.u[i][j];
                // }
                f.uPlus1[i][j] = f.u[i][j] - (dt / dy) * (f.fluxesY[i-nG+1][j-nG] - f.fluxesY[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update

            }
        }
    }
    f.u = f.uPlus1;
    setBCs(f);
    // print_arr(u,1);
    //std::cout << "updated u" << std::endl;
    // }
};

void solver::RK2(bool direction,fluid& f,EOS* eos){
    std::vector< std::vector< std::array<double,4> > > uStored = f.u;
    if (direction == XDIR){
        this->pointsUpdate(XDIR,f); // gives updated u
        f.uPlus1 = f.u; // stores new u as uPlus1
        this->fluxMethod(XDIR,f,eos); // calculates fluxes for new u

        #pragma omp parallel for collapse(2)
        for (size_t i = nG; i < f.u.size() - nG; ++i){
            for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
                f.uPlus1[i][j] = 0.5*(uStored[i][j]+f.uPlus1[i][j]) - 0.5 * (dt / dx) * (f.fluxesX[i-nG][j-nG+1] - f.fluxesX[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            }
        }
    } else {
        this->pointsUpdate(YDIR,f); // gives updated u
        f.uPlus1 = f.u; // stores new u as uPlus1
        this->fluxMethod(YDIR,f,eos); // calculates fluxes for new u

        #pragma omp parallel for collapse(2)
        for (size_t j = nG; j < f.u[0].size() - nG; ++j){
            for (size_t i = nG; i < f.u.size() - nG; ++i) {
                f.uPlus1[i][j] = 0.5*(uStored[i][j]+f.uPlus1[i][j]) - 0.5 * (dt / dx) * (f.fluxesY[i-nG+1][j-nG] - f.fluxesY[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            }
        }

    }
    f.u = f.uPlus1;
    setBCs(f);

}

// ---------------- Source Terms ---------------------------- //

/*
void solver::sourceUpdate(){
    setBCs();
    std::vector< std::array<double,4> > sourceArr = sourceTerm(u,alpha);
    //std::cout << "printing source array" << std::endl;
    // print_arr(sourceArr,0);

    for (std::vector<double>::size_type i = nG; i < u.size() - nG; ++i) {
        uPlus1[i] = u[i] + dt * sourceArr[i-nG]; // flux[i + 1] and flux[i] for the update

        //print_vect(uPlus1[i]);
    }
    
    //std::cout << std::endl << "-------------------------source terms above------------------------" << std::endl;
    // for (int i = 0; i < u.size(); ++i) {
    u = uPlus1;
    std::cout << "evolved source terms" << std::endl;
    // }
};
*/

// ---------------- Boundary Conditions --------------------- // 

void solver::transmissiveBC(fluid& fluid){
    for (int i = 0; i < nG; ++i){
        for (size_t j=0; j<fluid.uPlus1[0].size(); ++j){ // sets rows (i)
            fluid.uPlus1[i][j] = fluid.u[nG][j];
            fluid.u[i][j] = fluid.u[nG][j];
            fluid.uPlus1[nCellsY+2*nG-1-i][j] = fluid.u[nCellsY+nG-1][j];
            fluid.u[nCellsY+2*nG-1-i][j] = fluid.u[nCellsY+nG-1][j];
        }
    }
    for (size_t i=0; i<fluid.uPlus1.size(); ++i){
        for (int j = 0; j < nG; ++j){
            fluid.uPlus1[i][j] = fluid.u[i][nG];
            fluid.u[i][j] = fluid.u[i][nG];
            fluid.uPlus1[i][nCellsX+2*nG-1-j] = fluid.u[i][nCellsX+nG-1];
            fluid.u[i][nCellsX+2*nG-1-j] = fluid.u[i][nCellsX+nG-1];
        }
    }
}

// -------------------- Flux functions --------------------- //

std::array<double,4> solver::finvectLF(const std::array<double,4> x)const{
    double a=1;
    return a*x;
};

std::array<double,4> solver::fBurgersLF(const std::array<double,4> x)const{
    return 0.5*(x * x);
};

std::array<double,4> solver::fEuler(std::array<double,4> arr, EOS* eos, bool direction){ // takes primitive variables
    std::array<double,4> result;
    if (direction == XDIR){
        result[RHO] = arr[RHO]*arr[UX]; // rho1*v
        result[XMOM] = (arr[RHO])*arr[UX]*arr[UX] + arr[PRES]; // rho*v^2 + p
        result[YMOM] = (arr[RHO])*arr[UX]*arr[UY]; // rho*vx*vy
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UX]; // (E + p)v
        return result;
    } else {
        result[RHO] = arr[RHO]*arr[UY]; // rho1*v
        result[XMOM] = (arr[RHO])*arr[UX]*arr[UY]; // rho*vx*vy
        result[YMOM] = (arr[RHO])*arr[UY]*arr[UY] + arr[PRES]; // rho*v^2 + p
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UY]; // (E + p)v
        return result;
    }

}

// ------------------------- secondary fluxesX --------------------- //

std::array<double,4> solver::LF(std::array<double,4> v1, std::array<double,4> v2, EOS* eos,bool direction)const{
    if (direction == XDIR){
        return 0.5 * (flux(eos->consvToPrim(v1),eos,XDIR) + flux(eos->consvToPrim(v2),eos,XDIR)) + 0.5 * (dx / dt) * (v1 - v2);
    } else {
        return 0.5 * (flux(eos->consvToPrim(v1),eos,YDIR) + flux(eos->consvToPrim(v2),eos,YDIR)) + 0.5 * (dy / dt) * (v1 - v2);
    }
    
};

std::array<double,4> solver::laxFhalf(std::array<double,4> ui, std::array<double,4> uip1, EOS* eos, bool direction)const{
    if (direction == XDIR){
        return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos->consvToPrim(uip1),eos,XDIR)-flux(eos->consvToPrim(ui),eos,XDIR)));
    } else {
        return 0.5*(ui+uip1)-(0.5*(dt/dy)*(flux(eos->consvToPrim(uip1),eos,YDIR)-flux(eos->consvToPrim(ui),eos,YDIR)));
    }
};

// -------------- Level set ----------------------------------//

void solver::phiBC(){
    for (int i = 0; i < nG; ++i){
        for (size_t j=0; j<phi[0].size(); ++j){ // sets rows (i)
            phiPlus1[i][j] = phi[nG][j];
            phi[i][j] = phi[nG][j];
            phiPlus1[nCellsY+2*nG-1-i][j] = phi[nCellsY+nG-1][j];
            phi[nCellsY+2*nG-1-i][j] = phi[nCellsY+nG-1][j];
        }
    }
    for (size_t i=0; i<phi.size(); ++i){
        for (int j = 0; j < nG; ++j){
            phiPlus1[i][j] = phi[i][nG];
            phi[i][j] = phi[i][nG];
            phiPlus1[i][nCellsX+2*nG-1-j] = phi[i][nCellsX+nG-1];
            phi[i][nCellsX+2*nG-1-j] = phi[i][nCellsX+nG-1];
        }
    }
}

void solver::phiUpdate(int splitFlip){

    double uVal;

    phiOld = phi;
    if (splitFlip%2 == 0){
        for (std::vector<double>::size_type i = nG-1; i < phi.size()-nG+1; ++i){
            for (std::vector<double>::size_type j = nG; j < phi[0].size() - nG; ++j) {
                if (phi[i][j] <= 0){
                    uVal = eos[0]->consvToPrim(fluid1.u[i][j])[UX];
                } else {
                    uVal = eos[1]->consvToPrim(fluid2.u[i][j])[UX];
                }

                if (uVal >= 0){
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dx) * (phi[i][j] - phi[i][j-1]); // flux[i + 1] and flux[i] for the update
                } else {
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dx) * (phi[i][j+1] - phi[i][j]); // flux[i + 1] and flux[i] for the update
                }
            }
        }

        phi = phiPlus1;
        phiBC();
        //print_arr(u,RHO);
        //std::cout << "phiX updated" << std::endl;
        for (std::vector<double>::size_type i = nG; i < phi.size()-nG; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nG-1; j < phi[0].size() - nG+1; ++j) {
                if (phi[i][j] <= 0){
                    uVal = eos[0]->consvToPrim(fluid1.u[i][j])[UY];
                } else {
                    uVal = eos[1]->consvToPrim(fluid2.u[i][j])[UY];
                }
                
                if (uVal >= 0){
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dy) * (phi[i][j] - phi[i-1][j]); // flux[i + 1] and flux[i] for the update
                } else {
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dy) * (phi[i+1][j] - phi[i][j]); // flux[i + 1] and flux[i] for the update
                }
            }
        }

        phi = phiPlus1;
    } else {
        for (std::vector<double>::size_type i = nG; i < phi.size()-nG; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nG-1; j < phi[0].size() - nG+1; ++j) {
                if (phi[i][j] <= 0){
                    uVal = eos[0]->consvToPrim(fluid1.u[i][j])[UY];
                } else {
                    uVal = eos[1]->consvToPrim(fluid2.u[i][j])[UY];
                }
                
                if (uVal >= 0){
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dy) * (phi[i][j] - phi[i-1][j]); // flux[i + 1] and flux[i] for the update
                } else {
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dy) * (phi[i+1][j] - phi[i][j]); // flux[i + 1] and flux[i] for the update
                }
            }
        }

        phi = phiPlus1;
        phiBC();

        for (std::vector<double>::size_type i = nG-1; i < phi.size()-nG+1; ++i){
            for (std::vector<double>::size_type j = nG; j < phi[0].size() - nG; ++j) {
                if (phi[i][j] <= 0){
                    uVal = eos[0]->consvToPrim(fluid1.u[i][j])[UX];
                } else {
                    uVal = eos[1]->consvToPrim(fluid2.u[i][j])[UX];
                }

                if (uVal >= 0){
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dx) * (phi[i][j] - phi[i][j-1]); // flux[i + 1] and flux[i] for the update
                } else {
                    phiPlus1[i][j] = phi[i][j] - uVal * (dt / dx) * (phi[i][j+1] - phi[i][j]); // flux[i + 1] and flux[i] for the update
                }
            }
        }

        phi = phiPlus1;
    }

    // identify freshly cleared cells
    for (size_t i = 0; i<phiOld.size(); ++i){
        for (size_t j = 0; j < phiOld[0].size(); ++j){
            if (phi[i][j]*phiOld[i][j] < 0){
                std::array<int,2> coord = {i,j};
                freshlyCleared.push_back(coord);
            }
        }
    }

    phiBC();
}

void solver::neighbourAvg(fluid& f, int i, int j){

    std::array<double, 4> sum = {0.0, 0.0, 0.0, 0.0};
    int count = 0;

    // Get the sign of phi at (i, j)
    bool sign = phi[i][j] >= 0.0;

    // Offsets for the four neighbors
    int offsets[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

    // Check each neighbor
    for (int k = 0; k < 4; ++k) {
        int ni = i + offsets[k][0];
        int nj = j + offsets[k][1];

        // Boundary check
        if (ni >= 0 && ni < static_cast<int>(phi.size()) && nj >= 0 && nj < static_cast<int>(phi[0].size())) {
            // Check if the sign of phi is the same
            if ((phi[ni][nj] >= 0.0) == sign) {
                sum = sum + f.u[ni][nj];  // Sum the four-vectors
                ++count;
            }
        }
    }

    // Average and assign if there are valid neighbors
    if (count > 0) {
        f.u[i][j] = sum / count;  // Average the four-vectors
    } else {
        throw std::runtime_error("freshly cleared cell has no valid neighbours");
    }
}

std::array<double,4> solver::neighbourAvgInterp(fluid& f, int i, int j){

    std::array<double, 4> sum = {0.0, 0.0, 0.0, 0.0};
    int count = 0;

    // Get the sign of phi at (i, j)
    bool sign = phi[i][j] >= 0.0;

    // Offsets for the four neighbors
    int offsets[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

    // Check each neighbor
    for (int k = 0; k < 4; ++k) {
        int ni = i + offsets[k][0];
        int nj = j + offsets[k][1];

        // Boundary check
        if (ni >= 0 && ni < static_cast<int>(phi.size()) && nj >= 0 && nj < static_cast<int>(phi[0].size())) {
            // Check if the sign of phi is the same
            if ((phi[ni][nj] >= 0.0) == sign) {
                sum = sum + f.u[ni][nj];  // Sum the four-vectors
                ++count;
            }
        }
    }

    // Average and assign if there are valid neighbors
    if (count > 0) {
        return sum / count;  // Average the four-vectors
    } else {
        throw std::runtime_error("freshly cleared cell has no valid neighbours");
    }
}

void solver::fixFreshlyCleared(){
for (std::array<int,2> cell : freshlyCleared){
    neighbourAvg(fluid1,cell[0],cell[1]);
    neighbourAvg(fluid2,cell[0],cell[1]);
}
}

void solver::initialisePhi(){
    for (std::array<int,2> i : fluid1.interfaceCells){
        phi[i[0]][i[1]] = 0.0;
    }
    for (std::array<int,2> i : fluid2.interfaceCells){
        phi[i[0]][i[1]] = -0.0;
    }
    reinitPhiFastSweeping(20);
}

void solver::reinitPhiFastSweeping(int maxIter){
    //printInterfaceArrays(phi, fluid1, "interface1");
    //printInterfaceArrays(phi, fluid2, "interface2");

    for (int iter = 0; iter < maxIter; ++iter){

        for (int i = 1; i < phi.size()-1; ++i){
            for (int j = 1; j < phi[0].size()-1; ++j){
                solveEikonalPoint(i,j);
            }
        }
        for (int j = 1; j < phi[0].size()-1; ++j){
            for (int i = phi.size()-2; i >= 1; --i){
                solveEikonalPoint(i,j);
            }
        }
        for (int i = 1; i < phi.size()-1; ++i){
            for (int j = phi[0].size()-2; j >= 1; --j){
                solveEikonalPoint(i,j);
            }
        }
        for (int j = 1; j < phi[0].size()-1; ++j){
            for (int i = 1; i < phi.size()-1; ++i){
                solveEikonalPoint(i,j);
            }
        }
        phiBC();
    }
    phiBC();
}

void solver::solveEikonalPoint(int i, int j){
    double trialResult;
    if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
        && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
            
        double X{1.0/(dx*dx)}, Y{1.0/(dy*dy)};
        //std::cout << alpha << std::endl;
        double phiX,phiY;
        if (phi[i][j] <= 0){
            phiX = std::max(phi[i][j+1],phi[i][j-1]);
            phiY = std::max(phi[i+1][j],phi[i-1][j]);
        } else {
            phiX = std::min(phi[i][j+1],phi[i][j-1]);
            phiY = std::min(phi[i+1][j],phi[i-1][j]);
        }

        double A = (X+Y);
        double B = -2*(X*phiX + Y*phiY);
        double C = X*phiX*phiX + Y*phiY*phiY - 1;
        double D = B*B - 4*A*C;

        //std::cout << phiX << " " << phiY << " " << A << " " << B << " " << C << " " << D << std::endl;

        if (D<0){
            phiX = (phiX>phiY) ? 0 : phiX;  
            phiY = (phiY>phiX) ? 0 : phiY;
            A = (X+Y);
            B = -2*(X*phiX + Y*phiY);
            C = X*phiX*phiX + Y*phiY*phiY - 1;
            D = B*B - 4*A*C;

            if (D < 0){
                trialResult = (phi[i][j] <= 0) ? -1.0/sqrt(X+Y) : 1.0/sqrt(X+Y);
                //trialResult = phi[i][j];
            } else {
                trialResult = (phi[i][j] <= 0) ? ((-1.0)*sqrt(D)-B)/(2*A) : (sqrt(D)-B)/(2*A);
            }
        } else {
            trialResult = (phi[i][j] <= 0) ? ((-1.0)*sqrt(D)-B)/(2*A) : (sqrt(D)-B)/(2*A);
        }

        //std::cout << trialResult << std::endl;
        phi[i][j] = (std::abs(trialResult) < std::abs(phi[i][j])) ? trialResult : phi[i][j]; // "reinit 2"
        //phi[i][j] = trialResult; // "reinit 1"
    }
}


// Update the level-set value at an interior grid point (i,j).
// This routine computes a new candidate value for phi[i][j] by solving
// the quadratic equation derived from approximating |∇φ| = 1 using
// upwind differences. The neighbors are chosen based on the sign of φ.
void solver::updateLevelSetPoint(std::vector<std::vector<double>> &phi, int i, int j) {
    // Determine the sign of the current value
    int s = sign(phi[i][j]);

    // For reinitialization we choose neighboring values using upwinding.
    // For φ > 0, use the minimum of neighboring values;
    // for φ < 0, use the maximum.
    double a, b;
    if (s > 0) {
        b = std::min(phi[i - 1][j], phi[i + 1][j]);
        a = std::min(phi[i][j - 1], phi[i][j + 1]);
    } else {
        b = std::max(phi[i - 1][j], phi[i + 1][j]);
        a = std::max(phi[i][j - 1], phi[i][j + 1]);
    }

    // Compute coefficients for the quadratic equation:
    //   A * new_phi^2 + B * new_phi + C = 0,
    // where A, B, and C are derived from the finite-difference approximation:
    //   ( (φ_new - a)^2 / dx^2 ) + ( (φ_new - b)^2 / dy^2 ) = 1.
    double X = 1.0 / (dx * dx);
    double Y = 1.0 / (dy * dy);
    double A = X + Y;
    double B = -2.0 * (X * a + Y * b);
    double C = X * a * a + Y * b * b - 1.0;
    
    double disc = B * B - 4.0 * A * C;
    double phi_new;

    // If the discriminant is negative, fallback to a simple update.
    if (disc < 0) {
        // Here we simply take the neighbor with smaller absolute value and add a fixed increment.
        //phi_new = std::min(a, b) + s / sqrt(A);
        //phi_new = (s > 0) ? std::min(a, b) + 1.0 / sqrt(A) : std::max(a, b) - 1.0 / sqrt(A);
        if (fabs(a) > fabs(b)){
            phi_new = (s > 0) ? b+dy : b-dy;
        } else {
            phi_new = (s > 0) ? a+dx : a-dx;
        }

    } else {
        double sqrt_disc = sqrt(disc);
        // For φ > 0 choose the positive root; for φ < 0 choose the negative root.
        if (s > 0){
            phi_new = (-B + sqrt_disc) / (2.0 * A);
        } else {
            phi_new = (-B - sqrt_disc) / (2.0 * A);
        }
    }

    // Update only if the new value is closer to zero,
    // so that the solution converges toward a signed distance function.
    if (fabs(phi_new) < fabs(phi[i][j])) {
        phi[i][j] = phi_new;
    }
}

// Reinitializes the level set phi using a fast sweeping method.
// The function performs maxIter sweeps over the interior of the domain.
// Boundary values are kept fixed.
void solver::reinitialiseLevelSet(std::vector<std::vector<double>> &phi, int maxIter) {
    int m = phi.size();
    if (m == 0)
        return;
    int n = phi[0].size();

    #pragma omp parallel for collapse(2)
    for (int i = nG; i < m - nG; i++) {
        for (int j = nG; j < n - nG; j++) {
            if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
    && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
                phi[i][j] = (phi[i][j] <= 0) ? -1e9 : 1e9;
            }
        }
    }
    phiBC();



    // Loop over iterations
    for (int iter = 0; iter < maxIter; iter++) {
        // Sweep 1: top-left to bottom-right
        for (int i = nG; i < m - nG; i++) {
            for (int j = nG; j < n - nG; j++) {
                if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
        && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
                    updateLevelSetPoint(phi, i, j);
                }
            }
        }
        // Sweep 4: bottom-left to top-right
        for (int i = m - nG - 1; i >= nG; i--) {
            for (int j = nG; j < n - nG; j++) {
                if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
                && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
                    updateLevelSetPoint(phi, i, j);
                }
            }
        }
        // Sweep 2: bottom-right to top-left
        for (int i = m - nG - 1; i >= nG; i--) {
            for (int j = n - nG -1; j >= nG; j--) {
                if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
                    && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
                    updateLevelSetPoint(phi, i, j);
                }
            }
        }
        // Sweep 3: top-right to bottom-left
        for (int i = nG; i < m - nG; i++) {
            for (int j = n - nG - 1; j >= nG; j--) {
                if (std::signbit(phi[i][j]) == std::signbit(phi[i][j+1]) && std::signbit(phi[i][j]) == std::signbit(phi[i][j-1])
                && std::signbit(phi[i][j]) == std::signbit(phi[i+1][j]) && std::signbit(phi[i][j]) == std::signbit(phi[i-1][j])){
                    updateLevelSetPoint(phi, i, j);
                }
            }
        }

    }
    phiBC();
}


// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(bool prim, bool direction, fluid& f){
    if (direction == XDIR){
        #pragma omp parallel for collapse(2)
        for (size_t i=0; i < f.halfSlopesX.size(); ++i){
            for (size_t j = 0; j < f.halfSlopesX[0].size(); ++j){ // halfslopes[i] = halfslopes i-1/2
                f.halfSlopesX[i][j] = (prim == 0) ? f.u[i+nG][j+1]-f.u[i+nG][j] : f.uPrim[i+nG][j+1]-f.uPrim[i+nG][j];
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (size_t j=0; j < f.halfSlopesY[0].size(); ++j){
            for (size_t i = 0; i < f.halfSlopesY.size(); ++i){ // halfslopes[i] = halfslopes i-1/2
                f.halfSlopesY[i][j] = (prim == 0) ? f.u[i+1][j+nG]-f.u[i][j+nG] : f.uPrim[i+1][j+nG]-f.uPrim[i][j+nG];
            }
        }
    }

};

void solver::calcr(bool direction,fluid& f){
    if (direction == XDIR){
        #pragma omp parallel for collapse(2)
        for (size_t i=0; i < f.rX.size(); ++i){
            for (size_t j = 0; j < f.rX[0].size(); ++j){
                f.rX[i][j] = elementDivide(f.halfSlopesX[i][j],f.halfSlopesX[i][j+1]);
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (size_t j=0; j < f.rY[0].size(); ++j){
            for (size_t i = 0; i < f.rY.size(); ++i){
                f.rY[i][j] = elementDivide(f.halfSlopesY[i][j],f.halfSlopesY[i+1][j]);
            }
        }
    }
    //print_arr(r);
    //std::cout << "calcd r" << std::endl;
};


/*

+---------------+--------------+--------------+
|   u[0][0]     |   u[0][1]    |    u[0][2]   |  <- Top row (Row 0)
|               |              |  (u[0][nG])  |
+---------------+--------------+--------------+
|               |              |              |  <- Middle row (Row 1)
|               |              |              |
+---------------+--------------+--------------+
|               |  uBar[0][0]  |   u[nG][nG]  |  <- Bottom row (Row 2)
|        halfSlopesX[0][0]  f[0][0](&sStarsX)   |
|               |   r[0][0]    |              |
+---------------+--------------+--------------+

*/


void solver::calcUBars(bool prim, bool direction, fluid& f, EOS* eos){ // calculates for all u except leftmost cell
    if (direction == XDIR){
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < f.uBarLX.size(); ++i){
            for (size_t j = 0; j < f.uBarLX[0].size(); ++j){
                f.uBarLX[i][j] = (prim == 0) ? f.u[i+nG][j+nG-1] - 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1])) : f.uPrim[i+nG][j+nG-1] - 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1]));
                f.uBarRX[i][j] = (prim == 0) ? f.u[i+nG][j+nG-1] + 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1])) : f.uPrim[i+nG][j+nG-1] + 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1]));
            }
        }
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < f.sStarsX.size(); ++i){
            for (size_t j = 0; j < f.sStarsX[0].size(); ++j){
                f.sStarsX[i][j] = (prim == 0) ? 0.5*(eos->consvToPrim(f.uBarRX[i][j])[UX]+eos->consvToPrim(f.uBarLX[i][j+1])[UX]) : 0.5*(f.uBarRX[i][j][UX]+f.uBarLX[i][j+1][UX]);
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (size_t j = 0; j < f.uBarLY[0].size(); ++j){
            for (size_t i = 0; i < f.uBarLY.size(); ++i){
                f.uBarLY[i][j] = (prim == 0) ? f.u[i+nG-1][j+nG] - 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j])) : f.uPrim[i+nG-1][j+nG] - 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j]));
                f.uBarRY[i][j] = (prim == 0) ? f.u[i+nG-1][j+nG] + 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j])) : f.uPrim[i+nG-1][j+nG] + 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j]));
            }
        }
        #pragma omp parallel for collapse(2)
        for (size_t j = 0; j < f.sStarsY[0].size(); ++j){
            for (size_t i = 0; i < f.sStarsY.size(); ++i){
                f.sStarsY[i][j] = (prim == 0) ? 0.5*(eos->consvToPrim(f.uBarRY[i][j])[UY]+eos->consvToPrim(f.uBarLY[i+1][j])[UY]) : 0.5*(f.uBarRY[i][j][UY]+f.uBarLY[i+1][j][UY]);
            }
        }
    }

};

std::array<double,4> multiplyMatrixVector(
    const std::vector<std::vector<double>>& B,
    const std::array<double,4>& x) {
    
    std::array<double,4> result = {0.0, 0.0, 0.0, 0.0};
    
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            result[i] += B[i][j] * x[j];
        }
    }
    
    return result;
}

void solver::updUBars(bool prim,bool direction, fluid& f, EOS* eos){
    if (direction == XDIR){
        if (prim == false){
            #pragma omp parallel for collapse(2)
            for (size_t i = 0; i<f.uBarLupdX.size(); ++i){
                for (size_t j = 0; j<f.uBarLupdX[0].size(); ++j){
                    f.uBarLupdX[i][j] = f.uBarLX[i][j]-0.5*(dt/dx)*(flux(eos->consvToPrim(f.uBarRX[i][j]),eos,XDIR)-flux(eos->consvToPrim(f.uBarLX[i][j]),eos,XDIR));
                    f.uBarRupdX[i][j] = f.uBarRX[i][j]-0.5*(dt/dx)*(flux(eos->consvToPrim(f.uBarRX[i][j]),eos,XDIR)-flux(eos->consvToPrim(f.uBarLX[i][j]),eos,XDIR));
                }
            }
        } else {
            #pragma omp parallel for collapse(2)
            for (size_t i = 0; i<f.uBarLupdX.size(); ++i){
                for (size_t j = 0; j<f.uBarLupdX[0].size(); ++j){
                    double v = f.uPrim[i+nG][j+nG-1][UX]; double rho = f.uPrim[i+nG][j+nG-1][RHO];
                    std::vector< std::vector <double> > B = {   {v, rho,  0,  0     },
                            {0, v,  0,  1.0/rho },
                            {0, 0,  v,  0       },
                            {0, rho*pow(eos->calcSoundSpeed(f.uPrim[i+nG][j+nG-1]),2), 0, v}     };

                    std::array<double,4> deltaU = f.uBarRX[i][j]-f.uBarLX[i][j];
                    f.uBarLupdX[i][j] = eos->primToConsv(f.uBarLX[i][j]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
                    f.uBarRupdX[i][j] = eos->primToConsv(f.uBarRX[i][j]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
                }
            }

        }
        f.uBarLX = f.uBarLupdX;
        f.uBarRX = f.uBarRupdX;
    } else {
        if (prim == false){
            #pragma omp parallel for collapse(2)
            for (size_t j = 0; j<f.uBarLupdY[0].size(); ++j){
                for (size_t i = 0; i<f.uBarLupdY.size(); ++i){
                    f.uBarLupdY[i][j] = f.uBarLY[i][j]-0.5*(dt/dy)*(flux(eos->consvToPrim(f.uBarRY[i][j]),eos,YDIR)-flux(eos->consvToPrim(f.uBarLY[i][j]),eos,YDIR));
                    f.uBarRupdY[i][j] = f.uBarRY[i][j]-0.5*(dt/dy)*(flux(eos->consvToPrim(f.uBarRY[i][j]),eos,YDIR)-flux(eos->consvToPrim(f.uBarLY[i][j]),eos,YDIR));
                }
            }
        } else {
            #pragma omp parallel for collapse(2)
            for (size_t j = 0; j<f.uBarLupdY[0].size(); ++j){
                for (size_t i = 0; i<f.uBarLupdY.size(); ++i){
                    double v = f.uPrim[i+nG-1][j+nG][UY]; double rho = f.uPrim[i+nG-1][j+nG][RHO];
                    std::vector< std::vector <double> > B = {   {v, 0, rho, 0},
                            {0, v,  0,  0},
                            {0, 0,  v,  1.0/rho},
                            {0, 0,  rho*pow(eos->calcSoundSpeed(f.uPrim[i+nG-1][j+nG]),2), v}     };

                    std::array<double,4> deltaU = f.uBarRY[i][j]-f.uBarLY[i][j];
                    f.uBarLupdY[i][j] = eos->primToConsv(f.uBarLY[i][j]-0.5*(dt/dy)*multiplyMatrixVector(B,deltaU));
                    f.uBarRupdY[i][j] = eos->primToConsv(f.uBarRY[i][j]-0.5*(dt/dy)*multiplyMatrixVector(B,deltaU));
                }
            }

        }
        f.uBarLY = f.uBarLupdY;
        f.uBarRY = f.uBarRupdY;
    }
};

std::array<double,4> solver::calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,4> solver::minbee(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
    double minb;
    for (int k = 0; k<4; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minb = 0;
        } else if (slp[k] > 1){
            minb = std::min(1.0,2.0/(1+slp[k]));
        } else {
            minb = slp[k];
        }
        minbArr[k] = minb;
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    //return {minBval,minBval,minBval,minBval};
    return minbArr;
};

std::array<double,4> solver::superbee(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
    double minb;
    for (int k = 0; k<4; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minb = 0;
        } else if (slp[k] <= 0.5){
            minb = 2*slp[k];
        } else if (slp[k] <= 1.0){
            minb = 1.0;
        } else {
            minb = std::min(slp[k],std::min(2.0/(1+slp[k]),2.0));
        }
        minbArr[k] = minb;
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    //return {minBval,minBval,minBval};
    return minbArr;
};

std::array<double,4> solver::vanLeer(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
    for (int k = 0; k<4; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minbArr[k] = 0;
        } else {
            minbArr[k] = std::min(2.0*slp[k]/(1.0+slp[k]),2.0/(1+slp[k]));
        }
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    return minbArr;
    //return {minBval,minBval,minBval,minBval};
};


// ---------------------- Misc ------------------------------ //

std::array<double,4> solver::set_vals(int fluid1, double rho, double vx, double vy, double pres){
    std::array<double,4> result;
    if (fluid1 == 1){
        result = eos[0]->primToConsv(std::array<double,4>{rho,vx,vy,pres});
    } else if (fluid1 == 2) {
        result = eos[1]->primToConsv(std::array<double,4>{rho,vx,vy,pres});
    } else {
        throw std::runtime_error("fluid1 must be 1 or 2");
    }
    return result;
};

void solver::print_arr(std::vector< std::array<double,4> > arr, int var){
    for (std::vector<double>::size_type i = 0; i < arr.size(); ++i) {
        std::cout << arr[i][var] << " ";
    }
    std::cout << std::endl;
}

void solver::print_vect(std::array<double,4> v){
    for (int j = 0; j<4; ++j){
        std::cout << v[j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::setDt(){
    double aMax = -1e14;
    std::cout << "setting dt" << std::endl;
    for (size_t i = nG; i < fluid1.u.size()-nG; ++i) {
        for (size_t j = nG; j < fluid1.u[0].size()-nG; ++j) {
            //aMax = std::max(aMax, std::sqrt(pow(eos[0]->consvToPrim(f.u[i][j])[UX],2)+pow(eos[0]->consvToPrim(f.u[i][j])[UY],2)) + eos[0]->calcSoundSpeed(eos[0]->consvToPrim(f.u[i][j])));
            if (phi[i][j] <= 0){
                aMax = std::max(aMax, std::sqrt(pow(eos[0]->consvToPrim(fluid1.u[i][j])[UX],2)+pow(eos[0]->consvToPrim(fluid1.u[i][j])[UY],2))+eos[0]->calcSoundSpeed(eos[0]->consvToPrim(fluid1.u[i][j])));
            } else {
                aMax = std::max(aMax, std::sqrt(pow(eos[1]->consvToPrim(fluid2.u[i][j])[UX],2)+pow(eos[1]->consvToPrim(fluid2.u[i][j])[UY],2))+eos[1]->calcSoundSpeed(eos[1]->consvToPrim(fluid2.u[i][j])));
            }
        }
    }

    dtReducer = (splitFlip < 20) ? 0.5 : 1.0;

    dt = dtReducer * cour * std::min(dx,dy) / std::fabs(aMax);

    std::cout << "set dt" << std::endl;

    std::cout << "aMax =" << aMax << std::endl;


    double writeTime = writeInterval*(std::floor((time+1e-9) / writeInterval)+1);

    if (endTime - time <= dt && writeTime-time >= dt){
        dt = endTime - time;
        assert(dt != 0);
        time = endTime;
        std::cout << "last step:" << dt << std::endl;
    } else if (writeTime-time <= dt && writeTime-time > 0){ // manually reducing dt to hit the write times
        dt = writeTime - time;
        assert(dt != 0);
        time = writeTime;
    } else {
        time = time + dt;
    }

    if (std::fabs(writeTime - time) < 1e-9){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    // std::cout << "dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double y_0, double y_1, double t0, double t1, int n_CellsX, int n_CellsY, int n_Ghosts, double c)
    : fluid1(fluid(n_CellsX, n_CellsY, n_Ghosts)), fluid2(fluid(n_CellsX, n_CellsY, n_Ghosts)), bound(boundary(n_CellsX,n_CellsY,n_Ghosts)),
        nCellsX(n_CellsX), nCellsY(n_CellsY), nG(n_Ghosts), 
        x0(x_0), x1(x_1),
        y0(y_0), y1(y_1),
        startTime(t0), endTime(t1),
        
        cour(c){
    
    assert(t1>t0);
    timeMulti = ((endTime-startTime) < 1e-2) ? 1000 : 1; // 1.0/(endTime-startTime);

    dx = (x1 - x0)/nCellsX;
    dy = (y1 - y0)/nCellsY;
    
    resizeVector(phi,nCellsY+2*nG,nCellsX+2*nG);
    resizeVector(phiPlus1,nCellsY+2*nG,nCellsX+2*nG);
    
    resizeVector(phiNormals,nCellsY+2*nG,nCellsX+2*nG);
    resizeVector(phiGrads,nCellsY+2*nG,nCellsX+2*nG);

    slopeWeight = 1;

    cour = c;

    std::vector<std::string> variables = {"rho","vx","vy","p"};

    //varMap["rho"] = 0; varMap["v"] = 1; varMap["p"] = 2;
    for (size_t i = 0; i<variables.size(); ++i){
        varMap[variables[i]] = i;
    }
    
};

// ---------------------- utilities -------------------------- //

int solver::ghosts(){
    return nG;
};

std::array<double,2> solver::get_dxdy()const{
    std::array<double,2> result = {dx,dy};
    return result;
};


void solver::init(std::vector < std::vector< std::array<double,4> > > init, fluid& f) {
    assert(init.size() == f.u.size());
    assert(init[0].size() == f.u[0].size());
    f.u = init;
};

fluid& solver::get_fluid(bool f){
    if (f==0){
        return fluid1;
    } else if (f==1){
        return fluid2;
    } else {
        throw std::invalid_argument("Invalid fluid index");
    }
}

void solver::phiInit(std::vector< std::vector < double > > init) {
    assert(init.size() == phi.size());
    phi = init;
};

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData() const{
        // making data file
        std::string filename = dirName + "/" + std::to_string(time*timeMulti).substr(0,5);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        writeFile << "x y phi rho1 rho2 vx1 vx2 vy1 vy2 p1 p2" << std::endl;
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (size_t i=0; i < fluid1.u.size();i++){
                double y = y0 + (i-static_cast<double>(nG)+0.5)*dy;
                for (size_t j=0; j < fluid1.u[0].size();j++){
                    double x = x0 + (j-static_cast<double>(nG)+0.5)*dx;
                    std::array<double,4> prim1 = eos[0]->consvToPrim(fluid1.u[i][j]);
                    std::array<double,4> prim2 = eos[1]->consvToPrim(fluid2.u[i][j]);
                    writeFile << x << " " << y << " " << phi[i][j] << " "
                    << prim1[RHO] << " " << prim2[RHO] 
                    << " " << prim1[UX] << " " << prim2[UX]
                     << " " << prim1[UY] << " " << prim2[UY]
                      << " " << prim1[PRES] << " " << prim2[PRES]
                       << std::endl;
                    //std::cout << u[i] << " ";
                }
                writeFile << "\n";
            }
            writeFile.close();
            std::cout << filename << "written successfully" << std::endl;
        }else{
            std::cout << "failed to write." << std::endl;
        }

}



