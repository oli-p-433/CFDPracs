#include "solver.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>

#include <thread>
#include <chrono>


void solver::run(){
    //std::cout << "pre-bc" << std::endl;

    transmissiveBC(fluid1);
    transmissiveBC(fluid2);
    phiBC();

    //double realTime = 0.0;


    splitFlip = 0;
    do{

        if (splitFlip < (endTime-startTime)/20){
            dtReducer = 0.01;
        } else {
            dtReducer = 1.0;
        }
        //time = 0+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        
        
        this->findBoundary();
        //std::cout << "number of interface cells = " << interfaceCells.size() << std::endl;
        //std::cout << "boundary found" << std::endl;

        //printInterfaceArray(phi, interfaceCells, "interface.txt");

        this->calcPhiGrads();
        //std::cout << "phi grads calculated" << std::endl;
        this->calcInterface();
        //std::cout << "interface calculated" << std::endl;
        this->interpInterfaceStates();
        //std::cout << "interface states interpolated" << std::endl;
        this->resolveVelocities();
        //std::cout << "resolved velocities" << std::endl;
        this->interfaceRiem();
        //std::cout << "interface riem" << std::endl;
        this->calcStarStates();
        //std::cout << "star states calculated" << std::endl;

        //time = 1+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        this->setGhostFluids();
        //std::cout << "set ghost fluids" << std::endl; 

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);

        //time = 2+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        //std::this_thread::sleep_for(std::chrono::seconds(5));
        //time = realTime;
        
        this->setDt();
        
        //realTime = time;

        std::cout << "time is " << time << std::endl;

        //std::cout << "setting phiBC" << std::endl;
        phiBC();
        //std::cout << "updating phi" << std::endl;

        //time = 3+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        splitFlip++;
        phiUpdate(splitFlip);
        std::cout << "phi updated" << std::endl;

        //time = 4+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        fixFreshlyCleared();
        std::cout << "freshly cleared fixed" << std::endl;

        //time = 5+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        //time = 0.0;

        if (splitFlip % 10 == 0){
            reinitPhi();
            std::cout << "phi reinitialised" << std::endl;
        }
        
        
        std::cout << "Flux calc & point updates" << std::endl;
        if (splitFlip%2 == 0){
            //std::cout << "x first" << std::endl;
            direction = XDIR;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
            direction = YDIR;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
        } else {
            //std::cout << "y first" << std::endl;
            direction = YDIR;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
            direction = XDIR;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
        }
        //std::cout << "points updated" << std::endl;

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);

        this->interfaceCells={};
        this->freshlyCleared={};

        //time = 6+((splitFlip-1)*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        //throw std::runtime_error("simulation stopped");


        //this->sourceUpdate();

        if (checkWrite==1){
            std::cout << "------------------ writing " << time << " ------------------"<< std::endl;
            writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        }

    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
};

void solver::findBoundary(){ // checks for change in sign of level set. emplaces indices of cell with phi <= 0.
    for (size_t i = nGhost; i < phi.size()-nGhost; ++i) {  // Loop through all cells (except last)
        for (size_t j = nGhost; j < phi.size()-nGhost; ++j){
            if(phi[i][j] == 0){
                std::array<int,2> cell = {i,j};
                interfaceCells.push_back(cell);
            } else if ((std::signbit(phi[i][j]) != std::signbit(phi[i+1][j])) || 
                (std::signbit(phi[i][j]) != std::signbit(phi[i][j+1]))){  // Sign change detected
                //std::cout << phi[i] << " " << phi[i+1] << std::endl;
                if (std::signbit(phi[i][j]) == 1){
                    std::array<int,2> cell = {i,j};
                    interfaceCells.push_back(cell);
                } else if (std::signbit(phi[i][j]) == 0){
                    if (std::signbit(phi[i+1][j]) == 1){
                        std::array<int,2> cell = {i+1,j};
                        interfaceCells.push_back(cell);
                    }
                    if (std::signbit(phi[i][j+1]) == 1){
                        std::array<int,2> cell = {i,j+1};
                        interfaceCells.push_back(cell);
                    }
                }
                //std::cout << "interface at cell " << i << std::endl; 
            }
        }
    };
    int n = phi.size()-nGhost-1;
    if (phi[n-1][n] * phi[n][n] < 0){ // corner cases
        std::array<int,2> cell = (phi[n][n]>0) ? std::array<int,2>{n-1,n} : std::array<int,2>{n,n};
        interfaceCells.push_back(cell);
    }
    if (phi[n][n-1] * phi[n][n] < 0){
        std::array<int,2> cell = (phi[n][n]>0) ? std::array<int,2>{n,n-1} : std::array<int,2>{n,n};
        interfaceCells.push_back(cell);
    }
    interfaceStates.resize(interfaceCells.size());
    interfaceNormals.resize(interfaceCells.size());
    interfacePositions.resize(interfaceCells.size());
    resolvedVelocities.resize(interfaceCells.size());
    riemInterfaceStates.resize(interfaceCells.size());
    starStates.resize(interfaceCells.size());
    resize2Dc(nCells+2*nGhost,nCells+2*nGhost,uExtrap);

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
    printInterfaceArray(phi, interfaceCells, "interface.txt");

}

void solver::printInterfaceArray(std::vector<std::vector<double>> field, std::vector<std::array<int, 2>>& interfaceCells, const std::string& filename) {
    std::vector<std::vector<double>> interfaceArray = field;

    // Set all elements to zero
    for (auto& row : interfaceArray) {
        std::fill(row.begin(), row.end(), 0.0);
    }

    // Set the specified cells to 1
    for (const auto& cell : interfaceCells) {
        interfaceArray[cell[0]][cell[1]] = 1.0;
    }

    // Write the array to a file
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        for (const auto& row : interfaceArray) {
            for (const auto& value : row) {
                outFile << std::fixed << std::setprecision(2) << value << " ";
            }
            outFile << std::endl;
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void solver::calcPhiGrads(){
    //std::cout << std::endl;
    //printScalarField(phi);
    
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        int x = interfaceCells[i][1];
        int y = interfaceCells[i][0];
        double phi_y = (phi[y+1][x]-phi[y-1][x])/(2*dy);
        double phi_x = (phi[y][x+1]-phi[y][x-1])/(2*dx);

        std::array<double,2> grad = {phi_x/(sqrt(phi_x*phi_x+phi_y*phi_y)),phi_y*(sqrt(phi_x*phi_x+phi_y*phi_y))};
        interfaceNormals[i] = grad;
    }

    /*
    std::cout << std::endl; std::cout << "interface normals" << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        printPrec(interfaceNormals[i][0],2);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        printPrec(interfaceNormals[i][1],2);
    }
    */
    

    
}

void solver::calcInterface(){
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        double x_i = x0 + (interfaceCells[i][1]-1.5)*dx; // -2 for ghost cells, +0.5 for cell centre
        double y_i = y0 + (interfaceCells[i][0]-1.5)*dx;

        double x = x_i - phi[interfaceCells[i][0]][interfaceCells[i][1]]*interfaceNormals[i][0];
        double y = y_i - phi[interfaceCells[i][0]][interfaceCells[i][1]]*interfaceNormals[i][1];

        std::array<double,2> lPos = {x-1.5*dx*interfaceNormals[i][0],y-1.5*dy*interfaceNormals[i][1]};
        std::array<double,2> rPos = {x+1.5*dx*interfaceNormals[i][0],y+1.5*dy*interfaceNormals[i][1]};
        interfacePositions[i] = {lPos,rPos};
    }

    // optional - print positions
    /*
    std::cout << std::endl; std::cout << "left interface positions" << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        printPrec(interfacePositions[i][0][0],2);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        printPrec(interfacePositions[i][0][1],2);
    }
    */

}

void solver::interpInterfaceStates(){
    for (size_t i = 0; i < interfacePositions.size(); ++i){
        std::array<double,4> lState = bilinear(interfacePositions[i][0][0],interfacePositions[i][0][1],fluid1);
        std::array<double,4> rState = bilinear(interfacePositions[i][1][0],interfacePositions[i][1][1],fluid2);
        
        interfaceStates[i] = {lState,rState};
    }

    // optional - print states
    /*
    std::cout << std::endl; std::cout << "L/R rho" << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        printPrec(interfaceStates[i][0][0],3);
    }
    std::cout << std::endl;
    for (size_t i = 0; i < interfaceCells.size(); ++i){
        printPrec(interfaceStates[i][1][0],3);
    }
    */
    
}

void solver::resolveVelocities(){
    for (size_t i=0; i<interfaceStates.size(); ++i){
        double vNormL = interfaceStates[i][0][1]*interfaceNormals[i][0]+interfaceStates[i][0][2]*interfaceNormals[i][1];
        double vNormR = interfaceStates[i][1][1]*interfaceNormals[i][0]+interfaceStates[i][1][2]*interfaceNormals[i][1];
        double vTanL = -interfaceStates[i][0][1]*interfaceNormals[i][1]+interfaceStates[i][0][2]*interfaceNormals[i][0];
        double vTanR = -interfaceStates[i][1][1]*interfaceNormals[i][1]+interfaceStates[i][1][2]*interfaceNormals[i][0];

        std::array<double,2> vL = {vNormL,vTanL};
        std::array<double,2> vR = {vNormR,vTanR};
        resolvedVelocities[i] = {vL,vR};
    }
}

void solver::interfaceRiem(){
    for (size_t i=0; i<interfaceStates.size(); ++i){
        std::array<double,4> lState = {interfaceStates[i][0][0],resolvedVelocities[i][0][0],resolvedVelocities[i][0][1],interfaceStates[i][0][3]};  
        std::array<double,4> rState = {interfaceStates[i][1][0],resolvedVelocities[i][1][0],resolvedVelocities[i][1][1],interfaceStates[i][1][3]};

        riemann solution(eos[0]->get_gamma(),eos[1]->get_gamma(),eos[0]->consvToPrim(lState),eos[1]->consvToPrim(rState),1,-1.5*dx,1.5*dx,0,dt,2,0,0); // need p_inf setup
        //std::cout << "riemann solver object initialised" << std::endl;
        std::array<double,4> result = solution.exctRiemann();
        //if (sqrt(result[UX]*result[UX]) < 1e-4){result[UX] = 0.0;}
        //if (sqrt(result[UY]*result[UY]) < 1e-4){result[UY] = 0.0;}
        riemInterfaceStates[i] = result;
        //solution.exctRiemann();
    }
}

void solver::calcStarStates(){
    for (size_t i=0; i<interfaceStates.size(); ++i){
        // riemInterfaceStates[i] = {rho,v*_n,v*_t,p}
        std::array<double,4> sL = interfaceStates[i][0];
        std::array<double,4> sR = interfaceStates[i][1];
        double utLx = sL[UX]-(sL[UX]*interfaceNormals[i][0]+sL[UY]*interfaceNormals[i][1])*interfaceNormals[i][0];
        double utLy = sL[UY]-(sL[UX]*interfaceNormals[i][0]+sL[UY]*interfaceNormals[i][1])*interfaceNormals[i][1];
        double utRx = sR[UX]-(sR[UX]*interfaceNormals[i][0]+sR[UY]*interfaceNormals[i][1])*interfaceNormals[i][0];
        double utRy = sR[UY]-(sR[UX]*interfaceNormals[i][0]+sR[UY]*interfaceNormals[i][1])*interfaceNormals[i][1];
        
        double uLStarx = riemInterfaceStates[i][1]*interfaceNormals[i][0]+utLx;
        double uLStary = riemInterfaceStates[i][1]*interfaceNormals[i][1]+utLy;
        double uRStarx = riemInterfaceStates[i][1]*interfaceNormals[i][0]+utRx;
        double uRStary = riemInterfaceStates[i][1]*interfaceNormals[i][1]+utRy;

        std::array<double,4> uStarL = {riemInterfaceStates[i][0],uLStarx,uLStary,riemInterfaceStates[i][3]};
        std::array<double,4> uStarR = {riemInterfaceStates[i][0],uRStarx,uRStary,riemInterfaceStates[i][3]};
        std::array<std::array<double,4>,2> res = {uStarL,uStarR};
        starStates[i] = res;
        //std::cout << "star states left: " << uStarL[0] << " " << uStarL[1] << " " << uStarL[2] << " " << uStarL[3] << std::endl;
        //std::cout << "star states right: " << uStarR[0] << " " << uStarR[1] << " " << uStarR[2] << " " << uStarR[3] << std::endl;

    }
}

void solver::calcPhiNormals(){
    resize2Db(nCells+2*nGhost,nCells+2*nGhost,phiNormals);
    resize2D(nCells+2*nGhost,nCells+2*nGhost,phiGrads);
    for (int i = nGhost; i < nCells+nGhost; ++i){
        for (int j = nGhost; j < nCells+nGhost; ++j){
            double phi_y = (phi[i+1][j]-phi[i-1][j])/(2*dy);
            double phi_x = (phi[i][j+1]-phi[i][j-1])/(2*dx);
            std::array<double,2> norm = {phi_x/(sqrt(phi_x*phi_x+phi_y*phi_y)),phi_y*(sqrt(phi_x*phi_x+phi_y*phi_y))};
            phiNormals[i][j] = norm;
            phiGrads[i][j] = sqrt(phi_x*phi_x+phi_y*phi_y);
        }
    }
}

void solver::setInterface(){ // ind = 0 for fluid1, 1 for fluid2

    // setting uExtrap to the ghost fluids
    //std::cout << "no of star states = " << starStates.size() << std::endl;
    for (size_t i = 0; i<starStates.size(); ++i){
        fluid2.u[interfaceCells[i][0]][interfaceCells[i][1]] = eos[1]->primToConsv(starStates[i][1]); // interface states set to u*R
        if (phi[interfaceCells[i][0]+1][interfaceCells[i][1]] > 0){ // surrounding states set to u*L
            fluid1.u[interfaceCells[i][0]+1][interfaceCells[i][1]] = eos[0]->primToConsv(starStates[i][0]);
        }
        if (phi[interfaceCells[i][0]][interfaceCells[i][1]+1] > 0){
            fluid1.u[interfaceCells[i][0]][interfaceCells[i][1]+1] = eos[0]->primToConsv(starStates[i][0]);
        }
        if (phi[interfaceCells[i][0]-1][interfaceCells[i][1]] > 0){
            fluid1.u[interfaceCells[i][0]-1][interfaceCells[i][1]] = eos[0]->primToConsv(starStates[i][0]);
        }
        if (phi[interfaceCells[i][0]][interfaceCells[i][1]-1] > 0){
            fluid1.u[interfaceCells[i][0]][interfaceCells[i][1]-1] = eos[0]->primToConsv(starStates[i][0]);
        }
    }
    /*
    for (size_t i = 0; i<starStates.size(); ++i){
        //fluid1.u[interfaceCells[i][0]][interfaceCells[i][1]] = eos[0]->primToConsv(starStates[i][0]); // interface states set to u*L
        if (phi[interfaceCells[i][0]+1][interfaceCells[i][1]] <= 0){ // surrounding states set to u*R
            fluid2.u[interfaceCells[i][0]+1][interfaceCells[i][1]] = eos[1]->primToConsv(starStates[i][1]);
        } else if (phi[interfaceCells[i][0]][interfaceCells[i][1]+1] <= 0){
            fluid2.u[interfaceCells[i][0]][interfaceCells[i][1]+1] = eos[1]->primToConsv(starStates[i][1]);
        } else if (phi[interfaceCells[i][0]-1][interfaceCells[i][1]] <= 0){
            fluid2.u[interfaceCells[i][0]-1][interfaceCells[i][1]] = eos[1]->primToConsv(starStates[i][1]);
        } else if (phi[interfaceCells[i][0]][interfaceCells[i][1]-1] <= 0){
            fluid2.u[interfaceCells[i][0]][interfaceCells[i][1]-1] = eos[1]->primToConsv(starStates[i][1]);
        } else {
            std::cout << "no surrounding states found for interface cell " << interfaceCells[i][0] << " " << interfaceCells[i][1] << std::endl;
        }
    }
    */
}

void solver::calcnDotPhiNormals(){
    resize2Dc(nCells+2*nGhost,nCells+2*nGhost,nDotGradPhi);
    for (int i = nGhost; i < nCells+nGhost; ++i){
        for (int j = nGhost; j < nCells+nGhost; ++j){
            std::array<double,4> xCom{{0,0,0,0}}, yCom{{0,0,0,0}};
            for (int var = 0; var < 4; ++var){
                xCom[var] = (phiNormals[i][j][0] > 0) ? phiNormals[i][j][0]*(fluid1.u[i][j][var]-fluid1.u[i][j-1][var])/dx : phiNormals[i][j][0]*(fluid1.u[i][j+1][var]-fluid1.u[i][j][var])/dx;
                yCom[var] = (phiNormals[i][j][1] > 0) ? phiNormals[i][j][1]*(fluid1.u[i][j][var]-fluid1.u[i-1][j][var])/dy : phiNormals[i][j][1]*(fluid1.u[i+1][j][var]-fluid1.u[i][j][var])/dy;
            }
            fluid1.nDotGradPhi[i][j] = xCom + yCom;
        }
    }
    //std::cout << "ndotgradphi = " << fluid1.nDotGradPhi[10][12][0] << std::endl;

    for (int i = nGhost; i < nCells+nGhost; ++i){
        for (int j = nGhost; j < nCells+nGhost; ++j){
            std::array<double,4> xCom{{0,0,0,0}}, yCom{{0,0,0,0}};
            for (int var = 0; var < 4; ++var){
                xCom[var] = (phiNormals[i][j][0] > 0) ? phiNormals[i][j][0]*(fluid2.u[i][j+1][var]-fluid2.u[i][j][var])/dx : phiNormals[i][j][0]*(fluid2.u[i][j][var]-fluid2.u[i][j-1][var])/dx;
                yCom[var] = (phiNormals[i][j][1] > 0) ? phiNormals[i][j][1]*(fluid2.u[i+1][j][var]-fluid2.u[i][j][var])/dy : phiNormals[i][j][1]*(fluid2.u[i][j][var]-fluid2.u[i-1][j][var])/dy;
                //std::cout << fluid2.u[i][j+1][var] << " " << fluid2.u[i][j][var] << " " << fluid2.u[i][j-1][var] << std::endl;
                //std::cout << xCom[0] << " " << yCom[0] << std::endl;
            }   
            fluid2.nDotGradPhi[i][j] = xCom + yCom;
                //std::cout << "ndotgradphi (2) = " << fluid2.nDotGradPhi[i][j][0] << std::endl;

        }
    }
}



void solver::eikonalDt(){
    double aMax = 1e-15;
    for (size_t i = 0; i < phiNormals.size(); ++i) {
        for (size_t j = 0; j < phiNormals[0].size(); ++j) {
            aMax = std::max(aMax, std::sqrt(phiNormals[i][j][0]*phiNormals[i][j][0]*(phiGrads[i][j])/(dx*dx) + phiNormals[i][j][1]*phiNormals[i][j][1]*(phiGrads[i][j])/(dy*dy)));
        }
    }
    extrapDt = 0.5/aMax;
    //std::cout << "extrapDt = " << extrapDt << std::endl;
}

int sign(double value) {
    return (value >= 0) ? 1 : -1;
}

void solver::setGhostFluids(){
    if (splitFlip == 0) {
        maxIter = 101;
    } else {
        maxIter = 21;
    }

    calcPhiNormals();
    std::cout<< "phi normals calculated" << std::endl;
    for (int iter = 0; iter < maxIter; ++iter){
        setInterface();
        std::cout << "interface set" << std::endl;
        calcnDotPhiNormals();
        std::cout << "nDotPhiNormals calculated" << std::endl;
        eikonalDt();
        std::cout << "eikonal dt calculated" << std::endl;

        for (int i = nGhost; i < nCells+nGhost; ++i){
            for (int j = nGhost; j < nCells+nGhost; ++j){
                fluid1.u[i][j] = fluid1.u[i][j] - sign(phi[i][j])*extrapDt*(phi[i][j]>0)*fluid1.nDotGradPhi[i][j];
                fluid2.u[i][j] = fluid2.u[i][j] - sign(phi[i][j])*extrapDt*(phi[i][j]<=0)*fluid2.nDotGradPhi[i][j];
            }
        }
    }
    setInterface();

}

void solver::laxFriedrichs(fluid& fluid, EOS* eos){
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = LF(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost-1][j+nGhost],eos); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = LF(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost][j+nGhost-1],eos); // fluxes stored in conserved variable form
            }
        }
    }
}

void solver::richt(fluid& fluid, EOS* eos){
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = flux(eos->consvToPrim(laxFhalf(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost-1][j+nGhost],eos)),eos); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = flux(eos->consvToPrim(laxFhalf(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost][j+nGhost-1],eos)),eos); // fluxes stored in conserved variable form
            }
        }
    }
}

void solver::FORCE(fluid& fluid, EOS* eos){
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost-1][j+nGhost],eos)),eos)+LF(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost-1][j+nGhost],eos)); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost][j+nGhost-1],eos)),eos)+LF(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost][j+nGhost-1],eos)); // fluxes stored in conserved variable form
            }
        }
    }
}

void solver::SLIC(fluid& fluid, EOS* eos){

    calcHalfSlopes(fluid); 

    calcr(fluid);

    calcUBars(fluid);

    updUBars(fluid,eos); 

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1],eos)),eos)+LF(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1],eos)); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j],eos)),eos)+LF(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j],eos)); // fluxes stored in conserved variable form
            }
        }
    }
    /*
    std::cout << "u" << std::endl;
    print_arr(fluid.u,RHO);
    if (direction == YDIR){
        
        std::cout << "half slopes" << std::endl;
        print_arr(fluid.halfSlopesY,RHO);
        std::cout << "r" << std::endl;
        print_arr(fluid.rY,RHO);
        std::cout << "uBarsLY" << std::endl;
        print_arr(fluid.uBarLY,RHO);
        std::cout << "uBarsRY" << std::endl;
        print_arr(fluid.uBarRY,RHO);
        std::cout << "fluxesY" << std::endl;
        print_arr(fluid.fluxesY,RHO);
    } else {
        std::cout << "half slopes" << std::endl;
        print_arr(fluid.halfSlopesX,RHO);
        std::cout << "r" << std::endl;
        print_arr(fluid.rX,RHO);
        std::cout << "uBarsLX" << std::endl;
        print_arr(fluid.uBarLX,RHO);
        std::cout << "uBarsRX" << std::endl;
        print_arr(fluid.uBarRX,RHO);
        std::cout << "fluxesX" << std::endl;
        print_arr(fluid.fluxesX,RHO);

    }
        */
}

std::array<double,4> solver::riemannSolver(std::array<double,4> left,std::array<double,4> right, EOS* eos){
    riemann solution(eos->get_gamma(),eos->get_gamma(),eos->consvToPrim(left),eos->consvToPrim(right),direction,-0.5*dx,0.5*dx,0,dt,2,0,0);
    std::array<double,4> result = solution.exctRiemann();
    //if (sqrt(result[UX]*result[UX]) < 1e-4){result[UX] = 0.0;}
    //if (sqrt(result[UY]*result[UY]) < 1e-4){result[UY] = 0.0;}
    return result;
}

std::array<double,4> solver::HLLC(std::array<double,4> left,std::array<double,4> right, EOS* eos){ // takes conserved variables
    std::array<double,4> HLLFlux;
    double SL = 0.0, SR = 0.0;
    std::array<double,4> LPrim = eos->consvToPrim(left);
    std::array<double,4> RPrim = eos->consvToPrim(right);

    if (direction == XDIR){
        double SPlus = std::max(abs(LPrim[UX])+eos->calcSoundSpeed(LPrim),abs(RPrim[UX])+eos->calcSoundSpeed(RPrim));
        SL = -SPlus; SR = SPlus;
        //SL = LPrim[UX]-eos->calcSoundSpeed(LPrim);
        //SR = RPrim[UX]+eos->calcSoundSpeed(RPrim);

        if (0 <= SL){
            HLLFlux = flux(LPrim,eos);
        } else if (0 >= SR){
            HLLFlux = flux(RPrim,eos);
        } else {
            double Sdiff = SR - SL;
            if (std::abs(Sdiff) < 1e-8) Sdiff = (SL < 0) ? -1e-8 : 1e-8;
            HLLFlux = (1/(Sdiff))*(SR*flux(LPrim,eos)-SL*flux(RPrim,eos)+SL*SR*(right-left));
        }

    } else {
        double SPlus = std::max(abs(LPrim[UY])+eos->calcSoundSpeed(LPrim),abs(RPrim[UY])+eos->calcSoundSpeed(RPrim));
        SL = -SPlus; SR = SPlus;
        //SL = LPrim[UY]-eos->calcSoundSpeed(LPrim);
        //SR = RPrim[UY]+eos->calcSoundSpeed(RPrim);

        if (0 <= SL){
            HLLFlux = flux(LPrim,eos);
        } else if (0 >= SR){
            HLLFlux = flux(RPrim,eos);
        } else {
            double Sdiff = SR - SL;
            if (std::abs(Sdiff) < 1e-8) Sdiff = (SL < 0) ? -1e-8 : 1e-8;
            HLLFlux = (1/(Sdiff))*(SR*flux(LPrim,eos)-SL*flux(RPrim,eos)+SL*SR*(right-left));
        }
    }
    return HLLFlux;
}

void solver::MUSCL(fluid& fluid, EOS* eos){

    calcHalfSlopes(fluid); 

    calcr(fluid);

    calcUBars(fluid);

    updUBars(fluid,eos);

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = flux(riemannSolver(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1],eos),eos);
                //fluid.fluxesX[i][j] = HLLC(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1],eos);
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { // why am i getting a Y flux of vx? 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = flux(riemannSolver(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j],eos),eos);
                //fluid.fluxesY[i][j] = HLLC(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j],eos);
            }
        }
    }
}

void solver::godunov(fluid& fluid, EOS* eos){

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = flux(riemannSolver(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost-1][j+nGhost],eos),eos);
                //fluid.fluxesX[i][j] = HLLC(fluid.u[i+nGhost-1][j+1],fluid.u[i+nGhost-1][j+2],eos);
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = flux(riemannSolver(fluid.u[i+nGhost-1][j+nGhost-1],fluid.u[i+nGhost][j+nGhost-1],eos),eos);
                //fluid.fluxesY[i][j] = HLLC(fluid.u[i+1][j+nGhost-1],fluid.u[i+2][j+nGhost-1],eos);
            }
            
        }
    }
}



void solver::transmissiveBC(fluid& fluid){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<fluid.uPlus1.size(); ++j){ // sets rows (i)
            fluid.uPlus1[i][j] = fluid.u[nGhost][j];
            fluid.u[i][j] = fluid.u[nGhost][j];
            fluid.uPlus1[nCells+2*nGhost-1-i][j] = fluid.u[nCells+nGhost-1][j];
            fluid.u[nCells+2*nGhost-1-i][j] = fluid.u[nCells+nGhost-1][j];
        }
    }
    for (size_t i=0; i<fluid.uPlus1.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            fluid.uPlus1[i][j] = fluid.u[i][nGhost];
            fluid.u[i][j] = fluid.u[i][nGhost];
            fluid.uPlus1[i][nCells+2*nGhost-1-j] = fluid.u[i][nCells+nGhost-1];
            fluid.u[i][nCells+2*nGhost-1-j] = fluid.u[i][nCells+nGhost-1];
        }
    }
}

void solver::cylTransmissiveBC(fluid& fluid){
    for (int var =0; var<4;++var){
        if (var == 0 || var == 2 || var == 3){
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<fluid.uPlus1.size(); ++j){ // sets rows (ie z)
                    fluid.uPlus1[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.u[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.uPlus1[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                    fluid.u[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<fluid.uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    fluid.uPlus1[i][nGhost-1-j][var] = fluid.u[i][nGhost+j][var];
                    fluid.u[i][nGhost-1-j][var] = fluid.u[i][nGhost+j][var];
                    fluid.uPlus1[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                    fluid.u[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                }
            }
        } else {
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<fluid.uPlus1.size(); ++j){ // sets rows (ie z)
                    fluid.uPlus1[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.u[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.uPlus1[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                    fluid.u[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<fluid.uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    fluid.uPlus1[i][nGhost-1-j][var] = (-1)*fluid.u[i][nGhost+j][var];
                    fluid.u[i][nGhost-1-j][var] = (-1)*fluid.u[i][nGhost+j][var];
                    fluid.uPlus1[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                    fluid.u[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                }
            }
        }
    }
}

void solver::phiBC(){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<phi.size(); ++j){ // sets rows (i)
            phiPlus1[i][j] = phi[nGhost][j];
            phi[i][j] = phi[nGhost][j];
            phiPlus1[nCells+2*nGhost-1-i][j] = phi[nCells+nGhost-1][j];
            phi[nCells+2*nGhost-1-i][j] = phi[nCells+nGhost-1][j];
        }
    }
    for (size_t i=0; i<phi.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            phiPlus1[i][j] = phi[i][nGhost];
            phi[i][j] = phi[i][nGhost];
            phiPlus1[i][nCells+2*nGhost-1-j] = phi[i][nCells+nGhost-1];
            phi[i][nCells+2*nGhost-1-j] = phi[i][nCells+nGhost-1];
        }
    }
}

void solver::phiUpdate(int splitFlip){
    phiBC();

    double uVal;

    phiOld = phi;
    if (splitFlip%2 == 0){
        for (std::vector<double>::size_type i = nGhost-1; i < phi.size()-nGhost+1; ++i){
            for (std::vector<double>::size_type j = nGhost; j < phi[0].size() - nGhost; ++j) {
                if (phi[i][j] >= 0){
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
        //print_arr(u,RHO);
        //std::cout << "phiX updated" << std::endl;
        for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < phi[0].size() - nGhost+1; ++j) {
                if (phi[i][j] >= 0){
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
        for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < phi[0].size() - nGhost+1; ++j) {
                if (phi[i][j] >= 0){
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

        for (std::vector<double>::size_type i = nGhost-1; i < phi.size()-nGhost+1; ++i){
            for (std::vector<double>::size_type j = nGhost; j < phi[0].size() - nGhost; ++j) {
                if (phi[i][j] >= 0){
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
}

void solver::reinitPhi(){
    for (int iter = 0; iter < 11; iter++){
        //std::cout << "iteration " << iter << std::endl;
        for (int i = nGhost; i<nCells+nGhost; ++i){
            for (int j = nGhost; j < nCells+nGhost; ++j){
                phiPlus1[i][j] = phi[i][j] - extrapDt*sign(phi[i][j])*(phiGrads[i][j]-1.0);
            }
        }
    }
    phi = phiPlus1;
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

void solver::fixFreshlyCleared(){
    for (std::array<int,2> cell : freshlyCleared){
        neighbourAvg(fluid1,cell[0],cell[1]);
        neighbourAvg(fluid2,cell[0],cell[1]);
    }
}

void solver::pointsUpdate(fluid& fluid){ //second index = x, first index = y ; phiSign = 0 for negative phi, 1 for positive phi
    //cylTransmissiveBC();
    if (direction == XDIR){
        for (std::vector<double>::size_type i = nGhost-1 ; i < fluid.u.size()-nGhost+1; ++i){
            for (std::vector<double>::size_type j = nGhost; j < fluid.u[0].size() - nGhost; ++j) {
                fluid.uPlus1[i][j] = fluid.u[i][j] - (dt / dx) * (fluid.fluxesX[i-nGhost+1][j-nGhost+1] - fluid.fluxesX[i-nGhost+1][j-nGhost]); // flux[i + 1] and flux[i] for the update
            }
        }
        fluid.u = fluid.uPlus1;
        //print_arr(u,RHO);
    } else if (direction == YDIR){
        for (std::vector<double>::size_type i = nGhost; i < fluid.u.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < fluid.u[0].size() - nGhost+1; ++j) {
                fluid.uPlus1[i][j] = fluid.u[i][j] - (dt / dy) * (fluid.fluxesY[i-nGhost+1][j-nGhost+1] - fluid.fluxesY[i-nGhost][j-nGhost+1]); // flux[i + 1] and flux[i] for the update
            }
        }
        fluid.u = fluid.uPlus1;
        //print_arr(u,RHO);
    }
    
};

// -------------------- Flux functions --------------------- //

std::array<double,4> solver::finvectLF(const std::array<double,4> x)const{
    double a=1;
    return a*x;
};

std::array<double,4> solver::fBurgersLF(const std::array<double,4> x)const{
    return 0.5*(x * x);
};

std::array<double,4> solver::fEuler(std::array<double,4> arr, EOS* eos){ // takes primitive variables
    std::array<double,4> result;
    if (direction == XDIR){
        result[RHO] = arr[RHO]*arr[UX]; // rho*v
        result[XMOM] = arr[RHO]*arr[UX]*arr[UX]+arr[PRES]; // rho*v^2 + p
        result[YMOM] = arr[RHO]*arr[UX]*arr[UY]; // rho*v^2 + p
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UX]; // (E + p)v
    } else if (direction == YDIR){
        result[RHO] = arr[RHO]*arr[UY]; // rho*v
        result[XMOM] = arr[RHO]*arr[UX]*arr[UY]; // rho*v^2 + p
        result[YMOM] = arr[RHO]*arr[UY]*arr[UY]+arr[PRES]; // rho*v^2 + p
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UY]; // (E + p)v
    }

    return result;
}

// ------------------------- secondary fluxes --------------------- //

std::array<double,4> solver::LF(std::array<double,4> v1, std::array<double,4> v2, EOS* eos)const{
    if (direction == XDIR){
        return 0.5 * (flux(eos->consvToPrim(v1),eos) + flux(eos->consvToPrim(v2),eos)) + 0.5 * (dx / dt) * (v1 - v2);
    } else {
        return 0.5 * (flux(eos->consvToPrim(v1),eos) + flux(eos->consvToPrim(v2),eos)) + 0.5 * (dy / dt) * (v1 - v2);
    }
};

std::array<double,4> solver::laxFhalf(std::array<double,4> ui, std::array<double,4> uip1, EOS* eos)const{
    if (direction == XDIR){
        return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos->consvToPrim(uip1),eos)-flux(eos->consvToPrim(ui),eos)));
    } else if (direction == YDIR){
        return 0.5*(ui+uip1)-(0.5*(dt/dy)*(flux(eos->consvToPrim(uip1),eos)-flux(eos->consvToPrim(ui),eos)));
    } else {
        throw std::runtime_error("laxfhalf doesnt know which dimension");
    }
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(fluid& fluid){ // only using y-gradient - ie correspondng to y-faces; need extra column for x direction //DONE
    //std::cout << "calculating halfslopes n printing u" << std::endl;
    ////print_arr(u,RHO);
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.halfSlopesX.size(); ++i){
            for (size_t j = 0; j < fluid.halfSlopesX[0].size(); ++j){
                fluid.halfSlopesX[i][j] = fluid.u[i+nGhost-1][j+nGhost-1]-fluid.u[i+nGhost-1][j+nGhost-2];
                // i + nGhost - 1 is the second-outermost row/column
            }
        }
    } else {
        for (size_t i = 0; i < fluid.halfSlopesY.size(); ++i){
            for (size_t j = 0; j < fluid.halfSlopesY[0].size(); ++j){
                fluid.halfSlopesY[i][j] = fluid.u[i+nGhost-1][j+nGhost-1]-fluid.u[i+nGhost-2][j+nGhost-1];
            }
        }
    };
};

void solver::calcr(fluid& fluid){ // DONE
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.rX.size(); ++i){
            for (size_t j = 0; j < fluid.rX[0].size(); ++j){
                fluid.rX[i][j] = elementDivide(fluid.halfSlopesX[i][j],fluid.halfSlopesX[i][j+1]);
                //fluid.rX[i][j] = {0,0,0,0};
            }
        }
    } else {
        for (size_t i = 0; i < fluid.rY.size(); ++i){
            for (size_t j = 0; j < fluid.rY[0].size(); ++j){
                fluid.rY[i][j] = elementDivide(fluid.halfSlopesY[i][j],fluid.halfSlopesY[i+1][j]);
                //fluid.rY[i][j] = {0,0,0,0};
            }
        }
    };
    ////print_arr(r);
    //std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(fluid& fluid){ // calculates for all u except leftmost cell
    if (direction == XDIR){
        for (size_t i = 0; i<fluid.uBarLX.size(); ++i){
            for (size_t j = 0; j<fluid.uBarLX[0].size(); ++j){
            
                //print_vect(slopeLim(r[i-1]));
                //std::cout << std::endl << "slope ";
                //print_vect(calcSlope(slopeWeight, halfSlopes[i-1], halfSlopes[i]));
                //std::cout << std::endl;
                fluid.uBarLX[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] - 0.5 * ( slopeLim(fluid.rX[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesX[i][j], fluid.halfSlopesX[i][j+1]) );
                fluid.uBarRX[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] + 0.5 * ( slopeLim(fluid.rX[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesX[i][j], fluid.halfSlopesX[i][j+1]) );
                //std::cout << slopeLim(rX[i][j])[RHO] << " ";
            }
            //std::cout << std::endl;
        }
    } else {
        for (size_t i = 0; i<(fluid.uBarLY.size()); ++i){
            for (size_t j = 0; j<(fluid.uBarLY[0].size()); ++j){
            
                //print_vect(slopeLim(r[i-1]));
                //std::cout << std::endl << "slope ";
                //print_vect(calcSlope(slopeWeight, halfSlopes[i-1], halfSlopes[i]));
                //std::cout << std::endl;
                fluid.uBarLY[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] - 0.5 * ( slopeLim(fluid.rY[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesY[i][j], fluid.halfSlopesY[i+1][j]) );
                fluid.uBarRY[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] + 0.5 * ( slopeLim(fluid.rY[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesY[i][j], fluid.halfSlopesY[i+1][j]) );
                //std::cout << slopeLim(rY[i][j])[RHO] << " ";
            }
            //std::cout << std::endl;
        }
    };
};


void solver::updUBars(fluid& fluid, EOS* eos){
    if (direction == XDIR){
        for (size_t i = 0; i<(fluid.uBarLupdX.size()); ++i){
            for (size_t j = 0; j<(fluid.uBarLupdX[0].size()); ++j){
                fluid.uBarLupdX[i][j] = fluid.uBarLX[i][j]-0.5*(dt/dx)*(flux(eos->consvToPrim(fluid.uBarRX[i][j]),eos)-flux(eos->consvToPrim(fluid.uBarLX[i][j]),eos));
                fluid.uBarRupdX[i][j] = fluid.uBarRX[i][j]-0.5*(dt/dx)*(flux(eos->consvToPrim(fluid.uBarRX[i][j]),eos)-flux(eos->consvToPrim(fluid.uBarLX[i][j]),eos));
            }
        }
        fluid.uBarLX = fluid.uBarLupdX;
        fluid.uBarRX = fluid.uBarRupdX;
    } else if (direction == YDIR){
        for (size_t i = 0; i<(fluid.uBarLupdY.size()); ++i){
            for (size_t j = 0; j<(fluid.uBarLupdY[0].size()); ++j){
                fluid.uBarLupdY[i][j] = fluid.uBarLY[i][j]-0.5*(dt/dy)*(flux(eos->consvToPrim(fluid.uBarRY[i][j]),eos)-flux(eos->consvToPrim(fluid.uBarLY[i][j]),eos));
                fluid.uBarRupdY[i][j] = fluid.uBarRY[i][j]-0.5*(dt/dy)*(flux(eos->consvToPrim(fluid.uBarRY[i][j]),eos)-flux(eos->consvToPrim(fluid.uBarLY[i][j]),eos));
            }
        }
        fluid.uBarLY = fluid.uBarLupdY;
        fluid.uBarRY = fluid.uBarRupdY;
    }
};

std::array<double,4> solver::calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,4> solver::minbee(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
    double minb;
    for (int k = 0; k<4; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY || slp[k] == -INFINITY){
            minb = 0;
        } else if (slp[k] > 1){
            minb = std::min(1.0,2.0/(1.0+slp[k]));
        } else {
            minb = slp[k];
        }
        minbArr[k] = minb;
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[0],minbArr[0],minbArr[0],minbArr[0]};
    return {minBval,minBval,minBval,minBval};
    //return minbArr;
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

std::array<double,4> solver::vanAlbada(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
    double minb;
    for (int k = 0; k<4; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY || slp[k] == -INFINITY){
            minb = 0;
        } else {
            minb = std::min((slp[k]*(1.0+slp[k])/(1.0+slp[k]*slp[k])),2.0/(1.0+slp[k]));
        }
        minbArr[k] = minb;
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[0],minbArr[0],minbArr[0],minbArr[0]};
    return {minBval,minBval,minBval,minBval};
    //return minbArr;
};



// ---------------------- Misc ------------------------------ //

void solver::printPrec(double value, int precision){
    std::cout << std::fixed << std::setprecision(precision) << value << " ";
}

void solver::printScalarField(const std::vector<std::vector<double>>& field) {
    for (const auto& row : field) {
        for (const auto& value : row) {
            std::cout << std::fixed << std::setprecision(3) << value << " ";
        }
        std::cout << std::endl;
    }
}

void solver::print_arr(std::vector < std::vector< std::array<double,4> > > arr, int var) { 
    std::string row;
    for (size_t j = 0; j < arr.size(); ++j) {
        for (size_t i = 0; i < arr[j].size(); ++i) {
            row += " ";

            // Use std::ostringstream to format the number to 6 decimal places
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3) << arr[j][i][var];

            // Append the formatted string to the row
            row += oss.str();
        }
        std::cout << row << std::endl;
        row.clear(); // Clear the row for the next line
    }
};

void solver::print_vect(std::vector < std::array<double,4> > v){ // just prints a slice (0 by default)
    for (int j = 0; j<4; ++j){
        std::cout << v[0][j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::print_state(std::array<double,4> v){ // just prints a slice (0 by default)
    for (int j = 0; j<4; ++j){
        std::cout << v[j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::resize2D(int nx, int ny, std::vector< std::vector< double > >& arr){

    arr.resize(ny);
    for (size_t i = 0; i < arr.size(); ++i){
        arr[i].resize(nx);
    }

}

void solver::resize2Db(int nx, int ny, std::vector< std::vector< std::array <double,2> > >& arr){

    arr.resize(ny);
    for (size_t i = 0; i < arr.size(); ++i){
        arr[i].resize(nx);
    }

}

void solver::resize2Dc(int nx, int ny, std::vector< std::vector< std::array <double,4> > >& arr){

    arr.resize(ny);
    for (size_t i = 0; i < arr.size(); ++i){
        arr[i].resize(nx);
    }

}

std::array<double,4> solver::bilinear(double x, double y, fluid& fluid){
    int rows = nCells;
    int cols = nCells;

    // Compute indices of grid cell containing the point
    int j = ((x - x0) / dx) + nGhost - 0.5;
    int i = ((y - y0) / dy) + nGhost - 0.5;

    //std::cout << "i = " << i << " j = " << j << std::endl;

    //if (i < 0 || i > rows - 1 || j < 0 || j > cols - 1) {
        //std::cout << i << " " << j << std::endl;
        //return {NAN,NAN,NAN,NAN};
        //throw std::runtime_error("interpolation out of bounds");
    //};

    // Physical coordinates of the four surrounding points
    double xL = x0 + (j - nGhost + 0.5) * dx;
    double xR = x0 + (j + 1 - nGhost + 0.5) * dx;
    double yL = y0 + (i - nGhost + 0.5) * dy;
    double yR = y0 + (i + 1 - nGhost + 0.5) * dy;

    //std::cout << "xL = " << xL << " xR = " << xR << " yL = " << yL << " yR = " << yR << std::endl;

    // Grid values at the four points
    std::array<double,4> f00 = fluid.u[i][j];
    std::array<double,4> f10 = fluid.u[i + 1][j];
    std::array<double,4> f01 = fluid.u[i][j + 1];
    std::array<double,4> f11 = fluid.u[i + 1][j + 1];

    //std::cout << "neighbour values" << f00[0] << " " << f10[0] << " " << f01[0] << " " << f11[0] << std::endl;

    // Interpolation along x at y0 and y1
    std::array<double,4> f_x_y0 = (xR - x) / (xR - xL) * f00 + (x - xL) / (xR - xL) * f10;
    std::array<double,4> f_x_y1 = (xR - x) / (xR - xL) * f01 + (x - xL) / (xR - xL) * f11;

    // Interpolation along y
    return (yR - y) / (yR - yL) * f_x_y0 + (y - yL) / (yR - yL) * f_x_y1;
}


// ---------------------- Time Stepping ---------------------- //

void solver::setDt(){
    double aMax = -1e14;
    //std::cout << "setting dt" << std::endl;
    for (std::vector<double>::size_type i = 0; i < fluid1.u.size(); ++i) {
        for (std::vector<double>::size_type j = 0; j < fluid1.u.size(); ++j) {
            aMax = std::max(aMax, std::max(std::sqrt(pow(eos[1]->consvToPrim(fluid2.u[i][j])[UX],2)+pow(eos[1]->consvToPrim(fluid2.u[i][j])[UY],2))+eos[1]->calcSoundSpeed(eos[1]->consvToPrim(fluid2.u[i][j])), std::sqrt(pow(eos[0]->consvToPrim(fluid1.u[i][j])[UX],2)+pow(eos[0]->consvToPrim(fluid1.u[i][j])[UY],2))+eos[0]->calcSoundSpeed(eos[0]->consvToPrim(fluid1.u[i][j]))));

        }
    }
    //if (time > 0.16){writeInterval = 0.001;}
    dt = dtReducer * cour * std::min(dx,dy) / std::fabs(aMax);

    double writeTime = writeInterval*(std::floor((time+1e-9) / writeInterval)+1); //returing 0.29?

    if (endTime - time <= dt && writeTime-time >= dt){
        //std::cout << "option 0 " << std::endl;
        dt = endTime - time;
        assert(dt != 0);
        time = endTime;
        std::cout << "last step:" << dt << std::endl;
    } else if (writeTime-time <= dt && writeTime-time > 0){ // manually reducing dt to hit the write times
        //std::cout << "option 1 " << std::endl;

        dt = writeTime - time;
        assert(dt != 0);
        time = writeTime;
    } else {
        //std::cout << "option 2 " << std::endl;
        time = time + dt;
    }

    if (std::fabs(writeTime - time) < 1e-10){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    std::cout << "aMax = " << aMax << ", dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double y_0, double y_1, double t0, double t1, int n_Cells, int n_Ghosts, double c, double g)
    : fluid1(n_Cells, n_Ghosts), fluid2(n_Cells,n_Ghosts), gamma(g) {
    assert(t1>t0);
    assert(g > 0);

    x0 = x_0; x1 = x_1;
    y0 = y_0; y1 = y_1;
    startTime = t0; endTime = t1;
    nCells = n_Cells;
    nGhost = n_Ghosts;

    resize2D(nCells+2*nGhost,nCells+2*nGhost,phi);
    resize2D(nCells+2*nGhost,nCells+2*nGhost,phiPlus1);

    dx = (x1 - x0)/nCells;
    dy = (y1 - y0)/nCells;


    slopeWeight = 1;

    cour = c;

    variables = {"rho","vx","vy","p"};

    //varMap["rho"] = 0; varMap["v"] = 1; varMap["p"] = 2;
    for (size_t i = 0; i<variables.size(); ++i){
        varMap[variables[i]] = i;
    }
    
};

// ---------------------- utilities -------------------------- //

int solver::ghosts(){
    return nGhost;
};

double solver::get_dx()const{
    return dx;
};

double solver::get_dy()const{
    return dy;
};

void solver::init(std::vector< std::vector< std::array<double,4> > > init, fluid& fluid) {
    assert(init.size() == fluid.u.size());
    fluid.u = init;
};

void solver::phiInit(std::vector< std::vector < double > > init) {
    assert(init.size() == phi.size());
    phi = init;
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

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData(std::string varName) const{ // modify to write both fluids
        // making data file
        std::string filename = dirName + "/" + varName + "/" + std::to_string(time).substr(0,5);
        std::ofstream writeFile(filename);
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (std::vector<double>::size_type i=0; i<(fluid1.u.size())-0;i++){ // shoudl be from nGhost to size-nGhost
                double y = y0 + (i-static_cast<double>(nGhost)+0.5)*dy;
                for (std::vector<double>::size_type j=0; j<(fluid1.u[0].size())-0;j++){
                    double x = x0 + (j-static_cast<double>(nGhost)+0.5)*dx;
                    writeFile << x << " " << y << " " << eos[0]->consvToPrim(fluid1.u[i][j])[varMap.at(varName)] << " " << eos[1]->consvToPrim(fluid2.u[i][j])[varMap.at(varName)] << " " << phi[i][j] << std::endl;
                    //std::cout << u[i] << " ";
                }
                writeFile << "\n";
            }
            writeFile.close();
            //std::cout << filename << "written successfully" << std::endl;
        }else{
            std::cout << "failed to write." << std::endl;
        }

}

std::array<double,4> solver::set_vals(double rho, double vx, double vy, double pres){
        std::array<double,4> result = {rho,vx,vy,pres};
        return result;
};


// ------------------- Overloading operators ------------------ //

std::array<double, 4> operator-(const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = lhs[i] - rhs[i];
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << "RHS:" << rhs[0] << " " << rhs[1] << " " << rhs[2] << " " << rhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("element subtract error");
        }
    }
    return result;
}

std::array<double, 4> operator+(const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = lhs[i] + rhs[i];
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << "RHS:" << rhs[0] << " " << rhs[1] << " " << rhs[2] << " " << rhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("element add error");
        }
    }
    return result;
}

std::array<double, 4> operator*(const double scalar,const std::array<double, 4>& lhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = scalar*lhs[i];
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("scalar multiply error");
        }
    }
    return result;
}

std::array<double, 4> operator*(const std::array<double, 4>& rhs,const std::array<double, 4>& lhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = rhs[i]*lhs[i];
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << "RHS:" << rhs[0] << " " << rhs[1] << " " << rhs[2] << " " << rhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("element multiply error");
        }
    }
    return result;
}

std::array<double, 4> elementDivide(const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        if (rhs[i] == lhs[i]){
            result[i] = 1;
        }
        else if (rhs[i] == 0){
            result[i] = INFINITY;
        }
        else {
            result[i] = lhs[i] / rhs[i];
        }
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << "RHS:" << rhs[0] << " " << rhs[1] << " " << rhs[2] << " " << rhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("element divide error");
        }


    }
    return result;
}

std::array<double, 4> operator/(const std::array<double, 4>& lhs, const double scalar) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = lhs[i]/scalar;
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("scalar divide error");
        }
    }
    return result;
}



