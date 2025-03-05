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
    std::cout << "set BCs" << std::endl;

    //double realTime = 0.0;
    //print_arr(fluid2.u,UX);
    //print_arr(fluid2.u,UY);
    //print_arr(fluid2.u,PRES);

    splitFlip = 0;
    do{
        
        if (splitFlip < (endTime-startTime)/20){
            dtReducer = 0.5;
        } else {
            dtReducer = 1.0;
        }
        
        //time = 0+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        /*
        phiBC();
        std::cout << "finding boundary" << std::endl;
        this->findBoundary();
        //std::cout << "number of interface cells = " << interfaceCells.size() << std::endl;
        std::cout << "boundary found" << std::endl;

        //printInterfaceArray(phi, interfaceCells, "interface.txt");

        this->calcInterfaceNormals(fluid1);
        this->calcInterfaceNormals(fluid2);
        std::cout << "phi grads calculated" << std::endl;
        this->calcInterface(fluid1);
        this->calcInterface(fluid2);
        std::cout << "interface calculated" << std::endl;
        this->interpInterfaceStates(fluid1);
        this->interpInterfaceStates(fluid2);
        std::cout << "interface states interpolated" << std::endl;
        this->resolveVelocities(fluid1);
        this->resolveVelocities(fluid2);
        std::cout << "resolved velocities" << std::endl;
        this->interfaceRiem(fluid1);
        this->interfaceRiem(fluid2);
        std::cout << "interface riem" << std::endl;
        this->calcStarStates(fluid1);
        this->calcStarStates(fluid2);
        std::cout << "star states calcd" << std::endl;

        //print_arr(fluid1.u,UY);
        //std::cout << "u2 pre-ghost evolution" << std::endl;
        //print_arr(fluid2.u,UY);

        this->setGhostFluids();
        //std::cout << "u1 post-ghost evolution" << std::endl;
        //print_arr(fluid1.u,UY);
        //std::cout << "u2 post-ghost evolution" << std::endl;
        //print_arr(fluid2.u,UY);
        //std::cout << "set ghost fluids" << std::endl; 

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);
        */
        this->setDt();
        double origDt = dt;
        /*
        //realTime = time;

        std::cout << "time is " << time << std::endl;

        //std::cout << "setting phiBC" << std::endl;
        phiBC();
        //std::cout << "updating phi" << std::endl;

        //time = 3+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        
        phiUpdate(splitFlip);
        std::cout << "phi updated" << std::endl;

        fluid1.interfacePhis.resize(fluid1.interfaceCells.size());    
        fluid2.interfacePhis.resize(fluid2.interfaceCells.size());
        resize2Dc(nCellsX+2*nGhost,nCellsY+2*nGhost,fluid1.nDotGradPhi);
        resize2Dc(nCellsX+2*nGhost,nCellsY+2*nGhost,fluid2.nDotGradPhi);
    
        for (size_t i = 0; i<fluid1.interfaceCells.size(); ++i){
            fluid1.interfacePhis[i] = phi[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]];
        }
        for (size_t i = 0; i<fluid2.interfaceCells.size(); ++i){
            fluid2.interfacePhis[i] = phi[fluid2.interfaceCells[i][0]][fluid2.interfaceCells[i][1]];
        }

        //time = 4+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        //fixFreshlyCleared();
        //std::cout << "freshly cleared fixed" << std::endl;

        //time = 5+(splitFlip*7);
        //writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");

        //time = 0.0;
        
        //if (splitFlip % 10 == 0){
            //phiBC();
            //reinitPhi();
            //phiBC();
            //std::cout << "phi reinitialised" << std::endl;
        //}
        */
        splitFlip++;
        std::cout << "Flux calc & point updates" << std::endl;
        if (splitFlip%2 == 0){
            //std::cout << "x first" << std::endl;
            direction = XDIR;
            dt = origDt/2.0;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
            //std::cout << "u1" << std::endl;
            //print_arr(fluid1.u,UY);
            direction = YDIR;
            dt = origDt;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
            
            direction = XDIR;
            dt = origDt/2.0;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);

            //std::cout << "u1" << std::endl;
            //print_arr(fluid1.u,UY);
        } else {
            //std::cout << "x first" << std::endl;
            direction = YDIR;
            dt = origDt/2.0;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
            //std::cout << "u1" << std::endl;
            //print_arr(fluid1.u,UY);
            direction = XDIR;
            dt = origDt;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
            
            direction = YDIR;
            dt = origDt/2.0;
            this->fluxMethod(fluid1,eos[0]);
            this->pointsUpdate(fluid1);
            this->fluxMethod(fluid2,eos[1]);
            this->pointsUpdate(fluid2);
            transmissiveBC(fluid1);
            transmissiveBC(fluid2);
        }
        //std::cout << "points updated" << std::endl;

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);
        //setInterface();

        this->fluid1.interfaceCells={};
        this->fluid2.interfaceCells={};
        this->freshlyCleared={};

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);

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
        for (size_t j = nGhost; j < phi[0].size()-nGhost; ++j){
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
    int ny = phi.size()-nGhost-1;
    int nx = phi.size()-nGhost-1;
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
        resize2Dc(nCellsX+2*nGhost,nCellsY+2*nGhost,f.uExtrap);
    }
    resize2Db(nCellsX+2*nGhost,nCellsY+2*nGhost,phiNormals);
    resize2D(nCellsX+2*nGhost,nCellsY+2*nGhost,phiGrads);


    

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
    //printInterfaceArray(phi, fluid1.interfaceCells, "interface1.txt");
    //printInterfaceArray(phi, fluid2.interfaceCells, "interface2.txt");

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

void solver::calcInterfaceNormals(fluid& f){
    //std::cout << std::endl;
    //printScalarField(phi);
    
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        int x = f.interfaceCells[i][1];
        int y = f.interfaceCells[i][0];
        double phi_y;
        double phi_x;
        if (x == nGhost){
            phi_x = (phi[y][x+1]-phi[y][x])/dx;
        } else if (x == nCellsX+nGhost-1){
            phi_x = (phi[y][x]-phi[y][x-1])/dx;
        } else {
            phi_x = (phi[y][x+1]-phi[y][x-1])/(2*dx);
        }
        if (y == nGhost){
            phi_y = (phi[y+1][x]-phi[y][x])/dy;
        } else if (y == nCellsY+nGhost-1){
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
    for (size_t i = 0; i < f.interfaceCells.size(); ++i){
        double x_i = x0 + (f.interfaceCells[i][1]-1.5)*dx; // -2 for ghost cells, +0.5 for cell centre
        double y_i = y0 + (f.interfaceCells[i][0]-1.5)*dy;

        double x = x_i - phi[f.interfaceCells[i][0]][f.interfaceCells[i][1]]*f.interfaceNormals[i][0];
        double y = y_i - phi[f.interfaceCells[i][0]][f.interfaceCells[i][1]]*f.interfaceNormals[i][1];

        if (&f == &fluid1){
            std::array<double,2> lPos = {x-1.5*dx*f.interfaceNormals[i][0],y-1.5*dy*f.interfaceNormals[i][1]};
            std::array<double,2> rPos = {x+1.5*dx*f.interfaceNormals[i][0],y+1.5*dy*f.interfaceNormals[i][1]};
            f.interfacePositions[i] = {lPos,rPos};
        } else {
            std::array<double,2> lPos = {x-1.5*dx*f.interfaceNormals[i][0],y-1.5*dy*f.interfaceNormals[i][1]};
            std::array<double,2> rPos = {x+1.5*dx*f.interfaceNormals[i][0],y+1.5*dy*f.interfaceNormals[i][1]};
            f.interfacePositions[i] = {lPos,rPos};
        }
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
    for (size_t i = 0; i < f.interfacePositions.size(); ++i){
        std::array<double,4> lState = bilinear(f.interfacePositions[i][0][0],f.interfacePositions[i][0][1],fluid1, f.interfaceNormals[i]);
        std::array<double,4> rState = bilinear(f.interfacePositions[i][1][0],f.interfacePositions[i][1][1],fluid2, f.interfaceNormals[i]);
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

void solver::resolveVelocities(fluid& f){
    /*
    resolves normal and tangential components of velocities of the interpolated states
    */
    for (size_t i=0; i<f.interfaceStates.size(); ++i){
        double vNormL = f.interfaceStates[i][0][UX]*f.interfaceNormals[i][0]+f.interfaceStates[i][0][UY]*f.interfaceNormals[i][1];
        double vNormR = f.interfaceStates[i][1][UX]*f.interfaceNormals[i][0]+f.interfaceStates[i][1][UY]*f.interfaceNormals[i][1];
        double vTanL = sqrt((f.interfaceStates[i][0][UX]*f.interfaceStates[i][0][UX]+f.interfaceStates[i][0][UY]*f.interfaceStates[i][0][UY])-vNormL*vNormL);
        double vTanR = sqrt((f.interfaceStates[i][1][UX]*f.interfaceStates[i][1][UX]+f.interfaceStates[i][1][UY]*f.interfaceStates[i][1][UY])-vNormR*vNormR);
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

void solver::interfaceRiem(fluid& f){

    for (size_t i=0; i<f.interfaceStates.size(); ++i){
        std::array<double,4> lState = {f.interfaceStates[i][0][0],f.resolvedVelocities[i][0][0],f.resolvedVelocities[i][0][1],f.interfaceStates[i][0][3]};  
        std::array<double,4> rState = {f.interfaceStates[i][1][0],f.resolvedVelocities[i][1][0],f.resolvedVelocities[i][1][1],f.interfaceStates[i][1][3]};
        //std::cout << "riem" << std::endl;
        riemann solution(eos[0]->get_gamma(),eos[1]->get_gamma(),lState,rState,1,-1.5*(dx+dy)/2.0,1.5*(dx+dy)/2.0,0,dt,2,0,0); // need p_inf setup
        //std::cout << "riemann solver object initialised" << std::endl;
        std::array<double,4> result = solution.exctRiemann();
        if (result[RHO] <= 0){
            print_state(lState); print_state(rState);
            throw std::runtime_error("interface riemann solution rho < 0");
        }
        //if (sqrt(result[UX]*result[UX]) < 1e-4){result[UX] = 0.0;}
        //if (sqrt(result[UY]*result[UY]) < 1e-4){result[UY] = 0.0;}
        f.riemInterfaceStates[i] = result;
        //print_state(result);
        //solution.exctRiemann();
    }
}

void solver::calcStarStates(fluid& f){
    for (size_t i=0; i<f.interfaceStates.size(); ++i){
        if (&f == &fluid1){
            std::array<double,4> sL = f.interfaceStates[i][0];

            //double utLx = sL[UX]-(f.resolvedVelocities[i][0][0])*f.interfaceNormals[i][0];
            //double utLy = sL[UY]-(f.resolvedVelocities[i][0][0])*f.interfaceNormals[i][1];
            double uLStarx = (f.riemInterfaceStates[i][1]-(f.resolvedVelocities[i][0][0]))*f.interfaceNormals[i][0]+f.interfaceStates[i][0][UX];
            double uLStary = (f.riemInterfaceStates[i][1]-(f.resolvedVelocities[i][0][0]))*f.interfaceNormals[i][1]+f.interfaceStates[i][0][UY];
            //double uLStarx = f.riemInterfaceStates[i][1]*f.interfaceNormals[i][0]+utLx;
            //double uLStary = f.riemInterfaceStates[i][1]*f.interfaceNormals[i][1]+utLy;
            
            std::array<double,4> uStarL = {f.riemInterfaceStates[i][0],uLStarx,uLStary,f.riemInterfaceStates[i][3]};
            f.starStates[i] = uStarL;

            //std::cout << "star states left: " << uStarL[0] << " " << uStarL[1] << " " << uStarL[2] << " " << uStarL[3] << std::endl;

        
        } else if (&f == &fluid2){
            std::array<double,4> sR = f.interfaceStates[i][1]; // primitive variables

            //double utRx = sR[UX]-(f.resolvedVelocities[i][1][0])*f.interfaceNormals[i][0];
            //double utRy = sR[UY]-(f.resolvedVelocities[i][1][1])*f.interfaceNormals[i][1];
            //double uRStarx = f.riemInterfaceStates[i][1]*f.interfaceNormals[i][0]+utRx;
            //double uRStary = f.riemInterfaceStates[i][1]*f.interfaceNormals[i][1]+utRy;
            double uRStarx = (f.riemInterfaceStates[i][1]-(f.resolvedVelocities[i][1][0]))*f.interfaceNormals[i][0]+f.interfaceStates[i][1][UX];
            double uRStary = (f.riemInterfaceStates[i][1]-(f.resolvedVelocities[i][1][0]))*f.interfaceNormals[i][1]+f.interfaceStates[i][1][UY];

            std::array<double,4> uStarR = {f.riemInterfaceStates[i][0],uRStarx,uRStary,f.riemInterfaceStates[i][3]};
            f.starStates[i] = uStarR;
            
            //std::cout << "star states right: " << uStarR[0] << " " << uStarR[1] << " " << uStarR[2] << " " << uStarR[3] << std::endl;

        }

    }
}

void solver::calcPhiGrad(){
    phiBC();

    for (int y = nGhost; y < nCellsY+nGhost; ++y){
        for (int x = nGhost; x < nCellsX+nGhost; ++x){
            double phi_y;
            double phi_x;
            if (x == nGhost){
                phi_x = (phi[y][x+1]-phi[y][x])/dx;
            } else if (x == nCellsX+nGhost-1){
                phi_x = (phi[y][x]-phi[y][x-1])/dx;
            } else {
                phi_x = (phi[y][x+1]-phi[y][x-1])/(2*dx);
            }
            if (y == nGhost){
                phi_y = (phi[y+1][x]-phi[y][x])/dy;
            } else if (y == nCellsY+nGhost-1){
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
    transmissiveBC(fluid1);
    transmissiveBC(fluid2);
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

void solver::calcnDotPhiNormals(){
    //vectorTransBC(fluid1.u);
    //vectorTransBC(fluid2.u);


    for (int i = nGhost; i < nCellsY+nGhost; ++i){
        for (int j = nGhost; j < nCellsX+nGhost; ++j){
            std::array<double,4> xCom, yCom;
            
            xCom = (phiNormals[i][j][0] > 0) ? phiNormals[i][j][0]*(fluid1.u[i][j]-fluid1.u[i][j-1])/dx : phiNormals[i][j][0]*(fluid1.u[i][j+1]-fluid1.u[i][j])/dx;
            yCom = (phiNormals[i][j][1] > 0) ? phiNormals[i][j][1]*(fluid1.u[i][j]-fluid1.u[i-1][j])/dy : phiNormals[i][j][1]*(fluid1.u[i+1][j]-fluid1.u[i][j])/dy;

            fluid1.nDotGradPhi[i][j] = xCom + yCom;
        }
    }
    //vectorTransBC(fluid1.nDotGradPhi);

    //std::cout << "ndotgradphi = " << fluid1.nDotGradPhi[10][12][0] << std::endl;

    for (int i = nGhost; i < nCellsY+nGhost; ++i){
        for (int j = nGhost; j < nCellsX+nGhost; ++j){
            std::array<double,4> xCom, yCom;

            xCom = (phiNormals[i][j][0] > 0) ? phiNormals[i][j][0]*(fluid2.u[i][j+1]-fluid2.u[i][j])/dx : phiNormals[i][j][0]*(fluid2.u[i][j]-fluid2.u[i][j-1])/dx;
            yCom = (phiNormals[i][j][1] > 0) ? phiNormals[i][j][1]*(fluid2.u[i+1][j]-fluid2.u[i][j])/dy : phiNormals[i][j][1]*(fluid2.u[i][j]-fluid2.u[i-1][j])/dy;

            fluid2.nDotGradPhi[i][j] = xCom + yCom;
                //std::cout << "ndotgradphi (2) = " << fluid2.nDotGradPhi[i][j][0] << std::endl;

        }
    }
    //vectorTransBC(fluid2.nDotGradPhi);
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

void solver::dtExtrap(){
    double aMax = 1e-15;
    for (size_t i = 0; i < phiNormals.size(); ++i) {
        for (size_t j = 0; j < phiNormals[0].size(); ++j) {
            //print_state(fluid1.u[i][j]);
            //print_state(fluid2.u[i][j]);
            aMax = std::max(aMax, std::max(std::abs(phiNormals[i][j][0]*eos[0]->consvToPrim(fluid1.u[i][j])[UX] + phiNormals[i][j][1]*eos[0]->consvToPrim(fluid1.u[i][j])[UY]),std::abs(phiNormals[i][j][0]*eos[1]->consvToPrim(fluid2.u[i][j])[UX] + phiNormals[i][j][1]*eos[1]->consvToPrim(fluid2.u[i][j])[UY])));
        }
    }
    extrapDt = cour*std::min(dx,dy)/aMax;
    //std::cout << "extrapDt = " << extrapDt << std::endl;
}

int sign(double value) {
    return (value <= 0) ? -1 : 1;
}


void solver::fastSweeping(fluid& f){
    int f1 = (&f == &fluid1) ? 1 : 0;
    for (std::array<int,2> bPoint : f.interfaceCells){
        int i = bPoint[0];
        int j = bPoint[1];
        double fPhix = phiNormals[i][j][0];
        double fPhiy = phiNormals[i][j][1];
        if (fPhix <= 0){
            // for f1 sweep in +ve direction if fPhi < 0 and vice versa
            // for f2 sweep in -ve direction if fPhi < 0 and vice versa
            // if i,j+2 is not adjacent (else break):
            // start fast sweeping algorithm
            // equation: Q(nx/dx + ny/dy) = nxQx/dx + nyQy/dy
            for (int k = j; k < nCellsX+2*nGhost-1; ++k){
                if (std::signbit(phi[i][k+1])==std::signbit(phi[i][k]) && std::signbit(phi[i][k+1]) == std::signbit(phi[i][k+2]) && 
                    std::signbit(phi[i][k])==std::signbit(phi[i+1][k]) && std::signbit(phi[i][k])==std::signbit(phi[i-1][k])){
                    std::array<double,4> qX,qY;
                    if (phi[i][j]<= 0){
                        qX = std::max(f.u[i][k-1],f.u[i][k+1]);
                        qY = std::max(f.u[i-1][k],f.u[i+1][k]);
                    } else {
                        qX = std::min(f.u[i][k-1],f.u[i][k+1]);
                        qY = std::min(f.u[i-1][k],f.u[i+1][k]);
                    }
                    std::array<double,4> Q; 
                    Q = (((fPhix*qX)/dx) + ((fPhiy*qY)/dy))/(fPhix/dx + fPhiy/dy);
                    f.u[i][k] = Q;
                }
            }
        } else {
            for (int k = j; k > 0; --k){
                if (std::signbit(phi[i][k-1])==std::signbit(phi[i][k]) && std::signbit(phi[i][k-1]) == std::signbit(phi[i][k-2]) && 
                    std::signbit(phi[i][k])==std::signbit(phi[i+1][k]) && std::signbit(phi[i][k])==std::signbit(phi[i-1][k])){
                    std::array<double,4> qX,qY;
                    if (phi[i][j]<= 0){
                        qX = std::max(f.u[i][k-1],f.u[i][k+1]);
                        qY = std::max(f.u[i-1][k],f.u[i+1][k]);
                    } else {
                        qX = std::min(f.u[i][k-1],f.u[i][k+1]);
                        qY = std::min(f.u[i-1][k],f.u[i+1][k]);
                    }
                    std::array<double,4> Q;
                    Q = (((fPhix*qX)/dx) + ((fPhiy*qY)/dy))/(fPhix/dx + fPhiy/dy);
                    f.u[i][k] = Q;
                }
            }
        }
    }
    for (std::array<int,2> bPoint : f.interfaceCells){
        int i = bPoint[0];
        int j = bPoint[1];
        double fPhix = phiNormals[i][j][0];
        double fPhiy = phiNormals[i][j][1];
        if (fPhiy <= 0){
            // for f1 sweep in +ve direction if fPhi < 0 and vice versa
            // for f2 sweep in -ve direction if fPhi < 0 and vice versa
            // if i,j+2 is not adjacent (else break):
            // start fast sweeping algorithm
            // equation: Q(nx/dx + ny/dy) = nxQx/dx + nyQy/dy
            for (int k = i; k < nCellsY+2*nGhost-1; ++k){
                if (std::signbit(phi[k][j+1])==std::signbit(phi[k][j]) && std::signbit(phi[k][j+1]) == std::signbit(phi[k][j+2]) && 
                    std::signbit(phi[k][j])==std::signbit(phi[k+1][j]) && std::signbit(phi[k][j])==std::signbit(phi[k-1][j])){
                    std::array<double,4> qX,qY;
                    if (phi[i][j]<= 0){
                        qX = std::max(f.u[k][j-1],f.u[k][j+1]);
                        qY = std::max(f.u[k-1][j],f.u[k+1][j]);
                    } else {
                        qX = std::min(f.u[k][j-1],f.u[k][j+1]);
                        qY = std::min(f.u[k-1][j],f.u[k+1][j]);
                    }
                    std::array<double,4> Q;
                    Q = (((fPhix*qX)/dx) + ((fPhiy*qY)/dy))/(fPhix/dx + fPhiy/dy);
                    f.u[i][k] = Q;
                }
            }
        } else {
            for (int k = i; k < nCellsY+2*nGhost-1; ++k){
                if (std::signbit(phi[k][j+1])==std::signbit(phi[k][j]) && std::signbit(phi[k][j+1]) == std::signbit(phi[k][j+2]) && 
                    std::signbit(phi[k][j])==std::signbit(phi[k+1][j]) && std::signbit(phi[k][j])==std::signbit(phi[k-1][j])){
                    std::array<double,4> qX,qY;
                    if (phi[i][j]<= 0){
                        qX = std::max(f.u[k][j-1],f.u[k][j+1]);
                        qY = std::max(f.u[k-1][j],f.u[k+1][j]);
                    } else {
                        qX = std::min(f.u[k][j-1],f.u[k][j+1]);
                        qY = std::min(f.u[k-1][j],f.u[k+1][j]);
                    }
                    std::array<double,4> Q;
                    Q = (((fPhix*qX)/dx) + ((fPhiy*qY)/dy))/(fPhix/dx + fPhiy/dy);
                    f.u[i][k] = Q;
                }
            }
        }
    }
}

void solver::setGhostFluids(){

    if (splitFlip == 0) {
        maxIter = 301;
    } else {
        maxIter = 21;
    }

    int iter= 0;
    //double maxChange;

    calcPhiGrad();

    //std::cout<< "phi normals calculated" << std::endl;
    do {
        iter++;
        //maxChange = 0;
        setInterface();
        //std::cout << "interface set" << std::endl;
        calcnDotPhiNormals();
        //std::cout << "nDotPhiNormals calculated" << std::endl;
        dtExtrap();
        //std::cout << "dt calculated: " << extrapDt << std::endl;

        for (int i = nGhost; i < nCellsY+nGhost; ++i){
            for (int j = nGhost; j < nCellsX+nGhost; ++j){
                //fluid1.u[i][j] = fluid1.u[i][j] - sign(phi[i][j])*extrapDt*(phi[i][j]>0)*fluid1.nDotGradPhi[i][j];
                //fluid2.u[i][j] = fluid2.u[i][j] - sign(phi[i][j])*extrapDt*(phi[i][j]<=0)*fluid2.nDotGradPhi[i][j];
                std::array<double,4> c1,c2;
                c1 = 0.1*sign(phi[i][j])*extrapDt*(phi[i][j]>0)*fluid1.nDotGradPhi[i][j];
                c2 = 0.1*sign(phi[i][j])*extrapDt*(phi[i][j]<=0)*fluid2.nDotGradPhi[i][j];
                
                //maxChange = std::max(maxChange, *std::max_element(c1.begin(), c1.end()));
                //maxChange = std::max(maxChange, *std::max_element(c2.begin(), c2.end()));

                fluid1.u[i][j] = fluid1.u[i][j] - c1;
                fluid2.u[i][j] = fluid2.u[i][j] - c2;
            }
        }
        //std::cout << maxChange << " " << (maxChange < 1e-3) << std::endl;
        vectorTransBC(fluid1.u);
        vectorTransBC(fluid2.u);
    } while (iter < maxIter);

    setInterface();

    //throw std::runtime_error("ghost fluids set");
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
                fluid.fluxesX[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(fluid.uBarRX[i+1][j],fluid.uBarLX[i+1][j+1],eos)),eos)+LF(fluid.uBarRX[i+1][j],fluid.uBarLX[i+1][j+1],eos)); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = 0.5*(flux(eos->consvToPrim(laxFhalf(fluid.uBarRY[i][j+1],fluid.uBarLY[i+1][j+1],eos)),eos)+LF(fluid.uBarRY[i][j+1],fluid.uBarLY[i+1][j+1],eos)); // fluxes stored in conserved variable form
            }
        }
    }
    /*
    if (&fluid == &fluid2){
        std::cout << "u2" << std::endl;
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
    }
        */
}

std::array<double,4> solver::exactRiemann(std::array<double,4> left,std::array<double,4> right, EOS* eos){
    double d = (direction == XDIR) ? dx : dy;
    riemann solution(eos->get_gamma(),eos->get_gamma(),eos->consvToPrim(left),eos->consvToPrim(right),direction,-0.5*d,0.5*d,0,dt,2,0,0);
    std::array<double,4> result = solution.exctRiemann();
    if (result[RHO] <= 0){
        print_state(eos->consvToPrim(left)); print_state(eos->consvToPrim(right));
        throw std::runtime_error("flux riemann solution rho < 0");
    }
    //if (sqrt(result[UX]*result[UX]) < 1e-4){result[UX] = 0.0;}
    //if (sqrt(result[UY]*result[UY]) < 1e-4){result[UY] = 0.0;}
    return flux(result,eos);
}

std::array<double,4> solver::HLL(std::array<double,4> left,std::array<double,4> right, EOS* eos){ // takes conserved variables
    std::array<double,4> HLLFlux;
    double SL = 0.0, SR = 0.0;
    std::array<double,4> LPrim = eos->consvToPrim(left);
    std::array<double,4> RPrim = eos->consvToPrim(right);

    int dir = (direction == XDIR) ? 1 : 2;

    SL = std::min(LPrim[dir] - eos->calcSoundSpeed(LPrim), RPrim[dir] - eos->calcSoundSpeed(RPrim));
    SR = std::max(LPrim[dir] + eos->calcSoundSpeed(LPrim), RPrim[dir] + eos->calcSoundSpeed(RPrim));

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

    return HLLFlux;
}

std::array<double,4> solver::HLLC(std::array<double,4> left,std::array<double,4> right, EOS* eos){ // takes conserved variables
    std::array<double,4> HLLCFlux;
    double SL = 0.0, SR = 0.0;
    std::array<double,4> LPrim = eos->consvToPrim(left);
    std::array<double,4> RPrim = eos->consvToPrim(right);

    int dir = (direction == XDIR) ? 1 : 2;

    // Pressure - based wavespeed estimation
    
    double pEst = 0.5*(LPrim[PRES]+RPrim[PRES])-0.5*(RPrim[dir]-LPrim[dir])*0.25*(LPrim[RHO]+RPrim[RHO])*(eos->calcSoundSpeed(LPrim)+eos->calcSoundSpeed(RPrim));
    double qL,qR;
    
    qL = (pEst < LPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/LPrim[PRES])-1));
    qR = (pEst < RPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/RPrim[PRES])-1));

    SL = LPrim[dir] - eos->calcSoundSpeed(LPrim)*qL;
    SR = RPrim[dir] + eos->calcSoundSpeed(RPrim)*qR;
    

    // easy wavespeed
    //SL = std::min(LPrim[dir] - eos->calcSoundSpeed(LPrim), RPrim[dir] - eos->calcSoundSpeed(RPrim));
    //SR = std::max(LPrim[dir] + eos->calcSoundSpeed(LPrim), RPrim[dir] + eos->calcSoundSpeed(RPrim));
    //std::cout << SL << " " << SR << std::endl;
    if ((std::isnan(SL) == 1) || (std::isnan(SR) == 1)){
        std::cout << "SL SR " << SL << " " << SR << std::endl;
    }

    // Calculate the numerator of sStar
    double numerator = RPrim[PRES] - LPrim[PRES] + LPrim[RHO] * LPrim[dir] * (SL - LPrim[dir]) - RPrim[RHO] * RPrim[dir] * (SR - RPrim[dir]);
    
    // Calculate the denominator of sStar
    double denominator = LPrim[RHO] * (SL - LPrim[dir]) - RPrim[RHO] * (SR - RPrim[dir]);
    
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
        print_state(LPrim); print_state(RPrim);
        std::cout << "numerator " <<  (RPrim[PRES]-LPrim[PRES]+LPrim[RHO]*LPrim[dir]*(SL-LPrim[dir])-RPrim[RHO]*RPrim[dir]*(SR-RPrim[dir])) << std::endl;
        std::cout << "denominator " << (LPrim[RHO]*(SL-LPrim[dir]) - RPrim[RHO]*(SR-RPrim[dir])) << std::endl; 
        throw std::runtime_error("sStar is nan");
    }
    //

    
    std::array<double,4> Dstar;
    if (direction == XDIR) {
        Dstar = {0, 1, 0, sStar};
    } else {
        Dstar = {0, 0, 1, sStar};
    }

    double pLR = 0.5*(LPrim[PRES]+RPrim[PRES]+LPrim[RHO]*(SL-LPrim[dir])*(sStar-LPrim[dir])+RPrim[RHO]*(SR-RPrim[dir])*(sStar-RPrim[dir]));

    double SdiffL{SL-sStar},SdiffR{SR-sStar};
    if (std::abs(SdiffL) < 1e-8) SdiffL = (SL < 0) ? -1e-8 : 1e-8;
    if (std::abs(SdiffR) < 1e-8) SdiffR = (SR < 0) ? -1e-8 : 1e-8;

    //std::array<double,4> FhllcLeft = (sStar*(SL*left-flux(LPrim,eos))+SL*(LPrim[PRES]+LPrim[RHO]*(SL-LPrim[dir])*(sStar-LPrim[dir]))*Dstar)/(SdiffL);
    //std::array<double,4> FhllcRight = (sStar*(SR*right-flux(RPrim,eos))+SR*(RPrim[PRES]+RPrim[RHO]*(SR-RPrim[dir])*(sStar-RPrim[dir]))*Dstar)/(SdiffR);
    
    // using pLR:
    std::array<double,4> FhllcLeft = (sStar*(SL*left-flux(LPrim,eos))+SL*pLR*Dstar)/(SdiffL);
    std::array<double,4> FhllcRight = (sStar*(SR*right-flux(RPrim,eos))+SR*pLR*Dstar)/(SdiffR);

    
    if (0 <= SL){
        return flux(LPrim,eos);
    } else if (SL < 0 && 0 <= sStar){
        return FhllcLeft;
    } else if (sStar < 0 && 0 <= SR){
        return FhllcRight;
    } else if (SR < 0){
        return flux(RPrim,eos);
    } else {
        throw std::runtime_error("wavespeed condition invalid");
    }
}

void solver::MUSCL(fluid& fluid, EOS* eos){

    calcHalfSlopes(fluid); 

    calcr(fluid);

    calcUBars(fluid);

    updUBars(fluid,eos);

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = fluxRiemannSolver(fluid.uBarRX[i+1][j],fluid.uBarLX[i+1][j+1],eos);
                //fluid.fluxesX[i][j] = HLLC(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1],eos);
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { // why am i getting a Y flux of vx? 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = fluxRiemannSolver(fluid.uBarRY[i][j+1],fluid.uBarLY[i+1][j+1],eos);
                //fluid.fluxesY[i][j] = HLLC(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j],eos);
            }
        }
    }
}

void solver::godunov(fluid& fluid, EOS* eos){

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = fluxRiemannSolver(fluid.u[i+nGhost][j+nGhost-1],fluid.u[i+nGhost][j+nGhost],eos);
                //fluid.fluxesX[i][j] = HLLC(fluid.u[i+nGhost-1][j+1],fluid.u[i+nGhost-1][j+2],eos);
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = fluxRiemannSolver(fluid.u[i+nGhost-1][j+nGhost],fluid.u[i+nGhost][j+nGhost],eos);
                //fluid.fluxesY[i][j] = HLLC(fluid.u[i+1][j+nGhost-1],fluid.u[i+2][j+nGhost-1],eos);
            }
            
        }
    }
}



void solver::transmissiveBC(fluid& fluid){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<fluid.uPlus1[0].size(); ++j){ // sets rows (i)
            fluid.uPlus1[i][j] = fluid.u[nGhost][j];
            fluid.u[i][j] = fluid.u[nGhost][j];
            fluid.uPlus1[nCellsY+2*nGhost-1-i][j] = fluid.u[nCellsY+nGhost-1][j];
            fluid.u[nCellsY+2*nGhost-1-i][j] = fluid.u[nCellsY+nGhost-1][j];
        }
    }
    for (size_t i=0; i<fluid.uPlus1.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            fluid.uPlus1[i][j] = fluid.u[i][nGhost];
            fluid.u[i][j] = fluid.u[i][nGhost];
            fluid.uPlus1[i][nCellsX+2*nGhost-1-j] = fluid.u[i][nCellsX+nGhost-1];
            fluid.u[i][nCellsX+2*nGhost-1-j] = fluid.u[i][nCellsX+nGhost-1];
        }
    }
}

void solver::cylTransmissiveBC(fluid& fluid){
    for (int var =0; var<4;++var){
        if (var == 0 || var == 2 || var == 3){
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<fluid.uPlus1[0].size(); ++j){ // sets rows (ie z)
                    fluid.uPlus1[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.u[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.uPlus1[nCellsY+2*nGhost-1-i][j][var] = fluid.u[nCellsY+nGhost-1][j][var];
                    fluid.u[nCellsY+2*nGhost-1-i][j][var] = fluid.u[nCellsY+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<fluid.uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    fluid.uPlus1[i][nGhost-1-j][var] = fluid.u[i][nGhost+j][var];
                    fluid.u[i][nGhost-1-j][var] = fluid.u[i][nGhost+j][var];
                    fluid.uPlus1[i][nCellsX+2*nGhost-1-j][var] = fluid.u[i][nCellsX+nGhost-1][var];
                    fluid.u[i][nCellsX+2*nGhost-1-j][var] = fluid.u[i][nCellsX+nGhost-1][var];
                }
            }
        } else {
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<fluid.uPlus1[0].size(); ++j){ // sets rows (ie z)
                    fluid.uPlus1[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.u[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.uPlus1[nCellsY+2*nGhost-1-i][j][var] = fluid.u[nCellsY+nGhost-1][j][var];
                    fluid.u[nCellsY+2*nGhost-1-i][j][var] = fluid.u[nCellsY+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<fluid.uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    fluid.uPlus1[i][nGhost-1-j][var] = (-1)*fluid.u[i][nGhost+j][var];
                    fluid.u[i][nGhost-1-j][var] = (-1)*fluid.u[i][nGhost+j][var];
                    fluid.uPlus1[i][nCellsX+2*nGhost-1-j][var] = fluid.u[i][nCellsX+nGhost-1][var];
                    fluid.u[i][nCellsX+2*nGhost-1-j][var] = fluid.u[i][nCellsX+nGhost-1][var];
                }
            }
        }
    }
}

void solver::phiBC(){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<phi[0].size(); ++j){ // sets rows (i)
            phiPlus1[i][j] = phi[nGhost][j];
            phi[i][j] = phi[nGhost][j];
            phiPlus1[nCellsY+2*nGhost-1-i][j] = phi[nCellsY+nGhost-1][j];
            phi[nCellsY+2*nGhost-1-i][j] = phi[nCellsY+nGhost-1][j];
        }
    }
    for (size_t i=0; i<phi.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            phiPlus1[i][j] = phi[i][nGhost];
            phi[i][j] = phi[i][nGhost];
            phiPlus1[i][nCellsX+2*nGhost-1-j] = phi[i][nCellsX+nGhost-1];
            phi[i][nCellsX+2*nGhost-1-j] = phi[i][nCellsX+nGhost-1];
        }
    }
}


void solver::scalarTransBC(std::vector<std::vector<double>>& arr) {
    // Extrapolate in the Y-direction (top and bottom boundaries)
    for (int i = 0; i < nGhost; ++i) {
        for (size_t j = 0; j < arr[0].size(); ++j) { 
            // Bottom boundary (extrapolate using first two interior values)
            arr[i][j] = 2.0 * arr[nGhost][j] - arr[nGhost + 1][j];  
            
            // Top boundary
            arr[nCellsY + 2 * nGhost - 1 - i][j] = 
                2.0 * arr[nCellsY + nGhost - 1][j] - arr[nCellsY + nGhost - 2][j];
        }
    }

    // Extrapolate in the X-direction (left and right boundaries)
    for (size_t i = 0; i < arr.size(); ++i) {
        for (int j = 0; j < nGhost; ++j) {
            // Left boundary
            arr[i][j] = 2.0 * arr[i][nGhost] - arr[i][nGhost + 1];
            
            // Right boundary
            arr[i][nCellsX + 2 * nGhost - 1 - j] = 
                2.0 * arr[i][nCellsX + nGhost - 1] - arr[i][nCellsX + nGhost - 2];
        }
    }
}

template <size_t N>
void solver::vectorTransBC(std::vector<std::vector<std::array<double, N>>>& arr) {
    // Extrapolate in the Y-direction (top and bottom boundaries)
    for (int i = 0; i < nGhost; ++i) {
        for (size_t j = 0; j < arr[0].size(); ++j) { 
            // Bottom boundary (extrapolate using first two interior values)
            arr[i][j] = 2.0 * arr[nGhost][j] - arr[nGhost + 1][j];  
            
            // Top boundary
            arr[nCellsY + 2 * nGhost - 1 - i][j] = 
                2.0 * arr[nCellsY + nGhost - 1][j] - arr[nCellsY + nGhost - 2][j];
        }
    }

    // Extrapolate in the X-direction (left and right boundaries)
    for (size_t i = 0; i < arr.size(); ++i) {
        for (int j = 0; j < nGhost; ++j) {
            // Left boundary
            arr[i][j] = 2.0 * arr[i][nGhost] - arr[i][nGhost + 1];
            
            // Right boundary
            arr[i][nCellsX + 2 * nGhost - 1 - j] = 
                2.0 * arr[i][nCellsX + nGhost - 1] - arr[i][nCellsX + nGhost - 2];
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
        //print_arr(u,RHO);
        //std::cout << "phiX updated" << std::endl;
        for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < phi[0].size() - nGhost+1; ++j) {
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
        for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < phi[0].size() - nGhost+1; ++j) {
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

        for (std::vector<double>::size_type i = nGhost-1; i < phi.size()-nGhost+1; ++i){
            for (std::vector<double>::size_type j = nGhost; j < phi[0].size() - nGhost; ++j) {
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
}

void solver::reinitPhi(){

    for (int iter = 0; iter < 11; iter++){
        std::cout << "iteration " << iter << std::endl;
        //for (size_t i = 0; i<fluid1.interfacePhis.size(); ++i){
        //    phi[fluid1.interfaceCells[i][0]][fluid1.interfaceCells[i][1]] = fluid1.interfacePhis[i];
        //}
        //for (size_t i = 0; i<fluid2.interfacePhis.size(); ++i){
        //    phi[fluid2.interfaceCells[i][0]][fluid2.interfaceCells[i][1]] = fluid2.interfacePhis[i];
        //}
        calcPhiGrad();
        scalarTransBC(phiGrads);
        eikonalDt();
        for (int i = nGhost; i<nCellsY+nGhost; ++i){
            for (int j = nGhost; j < nCellsX+nGhost; ++j){
                phiPlus1[i][j] = phi[i][j] + std::min(dx,dy)*sign(phi[i][j])*(phiGrads[i][j]-1.0);
            }
            //std::cout << std::endl;
        }
        phi = phiPlus1;
        phiBC();
        
    }
    //printScalarField(phi);
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
        for (std::vector<double>::size_type i = nGhost ; i < fluid.u.size()-nGhost; ++i){
            for (std::vector<double>::size_type j = nGhost; j < fluid.u[0].size() - nGhost; ++j) {
                fluid.uPlus1[i][j] = fluid.u[i][j] - (dt / dx) * (fluid.fluxesX[i-nGhost][j-nGhost+1] - fluid.fluxesX[i-nGhost][j-nGhost]); // flux[i + 1] and flux[i] for the update
            }
        }
        fluid.u = fluid.uPlus1;
        //print_arr(fluid.u,RHO);
    } else if (direction == YDIR){
        for (std::vector<double>::size_type i = nGhost; i < fluid.u.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost; j < fluid.u[0].size() - nGhost; ++j) {
                fluid.uPlus1[i][j] = fluid.u[i][j] - (dt / dy) * (fluid.fluxesY[i-nGhost+1][j-nGhost] - fluid.fluxesY[i-nGhost][j-nGhost]); // flux[i + 1] and flux[i] for the update
            }
        }

        fluid.u = fluid.uPlus1;

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
    //return {minBval,minBval,minBval,minBval};
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
            std::cout << std::fixed << std::setprecision(2) << value << " ";
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
            oss << std::fixed << std::setprecision(2) << arr[j][i][var];

            // Append the formatted string to the row
            row += oss.str();
        }
        std::cout << row << std::endl;
        row.clear(); // Clear the row for the next line
    }
};

void solver::print_var(int var) { 
    std::string row;
    for (size_t i = 0; i < fluid1.u.size(); ++i) {
        for (size_t j = 0; j < fluid1.u[0].size(); ++j) {
            row += " ";

            // Use std::ostringstream to format the number to 3 decimal places
            std::ostringstream oss;
            double val = (phi[i][j] <= 0) ? fluid1.u[i][j][var] : fluid2.u[i][j][var];
            oss << std::fixed << std::setprecision(3) << val;

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

std::array<double,4> solver::bilinear(double x, double y, fluid& fluid, std::array<double,2> norm){
    int rows = nCellsY;
    int cols = nCellsX;


    // Compute indices of grid cell containing the point
    int j = ((x - x0) / dx) + nGhost - 0.5;
    int i = ((y - y0) / dy) + nGhost - 0.5;

    int jR,iR;
    if (norm[0] >= 0 && norm[1] >= 0){
        jR = j+1; iR = i+1;
    } else if (norm[0] >= 0){
        jR = j+1; iR = i-1;
    } else if (norm[0] < 0){
        jR = j-1; iR = i-1;
    } else {
        jR = j-1; iR = i+1;
    }

    if (i == 0 || i == (2*nGhost+nCellsY) || j==0 || j == (2*nGhost+nCellsX)){
        return fluid.u[i][j];
    }

    if ((iR >= 0 && iR < rows+2*nGhost) == 0){
        std::cout << "y x ij = " << y << " " << x << " " << i << " " << j << std::endl;
        throw std::runtime_error("i interpolation out of bounds");
    }
    if ((jR >= 0 && jR < cols+2*nGhost) == 0){
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
    double xL = x0 + (j - nGhost + 0.5) * dx;
    double xR = x0 + (j + 1 - nGhost + 0.5) * dx;
    double yL = y0 + (i - nGhost + 0.5) * dy;
    double yR = y0 + (i + 1 - nGhost + 0.5) * dy;
    */
    double xL = x0 + (j - nGhost + 0.5) * dx;
    double xR = x0 + (jR - nGhost + 0.5) * dx;
    double yL = y0 + (i - nGhost + 0.5) * dy;
    double yR = y0 + (iR - nGhost + 0.5) * dy;
    
    std::array<double,4> f00 = fluid.u[i][j];
    std::array<double,4> f10 = fluid.u[iR][j];
    std::array<double,4> f01 = fluid.u[i][jR];
    std::array<double,4> f11 = fluid.u[iR][jR];

    assert(f00[RHO] > 0 && f01[RHO] > 0 && f10[RHO] > 0 && f11[RHO] > 0);
    assert(f00[PRES] > 0 && f01[PRES] > 0 && f10[PRES] > 0 && f11[PRES] > 0);

    
    //std::cout << "interpolation points" << std::endl;
    //std::cout << f00[RHO] << " " << f00[UX]<< " "  << f00[UY]<< " "  << f00[PRES]<< " "  << std::endl;
    //std::cout << f01[RHO] << " " << f01[UX]<< " "  << f01[UY]<< " "  << f01[PRES]<< " "  << std::endl;
    //std::cout << f10[RHO] << " " << f10[UX]<< " "  << f10[UY]<< " "  << f10[PRES]<< " "  << std::endl;
    //std::cout << f11[RHO] << " " << f11[UX]<< " "  << f11[UY]<< " "  << f11[PRES]<< " "  << std::endl;
    
    
    

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


// ---------------------- Time Stepping ---------------------- //

void solver::setDt(){
    double aMax = -1e14;
    //std::cout << "setting dt" << std::endl;
    for (std::vector<double>::size_type i = 0; i < fluid1.u.size(); ++i) {
        for (std::vector<double>::size_type j = 0; j < fluid1.u[0].size(); ++j) {
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

    if (std::fabs(writeTime - time) < 1e-6){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    std::cout << "aMax = " << aMax << ", dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double y_0, double y_1, double t0, double t1, int n_CellsX, int n_CellsY, int n_Ghosts, double c, double g)
    : fluid1(n_CellsX, n_CellsY, n_Ghosts), fluid2(n_CellsX, n_CellsY, n_Ghosts), gamma(g) {
    assert(t1>t0);
    assert(g > 0);

    x0 = x_0; x1 = x_1;
    y0 = y_0; y1 = y_1;
    startTime = t0; endTime = t1;
    nCellsX = n_CellsX;
    nCellsY = n_CellsY;
    nGhost = n_Ghosts;

    resize2D(nCellsX+2*nGhost,nCellsY+2*nGhost,phi);
    resize2D(nCellsX+2*nGhost,nCellsY+2*nGhost,phiPlus1);

    dx = (x1 - x0)/nCellsX;
    dy = (y1 - y0)/nCellsY;


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

template<size_t N>
std::array<double, N> operator-(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
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

template<size_t N>
std::array<double, N> operator+(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
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

template<size_t N>
std::array<double, N> operator*(const double scalar,const std::array<double, N>& lhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = scalar*lhs[i];
        if (std::isnan(result[i]) == 1){
            std::cout << "scalar " << scalar << std::endl;
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("scalar multiply error");
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> operator*(const std::array<double, N>& rhs,const std::array<double, N>& lhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
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

template<size_t N>
std::array<double, N> elementDivide(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
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

template<size_t N>
std::array<double, N> operator/(const std::array<double, N>& lhs, const double scalar) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i]/scalar;
        if (std::isnan(result[i]) == 1){
            std::cout << "LHS:" << lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << std::endl;
            std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
            throw std::runtime_error("scalar divide error");
        }
    }
    return result;
}



