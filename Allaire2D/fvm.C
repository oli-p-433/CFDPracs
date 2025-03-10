#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>
#include "diffuseSolver.H"


// g++ -o fvm fvm.C solver.C


int main(){
    int nCellsX = 200;
    int nCellsY = 200;
    double x0{0}, x1{1}, y0{0}, y1{1};
    double startTime = 0.0, endTime = 0.3;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    //stiffenedGas stiffgas(7.15,3e8,0);
    idealGas idgas1(1.4,1.67); idealGas idgas2(1.67,1.67);
    std::array<EOS*,2> materials = {&idgas1,&idgas2};
    solver sim(x0,x1,y0,y1,startTime,endTime,nCellsX,nCellsY,2,cour);
    sim.setEOS(materials);

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    sim.dirName = dirname;


    sim.setWriteInterval(0.01);

    sim.setBCs = [&sim](fluid& f) {
        sim.get_boundary().transmissiveTopBC(f);
        sim.get_boundary().transmissiveBottomBC(f);
        sim.get_boundary().transmissiveLeftBC(f);
        sim.get_boundary().transmissiveRightBC(f);
        // 1st bool: left/right = IS reflective, 2nd bool: top/bottom = IS reflective
        sim.get_boundary().updateBottomLeftCorner(f,false,false);
        sim.get_boundary().updateTopRightCorner(f,false,false);
        sim.get_boundary().updateBottomRightCorner(f,false,false);
        sim.get_boundary().updateTopLeftCorner(f,false,false);
    };

    std::cout << "set BCs successfully" << std::endl;

    // Set the flux function
    sim.flux = [&sim](std::array<double,5> input, EOS* eos){
        return sim.fEuler(input,eos);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,5> input){
        return sim.minbee(input);
    }; // because fEuler isnt static

    sim.PRIM = true;
    sim.fluxMethod = [&sim]() {
        sim.MUSCL(); // Change this to sim.SLIC or sim.godunov as needed
    };

    std::cout << "setting initial conditions" << sim.ghosts() << std::endl;

    std::vector< std::vector<std::array<double,5> > > uInit; resizeVector(uInit,nCellsY+2*sim.ghosts(),nCellsX+2*sim.ghosts());

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (size_t i=sim.ghosts(); i<uInit.size()-sim.ghosts();++i){
        double y = y0 + (i-sim.ghosts()+0.5)*sim.get_dxdy()[1];
        for (size_t j=sim.ghosts(); j<uInit[0].size()-sim.ghosts();j++){
            double discPos = 0.5;
            double x = x0 + (j-sim.ghosts()+0.5)*sim.get_dxdy()[0];
            /*
            if (x <= 0.05){
                uInit[j] = sim.eos[0]->primToConsv({1,1.3333,0.3535*sqrt(1e5),1.5e5});
            } else if (x <= discPos && x > 0.05){
                uInit[j] = sim.eos[0]->primToConsv({1,1,0,1e5});
            } else {
                uInit[j] = sim.eos[0]->primToConsv({1,1,0,1});
            }
            */

            // Toro tests //

            /*
            double vf{1-(1e-6)}, rho1{1}, rho2{0.125};

            if (x < 0.5){
                uInit[j] = sim.eos[0]->primToConsv({vf,rho1*vf,rho2*(1-vf),0,1});
            } else {
                uInit[j] = sim.eos[0]->primToConsv({1-vf,rho1*(1-vf),rho2*vf,0,0.1});
            }
            */

            // air-helium (wangB)
            
            double vf{1-(1e-6)}, rho1{1.3765}, rho2{0.1380};

            if (x < 0.25){
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho1*vf,rho2*(1-vf),0.3948,1.57});
            } else if (x < 0.4){
                uInit[i][j] = sim.eos[0]->primToConsv({vf,vf,rho2*(1-vf),0,1});
            } else if (x < 0.6){
                uInit[i][j] = sim.eos[0]->primToConsv({1-vf,(1-vf),rho2*vf,0,1});
            } else {
                uInit[i][j] = sim.eos[0]->primToConsv({vf,vf,rho2*(1-vf),0,1});
            }
            
                

        }
    }

    std::ofstream initFile(dirname + "/0");
    initFile << "x y alpha rho1 rho2 v p" << std::endl;
    for (size_t i = 0; i < uInit.size(); ++i){
        double y = y0 + (i - static_cast<double>(sim.ghosts()) + 0.5) * sim.get_dxdy()[1];
        for (size_t j = 0; j < uInit[0].size(); ++j) {
            double x = x0 + (j - static_cast<double>(sim.ghosts()) + 0.5) * sim.get_dxdy()[0];
            std::array<double,5> prim = sim.eos[0]->consvToPrim(uInit[i][j]);
            initFile << x << " " << y << " " << prim[0] << " " << prim[1] << " " << prim[2] << " " << prim[3] << " " << prim[4] << std::endl;
        }
    }
    initFile.close();
    

    // Initialising u field: //

    sim.init(uInit);

    
    // for (int j = 0; j < uInit.size(); j++){
    //     std::cout << uInit[j] << std::endl;
    // }

    sim.run();
    std::cout << "run finished" << std::endl;

    // Calculating exact result with exact riemann solver

    std::array<double,3> sL = {1.333,0.3535*sqrt(1e5),1.5e5};//{1000,0,1e9};
    double gamma1 = 1.4;

    std::array<double,3> sR = {0.1379,0,1e5}; // {50,0,1e5};
    double gamma2 = 1.67;

    double time = 0.0002;
    double discPosition = 0.5;

    riemann solution(gamma1,gamma2,sL,sR,0,1,discPosition,time,100,0,0);
    solution.dirName = dirname;
    // initialise data locations
    std::vector<std::string> riemVariables = {"rhoExact","vExact","pExact"};
    for (size_t i = 0; i<riemVariables.size(); ++i){
        std::string name = dirname + "/" + riemVariables[i];
        std::filesystem::create_directory(name);
    }

    // Exact Riemann solver
    //double pStar = solution.exctRiemann()[2];
    //solution.exactRiemannSolution(pStar,discPosition);

    // end exact calculation

    return 0;
}