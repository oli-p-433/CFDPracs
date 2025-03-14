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
    int nCellsX = 325;
    int nCellsY = 90;
    double x0{0}, x1{0.325}, y0{-0.045}, y1{0.045};
    double startTime = 0.0, endTime = 900e-6;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    //stiffenedGas stiffgas(7.15,3e8,0);
    stiffenedGas idgas1(1.4,1.67,0,0); idealGas idgas2(1.67,1.67);
    std::array<EOS*,2> materials = {&idgas1,&idgas2};
    solver sim(x0,x1,y0,y1,startTime,endTime,nCellsX,nCellsY,2,cour);
    sim.setEOS(materials);

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    sim.dirName = dirname;


    sim.setWriteInterval(50e-6);

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
    sim.flux = [&sim](std::array<double,6> input, EOS* eos, bool direction){
        return sim.fEuler(input,eos,direction);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,6> input){
        return sim.minbee(input);
    }; // because fEuler isnt static

    sim.PRIM = true;
    sim.fluxMethod = [&sim](bool direction) {
        sim.MUSCL(direction); // Change this to sim.SLIC or sim.godunov as needed
    };
    sim.isGodunov = false; // should only be true for HLLCGodunov

    std::cout << "setting initial conditions" << sim.ghosts() << std::endl;

    std::vector< std::vector<std::array<double,6> > > uInit; resizeVector(uInit,nCellsY+2*sim.ghosts(),nCellsX+2*sim.ghosts());

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
            
            

            /*
            double vf{1-(1e-6)}, rho1{1.3765}, rho2{0.1380};
            if (x < 0.1){
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho1*vf,rho2*(1-vf),0.3948,0,1.57});
            } else if (x < 0.125){
                uInit[i][j] = sim.eos[0]->primToConsv({vf,vf,rho2*(1-vf),0,0,1});
            } else if (x < 0.175){
                uInit[i][j] = sim.eos[0]->primToConsv({1-vf,(1-vf),rho2*vf,0,0,1});
            } else {
                uInit[i][j] = sim.eos[0]->primToConsv({vf,vf,rho2*(1-vf),0,0,1});
            }
                */
                

            // Haas + Sturtevant Shock bubble test

            
            double vf{1-(1e-6)}, rho1{1.686}, rho2{1000}, rho3{1.225};
            // sod (x-1)*(x-1)+(y-1)*(y-1) < 0.4*0.4
            if ((x-0.175)*(x-0.175)+(y*y) < 0.025*0.025){
                //uInit1[i][j] = sim.set_vals(0.1380,0,0,1);
                uInit[i][j] = sim.eos[0]->primToConsv({1-vf,rho3*(1-vf),rho2*vf,0,0,101325});
                //uInit[i][j] = sim.eos[0]->primToConsv({vf,vf,rho2*(1-vf),0,0,1});
            } else if (x>=0.225){
                //uInit1[i][j] = sim.set_vals(1.3765,0.3948,0,1.57);
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho1*vf,rho1*(1-vf),-113.5,0,159060});
            } else {
                //uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho3*vf,rho3*(1-vf),0,0,101325});
            }
            

            // water+air shock bubble
            /*
            double vf{1-(1e-6)}, rho1{1.686}, rho2{0.223}, rho3{1.225};
            // sod (x-1)*(x-1)+(y-1)*(y-1) < 0.4*0.4
            if ((x-0.175)*(x-0.175)+(y*y) < 0.025*0.025){
                //uInit1[i][j] = sim.set_vals(0.1380,0,0,1);
                uInit[i][j] = sim.eos[0]->primToConsv({1-vf,rho3*(1-vf),rho2*vf,0,0,101325});
                //uInit[i][j] = sim.eos[0]->primToConsv({vf,vf,rho2*(1-vf),0,0,1});
            } else if (x>=0.225){
                //uInit1[i][j] = sim.set_vals(1.3765,0.3948,0,1.57);
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho1*vf,rho1*(1-vf),-113.5,0,159060});
            } else {
                //uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho3*vf,rho3*(1-vf),0,0,101325});
            }
            */

            // Murrone & Guillard water-air interface
            /*
            double vf{1-(1e-8)}, rho1{1000}, rho2{50};

            if (x < 0.7){
                uInit[i][j] = sim.eos[0]->primToConsv({vf,rho1*vf,rho2*(1-vf),0,0,1e9});
            } else {
                uInit[i][j] = sim.eos[0]->primToConsv({1-vf,rho1*(1-vf),rho2*vf,0,0,1e5});
            }
            */

                
            
            
                

        }
    }

    std::ofstream initFile(dirname + "/0");
    initFile << "x y alpha rho1 rho2 vx vy p" << std::endl;
    for (size_t i = 0; i < uInit.size(); ++i){
        double y = y0 + (i - static_cast<double>(sim.ghosts()) + 0.5) * sim.get_dxdy()[1];
        for (size_t j = 0; j < uInit[0].size(); ++j) {
            double x = x0 + (j - static_cast<double>(sim.ghosts()) + 0.5) * sim.get_dxdy()[0];
            std::array<double,6> prim = sim.eos[0]->consvToPrim(uInit[i][j]);
            initFile << x << " " << y << " " << prim[0] << " " << prim[1] << " " << prim[2] << " " << prim[3] << " " << prim[4] << " " << prim[5] << std::endl;
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