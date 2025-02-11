#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>
#include "solverGFM.H"


// g++ -o fvm fvm.C solver.C


int main(){
    int nCells = 500; int nGhost = 2;
    double x0{0}, x1{1};
    double startTime = 0.0, endTime = 0.05;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    //stiffenedGas stiffgas(7.15,3e8,0);
    
    idealGas idgas1(1.4); idealGas idgas2(1.4);
    fluid f1(nCells,nGhost);fluid f2(nCells,nGhost);
    f1.setEOS(&idgas1);f2.setEOS(&idgas2);
    solverGFM sim(x0,x1,startTime,endTime,nCells,nGhost,cour);

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    sim.dirName = dirname;

    // initialise data locations
    std::vector<std::string> variables = {"rho","v","p"};
    for (size_t i = 0; i<variables.size(); ++i){
        std::string name = dirname + "/" + variables[i];
        std::filesystem::create_directory(name);
    }

    sim.setWriteInterval(0.01);

    // Set the flux function
    sim.flux = [&sim](const EOS* eos,std::array<double,3> input){
        return sim.fEuler(eos,input);
    }; // because fEuler isnt static

    // Set the slope limiter
    sim.slopeLim = [&sim](std::array<double,3> input){
        return sim.minbee(input);
    }; // because fEuler isnt static

    std::vector<std::array<double,3>> u1Init; u1Init.resize(nCells+2*nGhost);
    std::vector<std::array<double,3>> u2Init; u2Init.resize(nCells+2*nGhost);
    std::vector<double> phiInit; phiInit.resize(nCells+2*nGhost);

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (std::vector<double>::size_type i=sim.ghosts(); i<u1Init.size()-sim.ghosts();i++){
        double discPos = 0.5;
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        if (x <= discPos){
            u1Init[i] = f1.eos->primToConsv({1,0,1});
            u2Init[i] = f2.eos->primToConsv({0.125,0,0.1});
        } else {
            u1Init[i] = f1.eos->primToConsv({0.125,0,0.1});
            u2Init[i] = f2.eos->primToConsv({1,0,1});
        }
        phiInit[i] = x-discPos;
    }


    for (std::size_t var = 0; var < variables.size(); ++var) {
        std::ofstream initFile(dirname + "/" + variables[var] + "/0");
        for (std::vector<double>::size_type i = sim.ghosts(); i < u1Init.size() - sim.ghosts(); i++) {
            double x = x0 + (i - sim.ghosts() + 0.5) * sim.get_dx();
            initFile << x << " " << phiInit[i] << " " << f1.eos->consvToPrim(u1Init[i])[var] << " " << f2.eos->consvToPrim(u2Init[i])[var] << std::endl;
        }
        initFile.close();
    }
    

    // Initialising u field: //

    f1.u = u1Init; f2.u = u2Init;
    std::vector<fluid*> fluids = {&f1,&f2};
    std::cout << "setting fluids:" << std::endl;
    sim.setFluids(fluids);
    sim.phiInit(phiInit);


    // Running simulation

    sim.run();
    std::cout << "run finished" << std::endl;

    return 0;
}