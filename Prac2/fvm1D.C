#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>
#include "solver1D.H"


// g++ -o fvm fvm.C solver.C


int main(){
    int nCells = 500;
    double x0{0}, x1{1};
    double startTime = 0.0, endTime = 0.25;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    //stiffenedGas stiffgas(7.15,3e8,0);
    idealGas idgas(1.4);
    solver sim(x0,x1,startTime,endTime,nCells,2,cour);
    sim.setEOS(&idgas);

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

    sim.setWriteInterval(0.05);

    // Set the flux function
    sim.flux = [&sim](std::array<double,3> input){
        return sim.fEuler(input);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,3> input){
        return sim.minbee(input);
    }; // because fEuler isnt static

    std::vector<std::array<double,3>> u1Init; u1Init.resize(sim.u1Plus1.size());
    std::vector<std::array<double,3>> u2Init; u2Init.resize(sim.u1Plus1.size());
    std::vector<double> phiInit; phiInit.resize(sim.phiPlus1.size());

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (std::vector<double>::size_type i=sim.ghosts(); i<u1Init.size()-sim.ghosts();i++){
        double discPos = 0.5;
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        if (x <= discPos){
            u1Init[i] = sim.eos->primToConsv({1,0,1});
            u2Init[i] = sim.eos->primToConsv({0.125,0,0.1});
        } else {
            u1Init[i] = sim.eos->primToConsv({0.125,0,0.1});
            u2Init[i] = sim.eos->primToConsv({1,0,1});
        }
        phiInit[i] = x-discPos;
    }


    for (std::size_t var = 0; var < variables.size(); ++var) {
        std::ofstream initFile(dirname + "/" + variables[var] + "/0");
        for (std::vector<double>::size_type i = sim.ghosts(); i < u1Init.size() - sim.ghosts(); i++) {
            double x = x0 + (i - sim.ghosts() + 0.5) * sim.get_dx();
            initFile << x << " " << sim.eos->consvToPrim(u1Init[i])[var] << " " << sim.eos->consvToPrim(u2Init[i])[var] << " " << phiInit[i] << std::endl;
        }
        initFile.close();
    }
    

    // Initialising u field: //

    sim.init(u1Init,u2Init);
    sim.phiInit(phiInit);

    
    // for (int j = 0; j < uInit.size(); j++){
    //     std::cout << uInit[j] << std::endl;
    // }

    sim.run();
    std::cout << "run finished" << std::endl;

    /*
    t = startTime;

    std::string errFileName = dirname + "/err";
    std::ofstream errFile(errFileName);
    if (errFile.is_open()){
        // writing data
        int counter = 0;
        do{
            errFile << t << " " << errs[counter] << std::endl;
            counter++;
            t = t + dt;
        } while (t < endTime);
    }else{
        errFile << "failed to write errors." << std::endl;
    }
    errFile.close();
    */

    return 0;
}