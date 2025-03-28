#ifndef FLUID_H
#define FLUID_H

#include <array>
#include <vector>

// Templated fluid class with compile-time dimensions
template <std::size_t NX, std::size_t NY, std::size_t NG>
class fluid {
public:
    // Define a shorthand for a 4-component vector
    using Vec4 = std::array<double, 4>;
    
    // u, uPlus1, uPrev, uPrim, fluxesX, fluxesY are defined on a grid 
    // of size (NY + 2*NG) x (NX + 2*NG), with 4 variables at each cell.
    using Array2DVec4 = std::array<std::array<Vec4, NX + 2 * NG>, NY + 2 * NG>;
    Array2DVec4 u, uPlus1, uPrev, uPrim, fluxesX, fluxesY;
    
    // sStarsX is defined on a grid of size NY x (NX+1)
    using Array2DDouble_SX = std::array<std::array<double, NX + 1>, NY>;
    Array2DDouble_SX sStarsX;
    
    // sStarsY is defined on a grid of size (NY+1) x NX
    using Array2DDouble_SY = std::array<std::array<double, NX>, NY + 1>;
    Array2DDouble_SY sStarsY;
    
    // halfSlopesX: grid size NY x (NX + 3) with 4 variables at each cell
    using Array2DVec4_HS_X = std::array<std::array<Vec4, NX + 3>, NY>;
    Array2DVec4_HS_X halfSlopesX;
    
    // halfSlopesY: grid size (NY + 3) x NX with 4 variables at each cell
    using Array2DVec4_HS_Y = std::array<std::array<Vec4, NX>, NY + 3>;
    Array2DVec4_HS_Y halfSlopesY;
    
    // rX: grid size NY x (NX + 2)
    using Array2DVec4_rX = std::array<std::array<Vec4, NX + 2>, NY>;
    Array2DVec4_rX rX;
    
    // rY: grid size (NY + 2) x NX
    using Array2DVec4_rY = std::array<std::array<Vec4, NX>, NY + 2>;
    Array2DVec4_rY rY;
    
    // uBarLX, uBarRX, uBarLupdX, uBarRupdX: grid size NY x (NX + 2)
    using Array2DVec4_uBarX = std::array<std::array<Vec4, NX + 2>, NY>;
    Array2DVec4_uBarX uBarLX, uBarRX, uBarLupdX, uBarRupdX;
    
    // uBarLY, uBarRY, uBarLupdY, uBarRupdY: grid size (NY + 2) x NX
    using Array2DVec4_uBarY = std::array<std::array<Vec4, NX>, NY + 2>;
    Array2DVec4_uBarY uBarLY, uBarRY, uBarLupdY, uBarRupdY;
    
    // uExtrap: grid size (NY + 2*NG) x (NX + 2*NG)
    using Array2DVec4_uExtrap = std::array<std::array<Vec4, NX + 2 * NG>, NY + 2 * NG>;
    Array2DVec4_uExtrap uExtrap;
    
    // The following interface-related containers are left as dynamic containers 
    // because their sizes may not be known at compile time.
    std::vector<std::array<int, 2>> interfaceCells;
    std::vector<std::array<double, 2>> interfaceNormals;
    std::vector<std::array<std::array<double, 2>, 2>> interfacePositions;
    std::vector<std::array<std::array<double, 4>, 2>> interfaceStates;
    std::vector<std::array<double, 4>> riemInterfaceStates;
    std::vector<std::array<std::array<double, 2>, 2>> resolvedVelocities;
    std::vector<std::array<double, 4>> starStates;
    std::vector<double> interfacePhis;
    
    // Constructor: with fixed sizes, no need for resizing
    fluid() {
        // Optionally initialize the arrays here if needed.
        // The memory for all fixed-size arrays is allocated automatically.
    }
};

#endif // FLUID_H