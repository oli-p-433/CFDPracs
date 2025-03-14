#!/usr/bin/gnuplot
# Example usage:
#   gnuplot -e "datafolder='my_data'; timestep='001.dat'; var='rho1'" plot_swap.gp

# Use defaults if variables aren’t passed
if (!exists("datafolder")) { datafolder = "." }
if (!exists("timestep")) { timestep = "data.dat" }
if (!exists("var")) { var = "rho1" }

# Build the full file path
filename = sprintf("%s/%s", datafolder, timestep)

# Map the input variable to the proper column indices.
# Assumes the data file has columns:
#   1: x, 2: y, 3: phi, then for each variable you have:
#      e.g. "rho1" in column 4 and "rho2" in column 5;
#           "vx1"  in column 6 and "vx2"  in column 7; etc.
if (var eq "rho1") {
    col1 = 4; col2 = 5; alt = "rho2"
} else {
    if (var eq "vx1") {
        col1 = 6; col2 = 7; alt = "vx2"
    } else {
        print "Error: Variable not recognized: " . var
        exit
    }
}

# Set up labels and title
set xlabel "x"
set ylabel "y"
set zlabel var
set title sprintf("Plot of %s (using %s when phi < 0, %s otherwise)", var, alt, var)

# Here the ternary operator is reversed:
#   When phi (column 3) is less than 0, use column col2 (e.g., rho2),
#   otherwise use column col1 (e.g., rho1)
splot filename using 1:2:( $3 < 0 ? column(col2) : column(col1) ) \
      with pm3d notitle

# Optional: pause the plot window until a key is pressed
pause -1 "Press any key to exit"