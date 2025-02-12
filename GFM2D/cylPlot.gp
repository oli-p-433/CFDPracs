# Set up the terminal and output file (adjust this based on your environment)
set terminal pngcairo enhanced font 'Arial,10' size 800, 600
set output 'cylindrical_field_plot.png'

# Set the labels for axes
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'

# Set the title of the plot
set title 'Cylindrical Data with Scalar Field φ'

# Use a 3D plot style (points or lines, depending on your preference)
set style data points

# Set color map for phi (adjust the range of phi if needed)
set palette defined ( 0 "blue", 1 "green", 2 "yellow", 3 "red" )

# Create the plot using cylindrical data in (r, z, phi) form
splot '0.10' using ($1*cos($3)):$($1*sin($3)):$2:3 with points palette title 'Scalar Field φ'
