# Inform the user what variable is being plotted
print "Variable to plot: " . variable_name

# Dynamically construct the file path based on the variable name
files = system("ls ./t/".variable_name."/[0-9]*")
print("Files: " . files)

# unset key

set style line 1 lt 1 lw 2  # Regular line (solid)
set style line 2 lt 2 lw 2 dashtype 2   # Dotted line (lt 2 is the line type for dotted)
set yrange [-5:5]

plot for [file in files] file using 1:2 with lines, for [file in files] file using 1:3 with lines
pause -1 "Press Enter to exit"
