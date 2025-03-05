# Inform the user what variable is being plotted
print "Variable to plot: " . variable_name

# Dynamically construct the file path based on the variable name
files = system("ls ./".directory_name."/".variable_name."/[0-9]*")
file = system("ls ./".directory_name."/".variable_name."/".timestamp."")
#print("Files: " . files)
last_file = word(files, words(files))
exact_file = system("ls ./".directory_name."/".variable_name."Exact/[0-9]*")
unset key

set style line 1 lt 1 lw 2  # Regular line (solid)
set style line 2 lt 2 lw 2 dashtype 2   # Dotted line (lt 2 is the line type for dotted)
#set yrange [0:1.5]
#set size 1.625,0.45
#set origin -0.25,0.25
#plot for [file in files] file using 1:($4 < 0 ? $2 : 1/0) with lines,for [file in files] file using  1:($4 > 0 ? $3 : 1/0) with lines, \
#for [file in files] file using 1:4 with lines
splot file using 1:2:($5 < 0 ? $3 : 1/0) with pm3d notitle, file using  1:2:($5 > 0 ? $4 : 1/0) with pm3d notitle, \
#exact_file with lines, \
#last_file using 1:4 with lines
pause -1 "Press Enter to exit"

