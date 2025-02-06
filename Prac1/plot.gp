# Inform the user what variable is being plotted
print "Variable to plot: " . variable_name

# Dynamically construct the file path based on the variable name
files = system("ls ./t/".variable_name."/[0-9]*")
print("Files: " . files)

plot for [file in files] file with lines
pause -1 "Press Enter to exit"
