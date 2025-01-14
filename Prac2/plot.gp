files = system("ls ./t/rho/[0-9]*")
print("Files: " . files)

plot for [file in files] file with lines
pause -1 "Press Enter to exit"
