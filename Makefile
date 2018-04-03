# Makefile to produce the executable manis_project1_exe
# -----------------------------------------------------
# Set compiler commands
# -----------------------------------------------------

FC90 = mpif90
FLAGS = -O3 -fcheck=all -fbacktrace

# -----------------------------------------------------
# Primary Rules
# -----------------------------------------------------
main_files = main.f90

all:	main_exe

main_exe: $(main_files)
		  $(FC90) $(FLAGS) $(main_files) -o $@

clean:
		rm *_exe

# END OF MAKEFILE
