# Voronoi Diagram based Dual Contouring Method
## File Structure
- build: The compiled executable of the c++ program will be in this folder

### Code files
- vdc.cpp: the main program
- vdc_func.h/cpp: functions or routines involved in the algorithm
- vdc_type.h: type definitions
- vdc_utilities.h/cpp: general utility functions
- vdc_voronoi.h/cpp: data structures and methods involved with voronoi diagrams used in the program
- vdc_delaunay.h/cpp: data structures and methods involved with delaunay triangulations used in the program
- vdc_debug.h/cpp: debug boolean variables and helper methods for this program
- vdc_grid.h/cpp: Related to the scalar grid data structure used 
- vdc_cube.h/cpp: data struct and methods for cube(centers) processing
- vdc_commandline.h/cpp: Component of reading and parsing the command line arguments
- vdc_globalvar.h/cpp : declaration of the global variables used, //To be improved
- vdc_io.h/cpp: Methods involved with reading input data and write output mesh 
- compExec.py is a python program that takes two executable of the dmr program and compare their output on some input datas
- CMakeList.txt: Needed for compilation if using CMake

