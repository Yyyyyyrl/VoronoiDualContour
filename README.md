# Delaunay Triangulation based Mesh Reconstruction Algorithm Implementation

## File Structure
- data: Place the nrrd data files here
- mesh: The resulting mesh of the algorithm will be placed in this folder
- temps: Holds some helping intermediate files during the process
- build: The compiled executable of the c++ program will be in this folder
- Codes are in the outer directory 'getCubeCenter.py' and 'dmr.cpp'

## To Run
- First run 'getCubeCenter.py', this will output the centers of active cubes into a .txt file in the temps folder
- Run the executable 'build/dmr' with the txt file generated in the previous step as input, recompile using cmake and make if any change is made