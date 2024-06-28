# Self-Organizing-Maps-for-ice-roughness-analysis
Self Organizing Maps code

%------ Created by R. Gaudioso on 23.12.2023 ------%

--> Unwrapping of iced airfoil geometries through the SOM algorithm proposed by McClain
--> Attempt to extend the formulation to a 3d description of the manifold

Folder structure:

--> som2d/: contains the finished 2d code, featuring functions that compute the manifold and the roughness properties, plus a bit of post-processing
--> som3d/: contains the ongoing development of the 3d code. While the SOM is working, the computation of the arc length in 3d is problematic
--> temp/: folder for testing any new script without touching the established code
--> config_template.m: baseline script to launch a simulation in Matlab. Stores the inputs and the options to control the code and to generate the outputs

%--------------------------------------------------%

som2d/ contains the current stable version. All the useful functions ans scripts are stored and described in this folder.
