# 2D FDTD code (Akash Dasgupta) 
### For submission in the Advanced Computing module

### Supervisor for project: 
* Simom Hannah

### University of Bristol, 2019-2020

* Contains source code for FDTD simulation engine, developed as part of high performance computing project
* Uses MPI for parallel processing
* Has built in support for gausian pulsed and contimous wave sources, in planewave and point source configurations.
* Has 3 built in scenirios, with tunable parameters:
    * Lensing: Both a planar convex and biconvex lense are provided, allows to control radius of curvature and diameter of lens
    * Double slit difraction: allows for control of slit width and seperation
    * Evanescent tunnelling: Places 2 triangles of high refractive index, seperated by a controlable gap, to span the grid. Waves can undergo total internal reflection, and tunnelling depending on the gap size.

## Building:

* Requirements: 
    * A vaugely unix-esque OS (it'll probably work on Windows too but I make no promises). Tested Os: Ubuntu 18.04.4 LTS (bionic), Scientific Linux release 6.4 (Carbon)
    * gcc and g++ compiler. Tested versions on which it'll definately work: gcc-7.4.0, gcc-4.8.5
    * MPI headers. Versions tested: 2.1.1
    * Burning passion for modeling optical systems
* Make file is provided, should build with 'make' command
* 'make clean' WILL DELETE CONTENTS OF DATA FOLDER, so be aware when rebuilding
* There are a lot of parameters that can be tweaked in main, both as variables and as defines

## Visualising 

* Can be visualised through provided vis.py, already set up for mentioned pre-built simulations
* Requires a whole load of modules, backends, etc, so best to just use an anaconda install
* Reads data from data folder, and saves video of the data
