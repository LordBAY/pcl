GSoC_PCL
========

Fitting superquadric project for PCL and GSoC 2014

Pre-requisites
-----------

1. Install the latest Ceres version (code was tested with version 1.9). Older versions 
might not have the solver with box constraints that is used in the code for
the minimization operations.

2. For levmar, you might be asked to install f2c, if it is not already installed. It is
   by default in the Ubuntu packages.
   
 
 Code:
 -----
 
1. Checkout the branch sq-fitting of this repo.
2. Go to apps/sq-fitting
3. mkdir build && cd ../build && cmake .. && make
4. In /bin, run demo:
   4.1. With an existing pointcloud: ./demo -r POINTCLOUD.pcd -t OPTIMIZER
   4.2. With the OpenNI grabber ./demo -t OPTIMIZER

   (By now only try the option 4.1. Option 4.2. works if I don't visualize the final results.)
   In OPTIMIZER: 0: CERES 1: LEVMAR.
   
5. The code also produce the cluster pointclouds (CLUSTER_#.pcd), mirrored pointcloud (MIRROR_#.pcd) and fitted SQ cloud (FIT_#.pcd)
