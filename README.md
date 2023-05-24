# BSplineSims
 Dynamic knitting simulations in Julia using bsplines

 Written by Sarah Gonzalez
 Start Date: May 22 2023
 Last Updated: 

 This code follows a similar framework to the static simulations published in [this paper](https://arxiv.org/abs/2302.13467). The stitch backbone is represented as a fifth order b-spline for C3 continuity and we consider an energy functional consisting of bending and compression potentials.

 Currently, the simulation follws one half of the stitch and minimizes the nergy functional using gradient descent to find the minimum energy stitch configuration.

 Things to do:
 1. Update the b-spline basis functions to represent fifth order b-splines with parameter t from 0:5
 2. Calculate the new functionals necessary to implement a stitch height constraint
 3. Implement the stitch height constraint
 4. Update the gradient descent algorithm
 5. Implement new plotting scheme to plot the final stitch as a 3x3 stitch grid
 6. Implement output of results
 7. Figure out how to import a starting configuration


 To run the code:
 1. Install the necessary packages listed in the top of the .jl file
 2. Introduce the initial stitch configuration
 3. Introduce the yarn/fabric parameters
 4. Run simulation in command line

 Outputs
 - energy
 - bending energy
 - compression energy
 - list of control points in x,y,z directions
 - stitch height
 - stitch width
 - stitch length
 - yarn/fabric parameters
 - figure of the final fabric backbone
