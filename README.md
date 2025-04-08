# BSplineSims
 Dynamic knitting simulations in Julia using bsplines

 Written by Sarah Gonzalez
 Start Date: May 22 2023
 Last Updated: Jul 18 2023

 NOTE AS OF 04/08/2025: The fixed height and fixed width julia files are my playground for each scenario and don't actually do what they promise in the name. Only stockinettefixed.jl actually runs as indicated in the title, with both dimensions fixed. You would need to accurately calculate the change in energy as the height/width changes to implement free boundaries. Additional issue: there is currently no constraint on the derivatives at the endpoints, which you should do to get accurate results. I completed testing without it, but only out of laziness. Realistically, the derivatives should be constrained to zero at the endpoints to make the entire fabric continuously differentiable.

 This code follows a similar framework to the static simulations published in [this paper](https://arxiv.org/abs/2302.13467). The stitch backbone is represented as a fifth order b-spline for C3 continuity and we consider an energy functional consisting of bending and compression potentials.

 Currently, the simulation follws one half of the stitch and minimizes the nergy functional using gradient descent to find the minimum energy stitch configuration.

 Things to do:
 1. Calculate the new functionals necessary to implement a stitch height constraint
 2. Implement the stitch height constraint
 3. Update the gradient descent algorithm
 4. Implement new plotting scheme to plot the final stitch as a 3x3 stitch grid
 5. Implement output of results
 6. Figure out how to import a starting configuration


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
