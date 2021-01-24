# Robust Dubins C++ Library

The ability to compensate for disturbances when following a Dubins path can be improved by artificially increasing the turn radius used in path planning. Consider a Dubins car with a nominal minimum turn radius subject to unknown disturbances of bounded magnitude. A Dubins path that is planned using an inflated turn radius (scaled as a function of the disturbance upper bound) is feasible, in the sense that it could be exactly followed if the disturbance were known. In practice, one typically uses feedback control to compensate for an unknown disturbance. In this case, Monte Carlo simulations suggest  that these feasible Dubins paths (planned using the inflated turn radius) minimize the mean, final cross-track error when the path is followed using a standard guidance algorithm. Furthermore, the simulation results indicate that there is no appreciable benefit in planning paths with a turn radius larger than the inflated turn radius.


<p align="center"> 
<img src="https://coefs.uncc.edu/awolek/files/2021/01/RobustDubins.png" width="300">
</p>

## Dependencies:
- Ubuntu 16.04
- cmake
- (Optional) Octave/MATLAB
  for generating .m files for plotting

## Build instructions:
- Run `./build.sh` to create a library in `./lib` and executables in `./bin`

## Usage instructions:
- Either link to and use the library in your own code, or run the command line DubinsSolver command with syntax:
`DubinsSolver x0 y0 h0 x1 y1 h1 R solnFile`
- See test programs under `./programs` for examples of how to define a `RobustDubins::Problem` to supply to the `RobustDubins::Solver` and obtain a `RobustDubins::Path`. Briefly, the syntax is:
```
// define the problem
RobustDubins::Problem problemStatement;
problemStatement.set_stateInitial(x1,y1,h1deg*M_PI/180.0);	
problemStatement.set_stateFinal(x1,y1,h1deg*M_PI/180.0);	

// run the solver
RobustDubins::Solver rds; 
rds.set_problemStatement(problemStatement);
rds.solve();
rds.print(); // print output to screen (all candidates)

// get the optimal path 
RobustDubins::Path optimalPath = rds.get_optimalPath();
optimalPath.print(); // print details of optimal path
```
- `clean.sh` removes all compiled files, build files, and temporary files 

## References:

1. Tang, G., Wang, Z., & Williams, A. L. (1998). On the construction of an optimal feedback control law for the shortest path problem for the Dubins car-like robot. In Proceedings of the 30th Southeastern Symposium on Systems Theory (pp. 280–284). 
https://doi.org/10.1109/SSST.1998.660075

2. Wolek, A., & Woolsey, C. (2015). Feasible Dubins paths in presence of unknown, unsteady velocity disturbances. Journal of Guidance, Control, and Dynamics, 38(4), 782–787. 
https://doi.org/10.2514/1.G000629

## Contact:
Artur Wolek, wolek@umd.edu
