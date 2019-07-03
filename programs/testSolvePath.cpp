// A. Wolek, 15-Mar-2015
// testDubins.cpp


#include<string>
#include<vector>
#include<cmath>
#include<iostream>

// custom
#include<RobustDubins_Problem.h>
#include<RobustDubins_Solver.h>
#include<RobustDubins_Path.h>
#include<MathTools.h>

int main ( int argc , char ** argv ) {
  printf("--------------------------------------------\n");
  printf("RobustDubins: Dubins Solver\n");
  printf("--------------------------------------------\n");
  printf("Use default start pt and radius? (0-no, 1-yes):");
  int defaultChoice; 
  std::cin >> defaultChoice;
  double R, x0, y0, h0deg;  
  if ( defaultChoice == 0){
    printf("Initial X-Position: ");
    std::cin >> x0;
    printf("Initial Y-Position: ");
    std::cin >> y0;
    printf("Initial H-Orientation (deg): ");
    std::cin >> h0deg;
    printf("Turn Radius: ");
    std::cin >> R;
  }
  double x1, y1, h1deg;  
  printf("Final X-Position: ");
  std::cin >> x1;
  printf("Final Y-Position: ");
  std::cin >> y1;
  printf("Final H-Orientation (deg): ");
  std::cin >> h1deg;

	// define the problem
	RobustDubins::Problem problemStatement;
	problemStatement.set_stateFinal(x1,y1,h1deg*M_PI/180.0);
  if ( defaultChoice == 0){
    problemStatement.set_stateInitial(x0,y0,h0deg*M_PI/180.0);
    problemStatement.set_minTurningRadius(R);
  }
  problemStatement.print(); 

	// run the solver
	RobustDubins::Solver rds;
  rds.set_problemStatement(problemStatement);
	rds.solve();
  rds.print();

  // get the optimal path 
  RobustDubins::Path optimalPath = rds.get_optimalPath();
  optimalPath.print(); 

	// optional generate m-file for octave/matlab
  std::string fileName = "dubinsTest.m";
	rds.writeOctaveCommandsToPlotSolution(fileName, 1);

  // optional run m-file in octave
	MathTools::runOctaveScript(fileName);
  std::cout << " end " << std::endl;
	return 0;
}




