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
  // Parse command line inputs
  // ---------------------------------------------------------------------------
  // 
  std::string solnFile;
  double x0, y0, h0, x1, y1, h1, R;
  if (argc == 9){ 
    x0 = std::atof(argv[1]);
    y0 = std::atof(argv[2]);
    h0 = std::atof(argv[3]);
    x1 = std::atof(argv[4]);
    y1 = std::atof(argv[5]);
    h1 = std::atof(argv[6]);
    R = std::atof(argv[7]);
    solnFile = argv[8]; 
  }
  else {
    std::cout << "Usage: DubinsSolver x0 y0 h0 x1 y1 h1 R solnFile" << std::endl;
    return 0;
  }

	//std::string fileName = "dubinsPath.m";

	// define the problem
	RobustDubins::Problem problemStatement;
  problemStatement.set_minTurningRadius(R);

  // Note: heading is in radians
  // Another example: (comment out previous one)
  problemStatement.set_stateInitial(x0,y0,h0); // (x, y, heading)
	problemStatement.set_stateFinal(x1,y1,h1); // (x, y, heading)
  problemStatement.print(); // 

	// run the solver
	RobustDubins::Solver rds;
  rds.set_problemStatement(problemStatement);
  //rds.set_numPts(200); // number of waypoints to return (approx.)
	rds.solve();

  // print results
  rds.print();

  // get waypoints/states of the solution
  //std::vector<double> x,y,h; // 
  //rds.get_optimalWaypoints(x,y,h);

  // alternately, get the RobustDubins::Path object
  // (this contains all the information about the optimal solution)
  RobustDubins::Path optimalPath = rds.get_optimalPath();
  optimalPath.print(); // prints info about this object

  double pathTypeVal;
  std::string pathType = optimalPath.get_pathType();
  if (pathType.compare("LSL")==0){
    pathTypeVal = 1;
  }
  else if (pathType.compare("LSR")==0){
    pathTypeVal = 2;
  }
  else if (pathType.compare("RSL")==0){
    pathTypeVal = 3;
  }
  else if (pathType.compare("RSR")==0){
    pathTypeVal = 4;
  }
  else if (pathType.compare("LRL")==0){
    pathTypeVal = 5;
  }
  else if (pathType.compare("RLR")==0){
    pathTypeVal = 6;
  }
	double aParamUnsigned = optimalPath.get_aParamUnsigned();
	double bParamUnsigned = optimalPath.get_bParamUnsigned();
	double cParamUnsigned = optimalPath.get_cParamUnsigned();
  std::cout << "aParamUnsigned" << aParamUnsigned  << std::endl;
  std::vector< std::vector<double> > data;
  std::vector<double> params = {aParamUnsigned, bParamUnsigned, cParamUnsigned, pathTypeVal};
  data.push_back(params);
  MathTools::writeMatrixData(solnFile, data, "data");

	// optional generate m-file for octave/matlab
	//rds.writeOctaveCommandsToPlotSolution(fileName, 1);

  // optional run m-file in octave
	//MathTools::runOctaveScript(fileName);
  std::cout << " end " << std::endl;
	return 0;
}




