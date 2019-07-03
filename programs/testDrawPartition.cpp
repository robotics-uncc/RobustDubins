#include<string>
#include<vector>
#include<cmath>
#include<iostream>

// custom
#include<RobustDubins_Problem.h>
#include<RobustDubins_Solver.h>
#include<RobustDubins_Path.h>
#include<MathTools.h>

// Intro:
// This function executes the Dubins path algorith over a square grid of
// final endpoints, with a fixed final heading. A plot of the optimal path 
// over such a grid is called a "partition" of the configuration
// space. See Ref. 2,3.

int main() {
  printf("--------------------------------------------\n");
  printf("RobustDubins: Dubins Partition Plotter\n");
  printf("--------------------------------------------\n");
  double hfDeg;
  printf("Final angle (deg.)? ");
  std::cin >> hfDeg;

	// file for plotting commands
	std::string fileName = "dubinsPartitionPlot.m";

	// final heading
	double hFinal = hfDeg*M_PI/180; // final heading
	double boxDim = 8.0; // the dimensions of the square grid of final endpoints
	int nPts = 501; // number of final endpoints per axis

  // vectors of final endpoints
	std::vector< double > xFinalPts = MathTools::linspace(-boxDim, boxDim, nPts);
	std::vector< double > yFinalPts = MathTools::linspace(-boxDim, boxDim, nPts);

	// results are stored in a matrix of size nPts x nPts
  // first intialize one dimension to be of size nPts
	std::vector < std::vector <int> > statusMatrix(nPts); 
	std::vector < std::vector <double> > xPts(nPts); 
	std::vector < std::vector <double> > yPts(nPts); 
  // resize the other dimension: 
	for (int i = 0; i < nPts; ++i){
    statusMatrix[i].resize(nPts);
	  xPts[i].resize(nPts);
	  yPts[i].resize(nPts);
	}

  // a total of npts^2 number of paths will be generated
	double xFinal, yFinal;
  // start for loop that adjusts final endpoint with each iteration
	for (int i = 0; i < nPts ; i++){
		xFinal = xFinalPts.at(i); // set final x coordinate
		for (int j = 0; j < nPts ; j++){
			yFinal = yFinalPts.at(j); // set final y coordinate

	    // define the problem
	    RobustDubins::Problem problemStatement;
			problemStatement.set_stateFinal(xFinal, yFinal, hFinal);	 // (x, y, heading)

			// run the solver
	    RobustDubins::Solver dubinsSolver;
      dubinsSolver.set_problemStatement(problemStatement);
			dubinsSolver.solve();
			
			// save optimal path type and corresponding coordinate
			statusMatrix[i][j] = dubinsSolver.get_optimalSolutionID();
			xPts[i][j] = xFinal;
			yPts[i][j] = yFinal;
		}
		std::cout << " outer loop iter : " << i << std::endl;
	}
	
  // generate octave/matlab .m file
	MathTools::writeMatrixData(fileName, xPts, "xPts");
	MathTools::writeMatrixData(fileName, yPts, "yPts");
	MathTools::writeMatrixData(fileName, statusMatrix, "statusMatrix");
	MathTools::writeVectorData(fileName, xFinalPts, "xFinalPts");
	MathTools::writeVectorData(fileName, yFinalPts, "yFinalPts");
//  MathTools::addCommand(fileName, "[XX,YY] = meshgrid(xFinalPts,yFinalPts);");
//  MathTools::addCommand(fileName, "figure 1; surf(XX,YY,statusMatrix, \'EdgeColor\', \'none\');");
//  MathTools::addCommand(fileName, "view (2)");
//  MathTools::addCommand(fileName, "colorbar;");
  MathTools::addCommand(fileName, "figure 2; surf(xPts,yPts,statusMatrix, \'EdgeColor\', \'none\');");
  MathTools::addCommand(fileName, "view (2)");
  MathTools::addCommand(fileName, "colorbar;");
  // run the .m file in octave
	MathTools::runOctaveScript(fileName);
	return 0;
}




