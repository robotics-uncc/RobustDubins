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

int main() {
  std::string pathType;
  double turnRadius, a, b, c;
	std::string fileName = "dubinsPath.m";
  RobustDubins::Path path;
  printf("--------------------------------------------\n");
  printf("RobustDubins: Dubins Path Plotter\n");
  printf("--------------------------------------------\n");
  printf("Specify path type: LSL, LSR, RSL, RSR, LRL, RLR? : ");
  std::cin  >> pathType;
  path.set_pathType(pathType);
  printf("Turn radius : ");
  std::cin  >> turnRadius;
  path.set_minTurnRadius(turnRadius);
  printf("Angle subtended by first segment (deg): ");
  std::cin  >> a;
  path.set_aParamUnsigned(a*M_PI/180.0);
  if ( pathType.compare("LRL")==0 || pathType.compare("RLR")==0 ){
    printf("Angle subtended by second segment (deg): ");
    std::cin  >> b;
    path.set_bParamUnsigned(b*M_PI/180.0);
  }
  else {
    printf("Length of second segment : ");
    std::cin  >> b;
    path.set_bParamUnsigned(b);
  }
  printf("Angle subtended by third segment (deg): ");
  std::cin  >> c;
  path.set_cParamUnsigned(c*M_PI/180.0);
  path.computeEndpoint();
  path.computePathHistory();
  path.print();
  path.writePathOctavePlotCommands(fileName, 1, "x", "y", "bo-");
  MathTools::runOctaveScript(fileName);
	return 0;
}




