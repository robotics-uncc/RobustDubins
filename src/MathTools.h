#ifndef MATH_TOOLS_VECTOR_OPERATORS_H
#define MATH_TOOLS_VECTOR_OPERATORS_H

// standard C++ headers
#include<vector> // std::vector
#include<cmath> // floor, cos, M_PI
#include<algorithm> //std::min_element
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<iomanip>      // std::setprecision

namespace MathTools {

// return the sign of a number
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// Finds minimum element(s) in a vector given a tolerance
std::vector<int> minIndicesWithTolerance(std::vector<double> testVector, 
										                     double tolerance);

// determine (smallest) angular angular distance between two angles (rad)
double polarDistance(double a, double b);

// Computes Euclidean distance between two vectors by taking the p-norm
double distance(std::vector<double> referenceVector, 
				std::vector<double> testVector);

// return the modulus of a number, e.g., mod(3pi/2, pi) = pi/2, forces positive
double        mod(double number, double base); 

// Modify a std::vector by appending a second vector to its end
void append(std::vector<double> &baseVector, std::vector<double> &appendVector);

// Returns the index of the minimum element in a vector
int minElement(std::vector<double> testVector);

// Returns the value of the minimum element in a vector
double minValue(std::vector<double> testVector);

// Generate a sequence of evenly spaced points within prescribed bounds
// if desired skip nSkip several points from the begining
std::vector<double> linspace(double xmin, double xmax, double spacing, 
							               int nSkip);
std::vector<double> linspace(double xmin, double xmax, double spacing);

// check if a double value is an integer
bool isInteger(double value);

// Generates octave commands to write vector data to .m file
void writeVectorData(std::string fileName, std::vector<double> data, 
	                   std::string varName);	

// Generates octave commands to plot a curve (with reference to existing vectrs)
void plot1DCurve(std::string fileName, int figureNumber, std::string xVarName, 
	               std::string yVarName, std::string style);

// Set axis properties e.g. "equal" , "tight"
void axisProperty(std::string fileName, std::string property);

// Generate octave commands to write matrix data to .m file
void writeMatrixData(std::string fileName, 
		                 std::vector< std::vector<double> > &data, 
                     std::string varName);
void writeMatrixData(std::string fileName, 
		                 std::vector< std::vector<int> > &data, 
                     std::string varName);

// Runs a (.m) script in octave
void runOctaveScript(std::string fileName);

// Add command
void addCommand(std::string fileName, std::string command);

} // namespace

#endif
