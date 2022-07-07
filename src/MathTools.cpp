#include<MathTools.h>

double MathTools::mod(double number, double base){
  return number - std::floor(number/base)*base;
}

int MathTools::mod(int number, int base){
	return number - std::floor(number/base)*base;
}

unsigned int MathTools::mod(unsigned int number, unsigned int base){
	return number - std::floor(number/base)*base;
}

double MathTools::distance(std::vector<double> referenceVector, 
		                       std::vector<double> testVector){
	std::vector<double>::size_type vectorSize = referenceVector.size();
	double distanceSquared = (referenceVector[0]-testVector[0]) * 
							             (referenceVector[0]-testVector[0]);
	for (std::vector<double>::size_type i = 1; i < vectorSize ; i++){
		distanceSquared = distanceSquared + (referenceVector[i]-testVector[i]) 
						                          * (referenceVector[i]-testVector[i]);
	}	
	return std::sqrt(distanceSquared);
}

double MathTools::polarDistance(double a, double b){
  // make both angles positive
  a = MathTools::mod(a, 2.0*M_PI);
  b = MathTools::mod(b, 2.0*M_PI);
  // calculate angular distance without wrapping around 2pi
  double delta = std::abs(b-a);
  // check if wrapping around is shorter
  if (delta > 2.0*M_PI - delta){
    return 2.0*M_PI - delta;
  }
  else {
    return delta;
  }
}

double MathTools::polarDistanceSigned(double a, double b){
  // make both angles positive
  a = MathTools::mod(a, 2.0*M_PI);
  b = MathTools::mod(b, 2.0*M_PI);
  if (b > a){
    if ( b-a > M_PI ){
      return -(2.0*M_PI - b + a);
    }
    else {
      return b - a;
    }
  }
  else {
    if ( b-a < -M_PI ){
      return 2.0*M_PI - a + b;
    }
    else {
      return b - a;
    }
  }
}

std::vector<int> MathTools::minIndicesWithTolerance(
                      std::vector<double> testVector,
                      // arma::vec testVector,
                                                double tolerance){
   // initialize vector of min indices
   std::vector<int> minIndices;
   // store index of smallest value in testVector
   //arma::uword minInd;
   int minInd;
   // store smallest value for reference
   double minValueRef = std::numeric_limits<double>::max();
   for(int i = 0; i < testVector.size(); i++)
     if(minValueRef > testVector[i]) { minValueRef = testVector[i]; minInd = i;}
   //double minValueRef = testVector.min(minInd);
   minIndices.push_back(minInd);
   // set this element to largest possible value since it has been recorded
   testVector[minInd] = std::numeric_limits<double>::max();
   //int i = 1;
   for(int i = 0; i < testVector.size(); i++){
     if (std::abs(testVector[i] - minValueRef) <= tolerance)
       minIndices.push_back(i);
   }
  
  //MathTools::debugPrintVector(minIndices,true);
   return minIndices;
}

void MathTools::append(std::vector<double> &baseVector, std::vector<double> &appendVector){
	// Modify a std::vector by appending a second vector to its end
	baseVector.insert(std::end(baseVector), std::begin(appendVector), 
	                  std::end(appendVector));
}

double MathTools::sum(std::vector<double> xValues){
	double currentSum = xValues[0];
	for (int i = 1; i < xValues.size(); i++){
		currentSum = currentSum + xValues[i];
	}
	return currentSum;
}


// -----------------------------------------------------------------------------
// plot1DCurve
void MathTools::plot1DCurve(std::string fileName, int figureNumber, 
		                        std::string xVarName, std::string yVarName, 
                            std::string style){	
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << "figure " << figureNumber << " \n ";
	outfile << "plot ("<< xVarName << " , " << yVarName <<",'" << style << 
		"', 'linewidth',3); hold on; \n";
	outfile << "grid on;" << "\n ";
	outfile << "set(gca, 'linewidth', 3, 'fontsize', 16) \n";
    	// close file
    	outfile.close();
	return;
}
// -----------------------------------------------------------------------------
// plot1DCurve
void MathTools::plot1DCurve(std::string fileName, int figureNumber, 
		                        std::string xVarName, std::string style){	
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << "figure " << figureNumber << " \n ";
	outfile << "plot ("<< xVarName << ",'" << style << 
		"', 'linewidth',3); hold on; \n";
	outfile << "grid on;" << "\n ";
	outfile << "set(gca, 'linewidth', 3, 'fontsize', 16) \n";
    	// close file
    	outfile.close();
	return;
}
// -----------------------------------------------------------------------------
// plot1DCurve 
void MathTools::plot1DCurve(std::string fileName, int figureNumber, 
		                        std::vector<double> xData, 
                            std::vector<double> yData, 
		                        std::string xVarName, std::string yVarName, 
                            std::string colorStyle){	
	MathTools::writeVectorData(fileName, xData, xVarName);
	MathTools::writeVectorData(fileName, yData, yVarName);
	MathTools::plot1DCurve(fileName,figureNumber,xVarName,yVarName,colorStyle);
}

// plot1DCurve for std::vector data, with data input
void MathTools::plot1DCurve(std::string fileName, int figureNumber, 
		                        std::vector<double> xData, 
		                        std::string xVarName, 
                            std::string colorStyle){	
	MathTools::writeVectorData(fileName, xData, xVarName);
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << "figure " << figureNumber << " \n ";
	outfile << "plot ("<< xVarName <<",'" << colorStyle
          << "', 'linewidth',3); hold on; \n";
	outfile << "grid on;" << "\n ";
	outfile << "set(gca, 'linewidth', 3, 'fontsize', 16) \n";
	// close file
	outfile.close();
	return;
}

void MathTools::axisProperty(std::string fileName, std::string property){
	// open file for writing
	std::ofstream outfile(fileName, std::ofstream::app);
	outfile << "axis(\"" << property << "\") \n";
	// close file
    	outfile.close();
}

// -----------------------------------------------------------------------------
// writeVectorData 
void MathTools::writeVectorData(std::string fileName, std::vector<double> data, 
		                            std::string varName){
  // open for writing
	std::ofstream outfile(fileName, std::ofstream::app); 
  // write coordinates
	outfile << std::setprecision(17) << varName << "= [ " << data[0]; 
    for (int i = 1; i < data.size(); i++){
        outfile << " , " << data[i];
    }	
	outfile << "]; \n";
  outfile.close();	   	// close file
}
// -----------------------------------------------------------------------------
// writeVectorData 
void MathTools::writeVectorData(std::string fileName, std::vector<int> data, 
		                            std::string varName){
  // open for writing
	std::ofstream outfile(fileName, std::ofstream::app); 
  // write coordinates
	outfile << std::setprecision(17) << varName << "= [ " << data[0]; 
    for (int i = 1; i < data.size(); i++){
        outfile << " , " << data[i];
    }	
	outfile << "]; \n";
  outfile.close();	   	// close file
}
// -----------------------------------------------------------------------------
// writeVectorData 
void MathTools::writeVectorData(std::string fileName, std::vector<double> data, 
		                            std::string varName, int precision){
  // open for writing
	std::ofstream outfile(fileName, std::ofstream::app); 
  // write coordinates
	outfile << std::setprecision(precision) << varName << "= [ " << data[0]; 
    for (int i = 1; i < data.size(); i++){
        outfile << " , " << data[i];
    }	
	outfile << "]; \n";
  outfile.close();	   	// close file
}
 
std::vector< std::vector<double> > MathTools::generateCircularArc( 
		                                              std::vector<double> arcParams, 
                                                  int numPts, int direction){
	double xc = arcParams[0];
	double yc = arcParams[1];
	double R = arcParams[2];
	double tInit = arcParams[3];
	double tFinal = arcParams[4];
	std::vector<double> t = MathTools::polarspace(tInit, tFinal, numPts, 
		direction);
	std::vector<double> xData, yData;
	for (int i = 0; i < numPts; i++){
		xData.push_back(xc + R*cos(t[i]));
		yData.push_back(yc + R*sin(t[i]));
	}
	std::vector< std::vector<double> > arcPoints;
	arcPoints.push_back(xData);
	arcPoints.push_back(yData);
	return arcPoints;
}
 
std::vector< std::vector<double> > MathTools::generateCircle( 
		std::vector<double> circleParams, int numPts){
	circleParams.push_back(0);
	circleParams.push_back(2.0*M_PI);
	return MathTools::generateCircularArc(circleParams, numPts, 1.0);
}

// linspace
std::vector<double> MathTools::linspace(double xmin, double xmax, 
		                                    double spacing, int nSkip){
	// Generate a squence of evenly spaced points within prescribed bounds
	// if desired skip nSkip several points from the begining
	// note: spacing cannot always be exactly enforced
	std::vector<double> outputVector;	
	double delx = xmax - xmin;
	double dir = MathTools::sign(delx); //arma direction of interval
	// check if request spacing is feasible
	if (std::abs(spacing) >= std::abs(delx)){ // not feasible, manually set 
		outputVector.push_back(xmin);
		outputVector.push_back(xmax);
		std::cout << "MathTools::linspace warning: requested spacing is "
				      << "larger than interval" << std::endl;
	}
	else { // spacing is feasible 
		double nptsDecimal = std::abs(delx/spacing)+1.0;
		int nptsInt = (int)nptsDecimal;
		// check if number of skip points is feasible
		if (nSkip >= nptsInt){
			std::cout << "MathTools::linspace error: number of skip points is "
					      << "not feasible, reduce spacing" << std::endl;
		}
		else {
			// check if spacing divides interval into a whole number
			if (MathTools::isInteger(nptsDecimal)){
				for (int i = nSkip; i < nptsInt ; i++){
					outputVector.push_back(xmin + dir*spacing*i);
				}
			}
			else { // set the last point to xmax, overriding spacing request
				for (int i = 0; i < nptsInt-1; i++){
					outputVector.push_back(xmin + dir*spacing*i);
				}
				outputVector.push_back(xmax);
			}
		}
	}
	return outputVector;
}

// linspace
std::vector<double> MathTools::linspace(double xmin, double xmax, 
		                                    double spacing){
	return MathTools::linspace(xmin, xmax, spacing, 0);
}

// linspace
std::vector<double> MathTools::linspace(double xmin, double xmax, int npts, 
		                                    int nSkip){
  // This function !MUST return npts 
  if (npts - nSkip < 2){
    std::cout << "MathTools::linspace error: insufficient num. pts requested" 
              << std::endl; 
    printf("(xmin,xmax,npts,nSkip) = (%3.3f, %3.3f, %d, %d) \n",xmin,xmax,npts,nSkip);
    std::vector<double> output(2);
    output.push_back(xmin);
    output.push_back(xmax);
  }
  else {
    std::vector<double> outputVector;
	  if (xmin == xmax){ // all the points are the same
		  for (int i = nSkip; i < npts; i++){
			  outputVector.push_back(xmin);
		  }
		  return outputVector;
	  }
	  else { 
	    double delx = xmax - xmin;
	    double dir = MathTools::sign(delx); // direction of interval
      double spacing = std::abs(delx/((double)npts-1.0));
      for (int i = nSkip; i < npts - 1; i++){
					outputVector.push_back(xmin + dir*spacing*i);
      }
      outputVector.push_back(xmax);
		  return outputVector;
	  }
  }
}

// linspace
std::vector<double> MathTools::linspace(double xmin, double xmax, int npts){
  int nSkip = 0;
	return MathTools::linspace(xmin, xmax, npts, nSkip);
}

// polarspace
std::vector<double> MathTools::polarspace(double angleInitial, 
		                                      double angleFinal, 
                                          int npts, int direction){
	// assumes less than one revolution (i.e. the interval requested is < 2PI)
	double angleInitialMod = MathTools::mod(angleInitial, 2.0*M_PI);
	double angleFinalMod;
	if (angleFinal == 2.0*M_PI){ // admit the case of a full circle
		angleFinalMod = angleFinal;
	}
	else {
		angleFinalMod = MathTools::mod(angleFinal, 2.0*M_PI);
	}
	if (direction == 1) { //counter-clockwise
		if (angleFinalMod < angleInitialMod){ 
			// the interval includes 2*PI, so we make the final angle > 2.0*PI
			angleFinalMod = angleFinalMod + 2.0*M_PI; 
		}
	}
	else if (direction == -1){ //clockwise
		if (angleFinalMod > angleInitialMod){
			// by making the final angle negative, we force the CW direction
			angleFinalMod = angleFinalMod - 2.0*M_PI;
		}
	}
	else {
		std::cout << "MathTools::polarspace, "
			 	      << "Error: unknown direction given" << std::endl;
    std::cout << "angleInitial : " << angleInitial << std::endl;
    std::cout << "angleFinal : " << angleFinal << std::endl;
    std::cout << "npts : " << npts << std::endl;
    std::cout << "direction : " << direction << std::endl;
	}
	return MathTools::linspace(angleInitialMod, angleFinalMod, npts);
}

// -----------------------------------------------------------------------------
// isInteger
bool MathTools::isInteger(double value){
	if (value == floor(value))
		return true;
	else 
		return false;
}

// runOctaveScript
void MathTools::runOctaveScript(std::string fileName){
	std::string runCommand = "octave --force-gui --persist " + fileName;
	system(runCommand.c_str());
}

void MathTools::runOctaveScript(std::string fileName, std::string options){
	std::string runCommand = "octave " + options + " " + fileName;
	system(runCommand.c_str());
}

// -----------------------------------------------------------------------------
// addCommand
void MathTools::addCommand(std::string fileName, std::string command){
  std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << command << "\n"; // write command
  outfile.close(); // close file
}

// -----------------------------------------------------------------------------
// randomNumberUniformDistribution
double MathTools::randomNumberUniformDistribution(double xmin, double xmax){
  srand (std::random_device()()); /* initialize random seed: */
  double delta = (xmax - xmin);    
  return (double)rand()/RAND_MAX*delta+xmin; // random number from (xmin,xmax)
}

// -----------------------------------------------------------------------------
// randomNumberUniformDistribution
std::vector<double> MathTools::randomNumberUniformDistribution(double xmin, 
		                                                           double xmax, 
                                                               int npts){
  std::vector<double> xvect;
  for (int i=0 ; i <npts; i++){
      xvect.push_back(MathTools::randomNumberUniformDistribution(xmin, xmax));
  }
  return xvect;
}

double MathTools::fmodPos( double val, double mod ){
  while ( val < 0 ){
    val = val + mod;
  }
  return fmod(val,mod);
}

// -----------------------------------------------------------------------------
// writeMatrixData
void MathTools::writeMatrixData(std::string fileName, 
		                            std::vector< std::vector<double> > &data, 
                                std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= [ ";
	// assume that data is given as a vector of row vectors
	std::vector< double > currentRow;
    for (int i = 0; i < data.size(); i++){ // for each row
		currentRow = data.at(i);
		  for (int j = 0; j < currentRow.size(); j++){
			  outfile << " , " << currentRow.at(j);
		  }
      outfile << " ; "; // end of row
    }	
	outfile << "]; \n";
    outfile.close();	   	// close file
}
// -----------------------------------------------------------------------------
// writeMatrixData
void MathTools::writeMatrixData(std::string fileName, 
		                            std::vector< std::vector<int> > &data, 
                                std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= [ ";
	// assume that data is given as a vector of row vectors
	std::vector< int > currentRow;
    for (int i = 0; i < data.size(); i++){ // for each row
		currentRow = data.at(i);
		for (int j = 0; j < currentRow.size(); j++){
			outfile << " , " << currentRow.at(j);
		}
        outfile << " ; "; // end of row
    }	
	outfile << "]; \n";
    outfile.close();	   	// close file
}
