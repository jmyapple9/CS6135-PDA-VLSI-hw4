#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Wrapper.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(wrapper::Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();
    double bellShapeFunc(unsigned i, unsigned binIdx, double oldX, double oldY);
    double thetaByConstrs(double ABSdX, double mW, double binW, double aX, double bX);
    double thetaGradByConstrs(double ABSdX, double mW, double binW, double aX, double bX, double signX, double thetaY, double c);

    double gamma, chipW, chipH, binW, binH, binArea, tarDensity;
    unsigned cutsize, lambda, numModules, binTotalNum;
private:
    wrapper::Placement &_placement;
};
#endif // EXAMPLEFUNCTION_H
