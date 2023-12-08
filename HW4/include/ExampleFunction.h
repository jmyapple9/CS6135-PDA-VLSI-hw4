#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Wrapper.hpp"
#include "PlaceData.h"
#include <string.h>
#include <cmath>
#include <iostream>
#include <vector>
class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(wrapper::Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();
    void increaseLambda();
    double bellShapeFunc(size_t i, size_t binIdx, double *grad, bool calcGrad);
    double thetaByConstrs(double ABSdX, double mW, double binW, double aX, double bX);
    double thetaGradByConstrs(double ABSdX, double mW, double binW, double aX, double bX, double signX, double thetaY, double c);

    unsigned numModules{0}, binCut{20}, binTotalNum{binCut * binCut};
    double gamma{90}, boundW{0.0}, boundH{0.0}, binW{0.0}, binH{0.0}, binArea{0.0}, avgDensity{0.0};

private:
    wrapper::Placement &_placement;
    double lambda = 1000;
};

#endif // EXAMPLEFUNCTION_H
