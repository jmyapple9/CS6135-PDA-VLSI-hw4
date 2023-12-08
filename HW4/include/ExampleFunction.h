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
    void increaseLambda(unsigned offset);
    double bellShapeFunc(size_t i, size_t binIdx, double oldX, double oldY);
    double thetaByConstrs(double ABSdX, double mW, double binW, double aX, double bX);
    double thetaGradByConstrs(double ABSdX, double mW, double binW, double aX, double bX, double signX, double thetaY, double c);

    double gamma, boundW, boundH, binW, binH, binArea, tarDensity;
    unsigned binCut{15}, lambda{0}, numModules, binTotalNum;
    // int binNumPerEdge;
    // int binNum;
    // double gamma;
    // double Width;
    // double Height;
    // double dBinW;
    // double dBinH;
    // double targetDensity;
    // double *g_temp;
    // double *dExp;
    // vector<double> binDensity;
private:
    wrapper::Placement &_placement;
};
#endif // EXAMPLEFUNCTION_H
