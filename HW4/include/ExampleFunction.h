#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Wrapper.hpp"
#include "PlaceData.h"
#include <string.h>
#include <cmath>
#include <iostream>

#define BINCUT 10
class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(wrapper::Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();
    void increaseLambda();

    unsigned numModules{0};
    int binCut{BINCUT}, binTotalNum{BINCUT * BINCUT};
    double eta{500}, boundW{0.0}, boundH{0.0}, binW{0.0}, binH{0.0}, binArea{0.0}, avgDensity{0.0};
    // double *grad{nullptr}, *xExp{nullptr}, *binDensity{nullptr};
    double grad[2 * BINCUT * BINCUT];
    double xExp[4 * BINCUT * BINCUT];
    double binDensity[BINCUT * BINCUT];
private:
    wrapper::Placement &_placement;
    unsigned lambda = 0;
};

#endif // EXAMPLEFUNCTION_H
