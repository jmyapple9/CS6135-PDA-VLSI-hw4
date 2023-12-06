#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Wrapper.hpp"
#include "PlaceData.h"

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(wrapper::Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();
    void increaseLambda();
private:
    wrapper::Placement &_placement;
    unsigned lambda = 0;
};

#endif // EXAMPLEFUNCTION_H
