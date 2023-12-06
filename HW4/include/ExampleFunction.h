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
private:
    wrapper::Placement &_placement;
};

#endif // EXAMPLEFUNCTION_H
