#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Wrapper.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <cfloat>
#include <climits>

class GlobalPlacer
{
public:
    GlobalPlacer(wrapper::Placement &placement);
    void randomPlace(std::vector<double> &x); // An example of random placement implemented by TA
    void place();
    void centerInit(std::vector<double> &sol);

private:
    wrapper::Placement &_placement;
};

#endif // GLOBALPLACER_H
