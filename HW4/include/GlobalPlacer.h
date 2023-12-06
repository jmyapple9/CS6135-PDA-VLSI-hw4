#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Wrapper.hpp"
#include "PlaceData.h"
#include <vector>
class GlobalPlacer
{
public:
    GlobalPlacer(wrapper::Placement &placement);

    void randomPlace(std::vector<double>& x); // An example of random placement implemented by TA
    void place();
    void initialPlacement(std::vector<double>& x);

private:
    wrapper::Placement &_placement;
    PlaceData util;
};

#endif // GLOBALPLACER_H
