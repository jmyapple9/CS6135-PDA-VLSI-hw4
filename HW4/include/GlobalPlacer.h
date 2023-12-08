#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Wrapper.hpp"
#include <iostream>
#include <vector>
class GlobalPlacer 
{
public:
    GlobalPlacer(wrapper::Placement &placement);
	void randomPlace(std::vector<double>& x); // An example of random placement implemented by TA
	void place();

private:
    wrapper::Placement &_placement;
};

#endif // GLOBALPLACER_H
