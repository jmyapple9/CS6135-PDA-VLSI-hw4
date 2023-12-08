#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <cstdlib>
#include <iostream>

GlobalPlacer::GlobalPlacer(wrapper::Placement &placement)
    : _placement(placement)
{
    // util = PlaceData(
    //     placement,
    //     30,
    //     100.0);
    // PlaceData util(30, 100.0);
}

void GlobalPlacer::randomPlace(vector<double> &sol)
{
    srand(0);
    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
    for (unsigned i = 0; i < _placement.numModules(); ++i)
    {
        auto mod = _placement.module(i);
        if (!mod.isFixed())
        {
            double width = mod.width();
            double height = mod.height();
            double x = rand() % static_cast<int>(coreWidth - width) + _placement.boundryLeft();
            double y = rand() % static_cast<int>(coreHeight - height) + _placement.boundryBottom();
            mod.setPosition(x, y);
        }
        sol[2 * i] = mod.x();
        sol[2 * i + 1] = mod.y();
    }
}

void GlobalPlacer::initialPlacement(vector<double> &sol)
{
    double
        bTop{_placement.boundryTop()},
        bBottom{_placement.boundryBottom()},
        bLeft{_placement.boundryLeft()},
        bRight{_placement.boundryRight()};

    double
        xcen{(bLeft + bRight) / 2.0},
        ycen{(bTop + bBottom) / 2.0};

    for (unsigned nID = 0; nID < _placement.numModules(); ++nID)
    {
        auto mod = _placement.module(nID);
        if (!mod.isFixed())
        {
            mod.setCenterPosition(xcen, ycen);
        }
        sol[2 * nID] = mod.x();
        sol[2 * nID + 1] = mod.y();
    }
}

void GlobalPlacer::place()
{
    ///////////////////////////////////////////////////////////////////
    // The following example is only for analytical methods.
    // if you use other methods, you can skip and delete it directly.
    //////////////////////////////////////////////////////////////////

    ExampleFunction ef(_placement); // require to define the object function and gradient function
    vector<double> sol(ef.dimension()); // solution vector, size: num_blocks*2
    // vector<double> sol2(ef.dimension()); // solution vector, size: num_blocks*2
    // vector<double> x(ef.dimension()); // solution vector, size: num_blocks*2
    // each 2 variables represent the X and Y dimensions of a block
    
    // initialize the solution vector
    // initialPlacement(sol);
    randomPlace(sol);
    double
        bTop{_placement.boundryTop()},
        bBottom{_placement.boundryBottom()},
        bLeft{_placement.boundryLeft()},
        bRight{_placement.boundryRight()};

    NumericalOptimizer no(ef);
    no.setX(sol);                  // set initial solution
    no.setStepSizeBound(3000);     // user-specified parameter
    // no.solve();

    unsigned numModules = _placement.numModules();
    unsigned EPOCH = 1;
    for (unsigned epoch = 0; epoch < EPOCH; ++epoch)
    {
        cout << "--------- epoch = " << epoch << "---------\n";
        unsigned numIter = (epoch==0)?300:50;
        no.setNumIteration(numIter);       // user-specified parameter
        no.solve(); // Conjugate Gradient solver
        for (unsigned nID = 0; nID < numModules; ++nID)
        {
            wrapper::Module mod = _placement.module(nID);
            double modW = mod.width();
            double modH = mod.height();

            double x_clip = min(bRight - modW, max(bLeft, no.x(2 * nID)));
            double y_clip = min(bTop - modH, max(bBottom, no.x(2 * nID + 1)));

            sol[2 * nID] = (mod.isFixed() ? mod.x() : x_clip);
            sol[2 * nID + 1] = (mod.isFixed() ? mod.y() : y_clip);

            _placement.module(nID).setPosition(sol[2 * nID], sol[2 * nID + 1]);
        }
        no.setX(sol);
        ef.increaseLambda();
    }

    cout << "Objective: " << no.objective() << "\n";
    ////////////////////////////////////////////////////////////////

    // An example of random placement implemented by TA.
    // If you want to use it, please uncomment the folllwing 1 line.
    // randomPlace();

    /* @@@ TODO
     * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
     * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
     * 3. For the bin density model, you could refer to the lecture notes
     * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN
     * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + beta*BinDensity()"
     * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
     * 7. Set the initial vector x in place(), set step size, set #iteration, and call the solver like above example
     * */
}
