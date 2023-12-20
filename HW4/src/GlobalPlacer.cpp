#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"

GlobalPlacer::GlobalPlacer(wrapper::Placement &placement)
    : _placement(placement)
{
}

// Randomly place modules implemented by TA
void GlobalPlacer::randomPlace(vector<double> &result)
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
            _placement.module(i).setPosition(x, y);
        }
        result[2 * i] = _placement.module(i).x();
        result[2 * i + 1] = _placement.module(i).y();
    }
}

void GlobalPlacer::place()
{
    ///////////////////////////////////////////////////////////////////
    // The following example is only for analytical methods.
    // if you use other methods, you can skip and delete it directly.
    //////////////////////////////////////////////////////////////////

    ExampleFunction ef(_placement);        // require to define the object function and gradient function
    vector<double> result(ef.dimension()); // solution vector, size: num_blocks*2
    // randomPlace(result);                   // initialize the solution vector
    centerInit(result);
    double
        bTop{_placement.boundryTop()},
        bBottom{_placement.boundryBottom()},
        bLeft{_placement.boundryLeft()},
        bRight{_placement.boundryRight()};

    unsigned stepSize;
    if ((bTop - bBottom) > 123000 or (bTop - bBottom) < 3000)
    {
        stepSize = (bTop - bBottom) * 6;
    }
    else
    {
        stepSize = (bTop - bBottom) * 5;
    }
    NumericalOptimizer no(ef);
    no.setX(result);               // set initial solution
    no.setStepSizeBound(stepSize); // user-specified parameter

    // no.solve();

    unsigned numModules, EPOCH, numIter;
    numModules = _placement.numModules();
    // numIter = 100;
    EPOCH = 4;
    ef.lambda = 4000;
    for (unsigned epoch = 0; epoch < EPOCH; ++epoch)
    {
        cout << "--------- epoch = " << epoch << "---------\n";
        numIter = 70;
        no.setNumIteration(numIter); // user-specified parameter
        no.solve();                  // Conjugate Gradient solver
        for (unsigned nID = 0; nID < numModules; ++nID)
        {
            wrapper::Module mod = _placement.module(nID);
            double modW = mod.width();
            double modH = mod.height();

            double x_clip = min(bRight - modW, max(bLeft, no.x(2 * nID)));
            double y_clip = min(bTop - modH, max(bBottom, no.x(2 * nID + 1)));

            result[2 * nID] = (mod.isFixed() ? mod.x() : x_clip);
            result[2 * nID + 1] = (mod.isFixed() ? mod.y() : y_clip);

            _placement.module(nID).setPosition(result[2 * nID], result[2 * nID + 1]);
        }
        no.setX(result);
        ef.lambda += 500;
        stepSize *= 1;
        no.setStepSizeBound(stepSize);
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
void GlobalPlacer::centerInit(vector<double> &result)
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
        wrapper::Module mod = _placement.module(nID);
        if (!mod.isFixed())
        {
            mod.setCenterPosition(xcen, ycen);
        }
        result[2 * nID] = mod.x();
        result[2 * nID + 1] = mod.y();
    }
}