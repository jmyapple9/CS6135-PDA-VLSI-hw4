#include "ExampleFunction.h"

ExampleFunction::ExampleFunction(wrapper::Placement &placement)
    : _placement(placement)
{
    boundW = _placement.boundryRight() - _placement.boundryLeft();
    boundH = _placement.boundryTop() - _placement.boundryBottom();
    numModules = _placement.numModules();

    binTotalNum = binCut * binCut;
    binW = boundW / binCut;
    binH = boundH / binCut;
    binArea = binW * binH;

    avgDensity = 0.0;
    for (unsigned i = 0; i < numModules; ++i)
        avgDensity += _placement.module(i).area();
    avgDensity /= (boundW * boundH);
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    // cout << "Enter evaluateFG\n";
    f = 0.0;
    fill(g.begin(), g.end(), 0.0);
    double *xExp = (double *)calloc(numModules * 4, sizeof(double));

    for (size_t i = 0; i < numModules; ++i)
    {
        xExp[4 * i] = exp(x[2 * i] / eta);
        xExp[4 * i + 1] = exp(-x[2 * i] / eta);
        xExp[4 * i + 2] = exp(x[2 * i + 1] / eta);
        xExp[4 * i + 3] = exp(-x[2 * i + 1] / eta);
    }

    for (size_t i = 0; i < _placement.numNets(); ++i)
    {
        double sumX1{0.0}, sumX2{0.0}, sumY1{0.0}, sumY2{0.0};
        auto Net = _placement.net(i);
        size_t numPins = Net.numPins();
        for (size_t pID = 0; pID < numPins; ++pID)
        {
            int mID = Net.pin(pID).moduleId();
            sumX1 += xExp[4 * mID];
            sumX2 += xExp[4 * mID + 1];
            sumY1 += xExp[4 * mID + 2];
            sumY2 += xExp[4 * mID + 3];
        }
        f += eta * (log(sumX1) + log(sumX2) + log(sumY1) + log(sumY2));

        for (size_t pID = 0; pID < Net.numPins(); ++pID)
        {
            int mID = Net.pin(pID).moduleId();
            if (!_placement.module(mID).isFixed())
            {
                g[2 * mID] += xExp[4 * mID] / (sumX1);
                g[2 * mID] -= xExp[4 * mID + 1] / (sumX2);
                g[2 * mID + 1] += xExp[4 * mID + 2] / (sumY1);
                g[2 * mID + 1] -= xExp[4 * mID + 3] / (sumY2);
            }
        }
    }

    if (lambda == 0)
        return;
    // f = 0.0;// debug!
    // fill(g.begin(), g.end(), 0.0);// debug!


    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));
    double *grad = (double *)calloc(numModules * 2, sizeof(double));

    for (size_t binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            if (!_placement.module(i).isFixed())
            {
                binDensity[binIdx] += bellShapeFunc(i, binIdx, grad, true);
                
            }
        }
        f += lambda * pow(binDensity[binIdx] - avgDensity, 2);

        for (size_t i = 0; i < numModules; ++i)
        {
            g[2 * i] += lambda * 2 * (binDensity[binIdx] - avgDensity) * grad[2 * i];
            g[2 * i + 1] += lambda * 2 * (binDensity[binIdx] - avgDensity) * grad[2 * i + 1];
        }
    }
    // cout << "Done evaluateFG\n";
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{

    f = 0.0;
    double *xExp = (double *)calloc(numModules * 4, sizeof(double));

    for (size_t i = 0; i < numModules; ++i)
    {
        xExp[4 * i] = exp(x[2 * i] / eta);
        xExp[4 * i + 1] = exp(-x[2 * i] / eta);
        xExp[4 * i + 2] = exp(x[2 * i + 1] / eta);
        xExp[4 * i + 3] = exp(-x[2 * i + 1] / eta);
    }

    for (size_t i = 0; i < _placement.numNets(); ++i)
    {
        double sumX1{0.0}, sumX2{0.0}, sumY1{0.0}, sumY2{0.0};
        auto Net = _placement.net(i);
        size_t numPins = Net.numPins();
        for (size_t j = 0; j < numPins; ++j)
        {
            int mID = Net.pin(j).moduleId();
            sumX1 += xExp[4 * mID];
            sumX2 += xExp[4 * mID + 1];
            sumY1 += xExp[4 * mID + 2];
            sumY2 += xExp[4 * mID + 3];
        }
        f += eta * (log(sumX1) + log(sumX2) + log(sumY1) + log(sumY2));
    }

    if (lambda == 0)
        return;
    // f = 0.0; // debug!

    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));

    for (size_t binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            if (!_placement.module(i).isFixed())
            {
                binDensity[binIdx] += bellShapeFunc(i, binIdx, nullptr, false);
            }
        }
        f += lambda * pow(binDensity[binIdx] - avgDensity, 2);
    }
}

unsigned ExampleFunction::dimension()
{
    return 2 * _placement.numModules(); // num_blocks*2
    // each two dimension represent the X and Y dimensions of each block
}

void ExampleFunction::increaseLambda()
{
    lambda += 1000;
}

double ExampleFunction::bellShapeFunc(size_t i, size_t binIdx, double *grad, bool calcGrad)
{
    auto mod = _placement.module(i);
    double a, b, c, mW, mH, thetaX, thetaY, dX, dY, ABSdX, ABSdY, aX, bX, aY, bY;
    a = binIdx % binCut;
    b = binIdx / binCut;
    c = mod.area() / binArea;
    mW = mod.width();
    mH = mod.height();

    dX = mod.centerX() - ((a + 0.5) * binW + _placement.boundryLeft());
    ABSdX = abs(dX);
    dY = mod.centerY() - ((b + 0.5) * binH + _placement.boundryBottom());
    ABSdY = abs(dY);
    aX = 4.0 / ((mW + 2 * binW) * (mW + 4 * binW));
    aY = 4.0 / ((mH + 2 * binH) * (mH + 4 * binH));
    bX = 2.0 / (binW * (mW + 4 * binW));
    bY = 2.0 / (binH * (mH + 4 * binH));

    thetaX = thetaByConstrs(ABSdX, mW, binW, aX, bX);
    thetaY = thetaByConstrs(ABSdY, mH, binH, aY, bY);

    if (calcGrad)
    {
        // c = 1.0; // debug!
        double signX{(dX >= 0) ? 1.0 : -1.0}, signY{(dY >= 0) ? 1.0 : -1.0};
        grad[2 * i] += thetaGradByConstrs(ABSdX, mW, binW, aX, bX, signX, thetaY, c);
        grad[2 * i + 1] += thetaGradByConstrs(ABSdY, mH, binH, aY, bY, signY, thetaX, c);
    }
    return c * thetaX * thetaY;
}

double ExampleFunction::thetaByConstrs(double ABSdX, double mW, double binW, double aX, double bX)
{
    double thetaX;
    if (ABSdX <= mW * 0.5 + binW)
        thetaX = 1 - aX * ABSdX * ABSdX;
    else if (ABSdX <= mW * 0.5 + binW * 2)
        thetaX = bX * pow(ABSdX - 2 * binW - 0.5 * mW, 2);
    else
        thetaX = 0;
    return thetaX;
}
double ExampleFunction::thetaGradByConstrs(double ABSdX, double mW, double binW, double aX, double bX, double signX, double thetaY, double c)
{
    if (ABSdX <= mW * 0.5 + binW)
        return -2 * signX * c * aX * ABSdX * thetaY;
    else if (ABSdX <= mW * 0.5 + binW * 2)
        return 2 * c * bX * signX * (ABSdX - 2 * binW - 0.5 * mW) * thetaY;
    else
        return 0;
}
