#include "ExampleFunction.h"

ExampleFunction::ExampleFunction(wrapper::Placement &placement)
    : _placement(placement)
{
    cout << "Enter ExampleFunction constructor\n";
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
    cout << "Done ExampleFunction constructor\n";
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{

    f = 0.0;
    fill(g.begin(), g.end(), 0.0);

    double *xExp = (double *)calloc(numModules * 4, sizeof(double));
    for (unsigned i = 0; i < numModules; ++i)
    {
        xExp[4 * i] = exp(x[2 * i] / eta);
        xExp[4 * i + 1] = exp(-x[2 * i] / eta);
        xExp[4 * i + 2] = exp(x[2 * i + 1] / eta);
        xExp[4 * i + 3] = exp(-x[2 * i + 1] / eta);
    }

    for (unsigned i = 0; i < _placement.numNets(); ++i)
    {
        double sumX1{0.0}, sumX2{0.0}, sumY1{0.0}, sumY2{0.0};
        for (unsigned j = 0; j < _placement.net(i).numPins(); ++j)
        {
            int mID = _placement.net(i).pin(j).moduleId();
            sumX1 += xExp[4 * mID];
            sumX2 += xExp[4 * mID + 1];
            sumY1 += xExp[4 * mID + 2];
            sumY2 += xExp[4 * mID + 3];
        }

        f += eta * (log(sumX1) + log(sumX2) + log(sumY1) + log(sumY2));

        for (unsigned j = 0; j < _placement.net(i).numPins(); ++j)
        {
            int mID = _placement.net(i).pin(j).moduleId();
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
    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));
    double *grad = (double *)calloc(numModules * 2, sizeof(double));

    double mW{0.0}, mH{0.0}, c{0.0};
    double thetaX{0.0}, thetaY{0.0}, dX{0.0}, dY{0.0}, ABSdX{0.0}, ABSdY{0.0}, aX{0.0}, bX{0.0}, aY{0.0}, bY{0.0};
    for (unsigned a = 0; a < binCut; ++a)
    {
        for (unsigned b = 0; b < binCut; ++b)
        {
            for (unsigned i = 0; i < numModules; ++i)
            {
                mW = _placement.module(i).width();
                mH = _placement.module(i).height();
                if (!_placement.module(i).isFixed())
                {
                    c = _placement.module(i).area() / binArea;

                    dX = _placement.module(i).centerX() - (((double)a + 0.5) * binW + _placement.boundryLeft());
                    ABSdX = abs(dX);
                    dY = _placement.module(i).centerY() - (((double)b + 0.5) * binH + _placement.boundryBottom());
                    ABSdY = abs(dY);
                    aX = 4.0 / ((mW + 2 * binW) * (mW + 4 * binW));
                    aY = 4.0 / ((mH + 2 * binH) * (mH + 4 * binH));
                    bX = 2.0 / (binW * (mW + 4 * binW));
                    bY = 2.0 / (binH * (mH + 4 * binH));

                    thetaX = (ABSdX <= mW * 0.5 + binW) ? (1 - aX * ABSdX * ABSdX) : (ABSdX <= mW * 0.5 + binW * 2) ? (bX * pow(ABSdX - 2 * mW - 2 * binW, 2))
                                                                                                                    : 0;
                    thetaY = (ABSdY <= mH * 0.5 + binH) ? (1 - aY * ABSdY * ABSdY) : (ABSdY <= mH * 0.5 + binH * 2) ? (bY * pow(ABSdY - 2 * mH - 2 * binH, 2))
                                                                                                                    : 0;

                    binDensity[a + binCut * b] += c * thetaX * thetaY;

                    if (!_placement.module(i).isFixed())
                    {
                        double signX{(dX >= 0) ? 1.0 : -1.0}, signY{(dY >= 0) ? 1.0 : -1.0};
                        grad[2 * i] += (ABSdX <= mW * 0.5 + binW) ? (-2 * c * aX * ABSdX * thetaY) : (ABSdX <= mW * 0.5 + binW * 2) ? (2 * c * bX * signX * (ABSdX - 2 * binW - 2 * mW) * thetaY)
                                                                                                                                    : 0;
                        grad[2 * i + 1] += (ABSdY <= mH * 0.5 + binH) ? (-2 * c * signY * aY * ABSdY * thetaX) : (ABSdY <= mH * 0.5 + binH * 2) ? (2 * c * bY * signY * (ABSdY - 2 * binH - 2 * mH) * thetaX)
                                                                                                                                                : 0;
                    }
                }
            }

            f += lambda * pow(binDensity[a + binCut * b] - avgDensity, 2);

            for (unsigned i = 0; i < numModules; ++i)
            {
                g[2 * i] += lambda * 2 * (binDensity[a + binCut * b] - avgDensity) * grad[2 * i];
                g[2 * i + 1] += lambda * 2 * (binDensity[a + binCut * b] - avgDensity) * grad[2 * i + 1];
            }
        }
    }
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

    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));

    for (size_t binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            if (!_placement.module(i).isFixed())
            {
                binDensity[binIdx] += bellShapeFunc(i, binIdx);
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
    lambda += 0.01;
}

double ExampleFunction::bellShapeFunc(size_t i, size_t binIdx)
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

    thetaX = thetaByConstraint(ABSdX, mW, binW, aX, bX);
    thetaY = thetaByConstraint(ABSdY, mH, binH, aY, bY);

    return c * thetaX * thetaY;
}

double ExampleFunction::thetaByConstraint(double ABSdX, double mW, double binW, double aX, double bX)
{
    double thetaX;
    if (ABSdX <= mW * 0.5 + binW)
        thetaX = 1 - aX * ABSdX * ABSdX;
    else if (ABSdX <= mW * 0.5 + binW * 2)
        thetaX = bX * pow(ABSdX - 2 * mW - 2 * binW, 2);
    else
        thetaX = 0;
    return thetaX;
}