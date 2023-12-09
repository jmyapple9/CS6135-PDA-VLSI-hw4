#include "ExampleFunction.h"

// minimize 3*x^2 + 2*x*y + 2*y^2 + 7

ExampleFunction::ExampleFunction(wrapper::Placement &placement)
    : _placement(placement)
{
    chipW = _placement.boundryRight() - _placement.boundryLeft();
    chipH = _placement.boundryTop() - _placement.boundryBottom();
    numModules = _placement.numModules();
    lambda = 0;
    binCut = 15;
    gamma = chipH / 700;
    binTotalNum = binCut * binCut;
    binW = chipW / binCut;
    binH = chipH / binCut;
    binArea = binW * binH;

    tarDensity = 0.0;
    for (size_t i = 0; i < numModules; ++i)
        tarDensity += _placement.module(i).area();
    tarDensity /= (chipW * chipH);
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{

    double sumX1, sumX2, sumY1, sumY2;
    double thetaX, thetaY, aX, bX, aY, bY, dX, dY, ABSdX, ABSdY, c;
    unsigned a, b;

    f = 0.0;
    fill(g.begin(), g.end(), 0.0);
    double *Exp = (double *)calloc(numModules * 4, sizeof(double));

    for (size_t i = 0; i < numModules; ++i)
    {
        Exp[4 * i] = exp(x[2 * i] / gamma);
        Exp[4 * i + 1] = exp(-x[2 * i] / gamma);
        Exp[4 * i + 2] = exp(x[2 * i + 1] / gamma);
        Exp[4 * i + 3] = exp(-x[2 * i + 1] / gamma);
    }

    // LSE
    for (size_t i = 0; i < _placement.numNets(); ++i)
    {
        auto Net = _placement.net(i);
        size_t numPins = Net.numPins();
        sumX1 = sumX2 = sumY1 = sumY2 = 0;
        for (size_t j = 0; j < numPins; ++j)
        {
            int mID = Net.pin(j).moduleId();
            sumX1 += Exp[4 * mID];
            sumX2 += Exp[4 * mID + 1];
            sumY1 += Exp[4 * mID + 2];
            sumY2 += Exp[4 * mID + 3];
        }

        f += log(sumX1) + log(sumX2) + log(sumY1) + log(sumY2);
        // LSE's gradient: sum on g
        for (size_t j = 0; j < numPins; j++)
        {
            wrapper::Module mod = _placement.module(i);
            int mID = Net.pin(j).moduleId();
            if (!mod.isFixed())
            {
                g[2 * mID] += Exp[4 * mID] / (gamma * sumX1);
                g[2 * mID] -= Exp[4 * mID + 1] / (gamma * sumX2);
                g[2 * mID + 1] += Exp[4 * mID + 2] / (gamma * sumY1);
                g[2 * mID + 1] -= Exp[4 * mID + 3] / (gamma * sumY2);
            }
            else
            {
                g[2 * mID] = g[2 * mID + 1] = 0;
            }
        }
    }

    if (lambda == 0)
        return;

    // density (bell-shape)
    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));
    double *grad = (double *)calloc(numModules * 2, sizeof(double));

    for (unsigned binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            wrapper::Module mod = _placement.module(i);
            if (!mod.isFixed())
            {
                aX = 4 / ((binW + mod.width()) * (2 * binW + mod.width()));
                bX = 4 / (binW * (2 * binW + mod.width()));
                aY = 4 / ((binH + mod.height()) * (2 * binH + mod.height()));
                bY = 4 / (binH * (2 * binH + mod.height()));

                a = binIdx % binCut;
                b = binIdx / binCut;
                c = mod.area() / binArea;

                dX = x[2 * i] - ((a + 0.5) * binW + _placement.boundryLeft());
                ABSdX = abs(dX);
                dY = x[2 * i + 1] - ((b + 0.5) * binH + _placement.boundryBottom());
                ABSdY = abs(dY);

                thetaX = thetaByConstrs(ABSdX, mod.width(), binW, aX, bX);
                thetaY = thetaByConstrs(ABSdY, mod.height(), binH, aY, bY);

                binDensity[binIdx] += c * thetaX * thetaY;

                grad[2 * i] = thetaGradByConstrs(ABSdX, mod.width(), binW, aX, bX, dX, thetaY, c);
                grad[2 * i + 1] = thetaGradByConstrs(ABSdY, mod.height(), binH, aY, bY, dY, thetaX, c);
            }
        }

        f += lambda * (binDensity[binIdx] - tarDensity) * (binDensity[binIdx] - tarDensity);
        // bell shape's gradient: sum on g
        for (size_t j = 0; j < numModules; ++j)
        {
            g[2 * j] += lambda * 2 * (binDensity[binIdx] - tarDensity) * grad[2 * j];
            g[2 * j + 1] += lambda * 2 * (binDensity[binIdx] - tarDensity) * grad[2 * j + 1];
        }
    }
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{

    double sumX1, sumX2, sumY1, sumY2;

    f = 0.0;
    double *Exp = (double *)calloc(numModules * 4, sizeof(double));

    // LSE
    for (size_t i = 0; i < numModules; ++i)
    {
        Exp[4 * i] = exp(x[2 * i] / gamma);
        Exp[4 * i + 1] = exp(-x[2 * i] / gamma);
        Exp[4 * i + 2] = exp(x[2 * i + 1] / gamma);
        Exp[4 * i + 3] = exp(-x[2 * i + 1] / gamma);
    }

    for (size_t i = 0; i < _placement.numNets(); ++i)
    {
        auto Net = _placement.net(i);
        size_t numPins = Net.numPins();
        sumX1 = sumX2 = sumY1 = sumY2 = 0;
        for (size_t j = 0; j < numPins; ++j)
        {
            int mID = Net.pin(j).moduleId();
            sumX1 += Exp[4 * mID];
            sumX2 += Exp[4 * mID + 1];
            sumY1 += Exp[4 * mID + 2];
            sumY2 += Exp[4 * mID + 3];
        }

        f += log(sumX1) + log(sumX2) + log(sumY1) + log(sumY2);
    }
    if (lambda == 0)
        return;

    // density (bell-shape)
    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));

    for (unsigned binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            wrapper::Module mod = _placement.module(i);
            if (!mod.isFixed())
            {
                binDensity[binIdx] += bellShapeFunc(i, binIdx, x[2 * i], x[2 * i + 1]);
            }
        }
        f += lambda * (binDensity[binIdx] - tarDensity) * (binDensity[binIdx] - tarDensity);
    }
}

unsigned ExampleFunction::dimension()
{
    return _placement.numModules() * 2; // num_blocks*2
    // each two dimension represent the X and Y dimensions of each block
}

void ExampleFunction::increaseLambda(unsigned offset)
{
    lambda += offset;
}
double ExampleFunction::bellShapeFunc(unsigned i, unsigned binIdx, double oldX, double oldY)
{
    auto mod = _placement.module(i);
    double mW, mH, thetaX, thetaY, dX, dY, ABSdX, ABSdY, aX, bX, aY, bY;
    mW = mod.width(), mH = mod.height();

    aX = 4 / ((binW + mW) * (2 * binW + mW));
    bX = 4 / (binW * (2 * binW + mW));
    aY = 4 / ((binH + mH) * (2 * binH + mH));
    bY = 4 / (binH * (2 * binH + mH));

    unsigned a{binIdx % binCut}, b{binIdx / binCut};
    double c = mod.area() / binArea;

    dX = oldX - ((a + 0.5) * binW + _placement.boundryLeft());
    ABSdX = abs(dX);
    dY = oldY - ((b + 0.5) * binH + _placement.boundryBottom());
    ABSdY = abs(dY);

    thetaX = thetaByConstrs(ABSdX, mW, binW, aX, bX);
    thetaY = thetaByConstrs(ABSdY, mH, binH, aY, bY);
    return c * thetaX * thetaY;
}
/* return: Theta (density) of bell-shape func */
double ExampleFunction::thetaByConstrs(double ABSdX, double mW, double binW, double aX, double bX)
{
    // v1: from teacher's slide
    if (ABSdX <= mW * 0.5 + binW * 0.5)
        return 1 - aX * ABSdX * ABSdX;
    else if (ABSdX <= mW * 0.5 + binW)
        return bX * pow(ABSdX - binW - 0.5 * mW, 2);
    else
        return 0;
    // v2: from NTUplace's paper
    // if (ABSdX <= mW * 0.5 + binW)
    //     return 1 - aX * ABSdX * ABSdX;
    // else if (ABSdX <= mW * 0.5 + binW * 2)
    //     return bX * pow(ABSdX - 2 * binW - 0.5 * mW, 2);
    // else
    //     return 0;
}

/* return: gradient of bell-shape func */
double ExampleFunction::thetaGradByConstrs(double ABSdX, double mW, double binW, double aX, double bX, double dX, double thetaY, double c)
{
    double signX = (dX > 0) ? 1.0 : -1.0;
    // v1: from teacher's slide
    if (ABSdX <= mW * 0.5 + binW * 0.5)
        return c * (-2 * aX * signX * ABSdX) * thetaY;
    else if (ABSdX <= mW * 0.5 + binW)
        return c * 2 * bX * signX * (ABSdX - (binW + 0.5 * mW)) * thetaY;
    else
        return 0;
    // v2: from NTUplace's paper
    // signX = 1;
    // if (ABSdX <= mW * 0.5 + binW)
    //     return c * (-2 * aX * signX * ABSdX) * thetaY;
    // else if (ABSdX <= mW * 0.5 + binW * 2)
    //     return c * 2 * bX * (ABSdX - signX * (2 * binW - 0.5 * mW)) * thetaY;
    // else
    //     return 0;
}
