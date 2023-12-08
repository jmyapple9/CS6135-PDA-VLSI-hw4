#include "ExampleFunction.h"

// minimize 3*x^2 + 2*x*y + 2*y^2 + 7

ExampleFunction::ExampleFunction(wrapper::Placement &placement)
    : _placement(placement)
{
    boundW = _placement.boundryRight() - _placement.boundryLeft();
    boundH = _placement.boundryTop() - _placement.boundryBottom();
    numModules = _placement.numModules();
    gamma = boundH / 500;
    binTotalNum = binCut * binCut;
    binW = boundW / binCut;
    binH = boundH / binCut;
    binArea = binW * binH;

    tarDensity = 0.0;
    for (size_t i = 0; i < numModules; ++i)
        tarDensity += _placement.module(i).area();
    tarDensity /= (boundW * boundH);
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{

    double sumX1, sumX2, sumY1, sumY2;
    double overlapX, overlapY, alphaX, betaX, alphaY, betaY, dX, dY, ABSdX, ABSdY;
    double densityRatio;

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
        // gradient: g
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

    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));
    double *grad = (double *)calloc(numModules * 2, sizeof(double));

    for (unsigned binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            wrapper::Module mod = _placement.module(i);
            if (!mod.isFixed())
            {
                alphaX = 4 / ((binW + mod.width()) * (2 * binW + mod.width()));
                betaX = 4 / (binW * (2 * binW + mod.width()));
                alphaY = 4 / ((binH + mod.height()) * (2 * binH + mod.height()));
                betaY = 4 / (binH * (2 * binH + mod.height()));

                densityRatio = mod.area() / (binW * binH);
                unsigned a{binIdx % binCut}, b{binIdx / binCut};
                dX = x[2 * i] - ((a + 0.5) * binW + _placement.boundryLeft());
                ABSdX = abs(dX);
                dY = x[2 * i + 1] - ((b + 0.5) * binH + _placement.boundryBottom());
                ABSdY = abs(dY);

                overlapX = thetaByConstrs(ABSdX, mod.width(), binW, alphaX, betaX);
                overlapY = thetaByConstrs(ABSdY, mod.height(), binH, alphaY, betaY);

                binDensity[binIdx] += densityRatio * overlapX * overlapY;

                grad[2 * i] = thetaGradByConstrs(ABSdX, mod.width(), binW, alphaX, betaX, dX, overlapY, densityRatio);
                grad[2 * i + 1] = thetaGradByConstrs(ABSdY, mod.height(), binH, alphaY, betaY, dY, overlapX, densityRatio);
            
                // calculate overlap length
                // binDensity[binIdx] += bellShapeFunc(i, binIdx, grad, true);

            }
        }
        // calculate (lambda * Σ((Db(x, y)-Tb)^2)) part in objective function
        f += lambda * (binDensity[binIdx] - tarDensity) * (binDensity[binIdx] - tarDensity);

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
    double overlapX, overlapY, alphaX, betaX, alphaY, betaY, dX, dY, ABSdX, ABSdY;
    double densityRatio;
    // f is the objective cost function
    f = 0.0;
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
    }
    if (lambda == 0)
        return;

    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));

    for (unsigned binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t i = 0; i < numModules; ++i)
        {
            wrapper::Module mod = _placement.module(i);
            if (!mod.isFixed())
            {
                alphaX = 4 / ((binW + mod.width()) * (2 * binW + mod.width()));
                betaX = 4 / (binW * (2 * binW + mod.width()));
                alphaY = 4 / ((binH + mod.height()) * (2 * binH + mod.height()));
                betaY = 4 / (binH * (2 * binH + mod.height()));

                densityRatio = mod.area() / (binW * binH);
                unsigned a{binIdx % binCut}, b{binIdx / binCut};

                dX = x[2 * i] - ((a + 0.5) * binW + _placement.boundryLeft());
                ABSdX = abs(dX);
                dY = x[2 * i + 1] - ((b + 0.5) * binH + _placement.boundryBottom());
                ABSdY = abs(dY);

                overlapX = thetaByConstrs(ABSdX, mod.width(), binW, alphaX, betaX);
                overlapY = thetaByConstrs(ABSdY, mod.height(), binH, alphaY, betaY);
                // calculate overlap length
                // binDensity[binIdx] += densityRatio * overlapX * overlapY;

                binDensity[binIdx] += bellShapeFunc(i, binIdx, x[2 * i], x[2 * i + 1]);

            }
        }
        // calculate (lambda * Σ((Db(x, y)-Tb)^2)) part in objective function
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
double ExampleFunction::bellShapeFunc(size_t i, size_t binIdx, double oldX, double oldY)
{
    
    auto mod = _placement.module(i);
    double a, b, c, mW, mH, thetaX, thetaY, dX, dY, ABSdX, ABSdY, aX, bX, aY, bY;
    aX = 4 / ((binW + mod.width()) * (2 * binW + mod.width()));
    bX = 4 / (binW * (2 * binW + mod.width()));
    aY = 4 / ((binH + mod.height()) * (2 * binH + mod.height()));
    bY = 4 / (binH * (2 * binH + mod.height()));

    a = binIdx % binCut;
    b = binIdx / binCut;
    c = mod.area() / (binW * binH);
    // unsigned a{binIdx % binCut}, b{binIdx / binCut};

    dX = oldX - ((a + 0.5) * binW + _placement.boundryLeft());
    ABSdX = abs(dX);
    dY = oldY - ((b + 0.5) * binH + _placement.boundryBottom());
    ABSdY = abs(dY);

    thetaX = thetaByConstrs(ABSdX, mod.width(), binW, aX, bX);
    thetaY = thetaByConstrs(ABSdY, mod.height(), binH, aY, bY);
    return c * thetaX * thetaY;

}

double ExampleFunction::thetaByConstrs(double ABSdX, double mW, double binW, double aX, double bX)
{
    if (ABSdX <= (binW / 2 + mW / 2))
    {
        return 1 - aX * ABSdX * ABSdX;
    }
    else if (ABSdX <= (binW + mW / 2))
    {
        return bX * (ABSdX - (binW + mW / 2)) * (ABSdX - (binW + mW / 2));
    }
    else
    {
        return 0;
    }
    // v1: from teacher's slide
    // if (ABSdX <= mW * 0.5 + binW * 0.5)
    //     return 1 - aX * ABSdX * ABSdX;
    // else if (ABSdX <= mW * 0.5 + binW)
    //     return bX * pow(ABSdX - binW - 0.5 * mW, 2);
    // else
    //     return 0;
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
    if (ABSdX <= (binW / 2 + mW / 2))
    {
        return c * ((-2) * aX * dX) * thetaY;
    }
    else if (ABSdX <= (binW + mW / 2))
    {
        if (dX > 0)
        {
            return c * 2 * bX * (dX - (binW + mW / 2)) * thetaY;
        }
        else
        {
            return c * 2 * bX * (dX + (binW + mW / 2)) * thetaY;
        }
    }
    else
    {
        return 0;
    }
    // v1: from teacher's slide
    // if (ABSdX <= mW * 0.5 + binW * 0.5)
    //     return c * (-2 * aX * signX * ABSdX) * thetaY;
    // else if (ABSdX <= mW * 0.5 + binW)
    //     return c * 2 * bX * signX * (ABSdX - (binW + 0.5 * mW)) * thetaY;
    // else
    //     return 0;
    // v2: from NTUplace's paper
    // signX = 1;
    // if (ABSdX <= mW * 0.5 + binW)
    //     return c * (-2 * aX * signX * ABSdX) * thetaY;
    // else if (ABSdX <= mW * 0.5 + binW * 2)
    //     return c * 2 * bX * (ABSdX - signX * (2 * binW - 0.5 * mW)) * thetaY;
    // else
    //     return 0;
}
