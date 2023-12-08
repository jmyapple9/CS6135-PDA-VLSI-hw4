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

                // bell-shaped function of bin density smoothing
                if (ABSdX <= (binW / 2 + mod.width() / 2))
                {
                    overlapX = 1 - alphaX * ABSdX * ABSdX;
                }
                else if (ABSdX <= (binW + mod.width() / 2))
                {
                    overlapX = betaX * (ABSdX - (binW + mod.width() / 2)) * (ABSdX - (binW + mod.width() / 2));
                }
                else
                {
                    overlapX = 0;
                }
                // bell-shaped function of bin density smoothing
                if (ABSdY <= (binH / 2 + mod.height() / 2))
                {
                    overlapY = 1 - alphaY * ABSdY * ABSdY;
                }
                else if (ABSdY <= (binH + mod.height() / 2))
                {
                    overlapY = betaY * (ABSdY - (binH + mod.height() / 2)) * (ABSdY - (binH + mod.height() / 2));
                }
                else
                {
                    overlapY = 0;
                }

                if (ABSdX <= (binW / 2 + mod.width() / 2))
                {
                    grad[2 * i] = densityRatio * ((-2) * alphaX * dX) * overlapY;
                }
                else if (ABSdX <= (binW + mod.width() / 2))
                {
                    if (dX > 0)
                    {
                        grad[2 * i] = densityRatio * 2 * betaX * (dX - (binW + mod.width() / 2)) * overlapY;
                    }
                    else
                    {
                        grad[2 * i] = densityRatio * 2 * betaX * (dX + (binW + mod.width() / 2)) * overlapY;
                    }
                }
                else
                {
                    grad[2 * i] = 0;
                }

                if (ABSdY <= (binH / 2 + mod.height() / 2))
                {
                    grad[2 * i + 1] = densityRatio * ((-2) * alphaY * dY) * overlapX;
                }
                else if (ABSdY <= (binH + mod.height() / 2))
                {
                    if (dY > 0)
                    {
                        grad[2 * i + 1] = densityRatio * 2 * betaY * (dY - (binH + mod.height() / 2)) * overlapX;
                    }
                    else
                    {
                        grad[2 * i + 1] = densityRatio * 2 * betaY * (dY + (binH + mod.height() / 2)) * overlapX;
                    }
                }
                else
                {
                    grad[2 * i + 1] = 0;
                }
                // calculate overlap length
                binDensity[binIdx] += densityRatio * overlapX * overlapY;
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

                // bell-shaped function of bin density smoothing
                if (ABSdX <= (binW / 2 + mod.width() / 2))
                {
                    overlapX = 1 - alphaX * ABSdX * ABSdX;
                }
                else if (ABSdX <= (binW + mod.width() / 2))
                {
                    overlapX = betaX * (ABSdX - (binW + mod.width() / 2)) * (ABSdX - (binW + mod.width() / 2));
                }
                else
                {
                    overlapX = 0;
                }
                // bell-shaped function of bin density smoothing
                if (ABSdY <= (binH / 2 + mod.height() / 2))
                {
                    overlapY = 1 - alphaY * ABSdY * ABSdY;
                }
                else if (ABSdY <= (binH + mod.height() / 2))
                {
                    overlapY = betaY * (ABSdY - (binH + mod.height() / 2)) * (ABSdY - (binH + mod.height() / 2));
                }
                else
                {
                    overlapY = 0;
                }
                // calculate overlap length
                binDensity[binIdx] += densityRatio * overlapX * overlapY;
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