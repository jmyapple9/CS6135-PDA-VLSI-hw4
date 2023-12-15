#include "ExampleFunction.h"

// minimize 3*x^2 + 2*x*y + 2*y^2 + 7

ExampleFunction::ExampleFunction(wrapper::Placement &placement)
    : _placement(placement)
{
    double
        bTop{_placement.boundryTop()},
        bBottom{_placement.boundryBottom()},
        bRight{_placement.boundryRight()},
        bLeft{_placement.boundryLeft()};
    
    chipW = bRight - bLeft;
    chipH = bTop - bBottom;
    numModules = _placement.numModules();

    lambda = 0;
    cutsize = 15;
    gamma = chipH / 700;
    
    binTotalNum = cutsize * cutsize;
    binW = chipW / cutsize;
    binH = chipH / cutsize;
    binArea = binW * binH;

    tarDensity = 0.0;
    for (size_t mID = 0; mID < numModules; ++mID)
        tarDensity += _placement.module(mID).area();
    tarDensity /= (chipW * chipH);

    /* gamma = 1.0;

    cout << "looking for best gamma...\n";
    while (1)
    {
        double xMaxSum{0}, xMinSum{0}, yMaxSum{0}, yMinSum{0};

        for (size_t nID = 0; nID < _placement.numNets(); ++nID)
        {
            size_t numPins = _placement.net(nID).numPins();
            xMaxSum = xMinSum = yMaxSum = yMinSum = 0;
            for (size_t pID = 0; pID < numPins; ++pID)
            {
                xMaxSum += exp(bTop / gamma);
                xMinSum += exp(bBottom / gamma);
                yMaxSum += exp(bRight / gamma);
                yMinSum += exp(bLeft / gamma);
            }
        }
        if (isinf(xMaxSum) or isnan(xMaxSum)
        or isinf(xMinSum) or isnan(xMinSum)
        or isinf(yMaxSum) or isnan(yMaxSum)
        or isinf(yMinSum) or isnan(yMinSum))
        {
            gamma += 0.1;
        }
        else
        {
            cout << "Find valid gamma = " << gamma << endl;
            // cout << "gamma in origin version" << (_placement.boundryTop() - _placement.boundryBottom())/700 << endl;
            break;
        }
    } */

}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{

    double xMaxSum, xMinSum, yMaxSum, yMinSum;
    double thetaX, thetaY, aX, bX, aY, bY, dX, dY, ABSdX, ABSdY, c;
    unsigned a, b;
    // Initialize f and g
    f = 0.0;
    fill(g.begin(), g.end(), 0.0);
    double *Exp = (double *)calloc(numModules * 4, sizeof(double));
    // LSE
    for (size_t mID = 0; mID < numModules; ++mID)
    {
        Exp[4 * mID] = exp(x[2 * mID] / gamma);
        Exp[4 * mID + 1] = exp(-x[2 * mID] / gamma);
        Exp[4 * mID + 2] = exp(x[2 * mID + 1] / gamma);
        Exp[4 * mID + 3] = exp(-x[2 * mID + 1] / gamma);
    }

    for (size_t nID = 0; nID < _placement.numNets(); ++nID)
    {
        auto Net = _placement.net(nID);
        size_t numPins = Net.numPins();
        xMaxSum = xMinSum = yMaxSum = yMinSum = 0;
        for (size_t pID = 0; pID < numPins; ++pID)
        {
            int mID = Net.pin(pID).moduleId();
            xMaxSum += Exp[4 * mID];
            xMinSum += Exp[4 * mID + 1];
            yMaxSum += Exp[4 * mID + 2];
            yMinSum += Exp[4 * mID + 3];
        }

        f += log(xMaxSum) + log(xMinSum) + log(yMaxSum) + log(yMinSum);
        // LSE's gradient: sum on g
        for (size_t pID = 0; pID < numPins; ++pID)
        {
            int mID = Net.pin(pID).moduleId();
            if (!_placement.module(mID).isFixed())
            {
                g[2 * mID] += Exp[4 * mID] / (gamma * xMaxSum);
                g[2 * mID] -= Exp[4 * mID + 1] / (gamma * xMinSum);
                g[2 * mID + 1] += Exp[4 * mID + 2] / (gamma * yMaxSum);
                g[2 * mID + 1] -= Exp[4 * mID + 3] / (gamma * yMinSum);
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
        for (size_t mID = 0; mID < numModules; ++mID)
        {
            wrapper::Module mod = _placement.module(mID);
            if (!mod.isFixed())
            {
                aX = 4 / ((binW + mod.width()) * (2 * binW + mod.width()));
                bX = 4 / (binW * (2 * binW + mod.width()));
                aY = 4 / ((binH + mod.height()) * (2 * binH + mod.height()));
                bY = 4 / (binH * (2 * binH + mod.height()));

                a = binIdx % cutsize;
                b = binIdx / cutsize;
                c = mod.area() / binArea;

                dX = x[2 * mID] - ((a + 0.5) * binW + _placement.boundryLeft());
                ABSdX = abs(dX);
                dY = x[2 * mID + 1] - ((b + 0.5) * binH + _placement.boundryBottom());
                ABSdY = abs(dY);

                thetaX = thetaByConstrs(ABSdX, mod.width(), binW, aX, bX);
                thetaY = thetaByConstrs(ABSdY, mod.height(), binH, aY, bY);

                binDensity[binIdx] += c * thetaX * thetaY;

                grad[2 * mID] = thetaGradByConstrs(ABSdX, mod.width(), binW, aX, bX, dX, thetaY, c);
                grad[2 * mID + 1] = thetaGradByConstrs(ABSdY, mod.height(), binH, aY, bY, dY, thetaX, c);
            }
        }

        f += lambda * (binDensity[binIdx] - tarDensity) * (binDensity[binIdx] - tarDensity);
        // bell shape's gradient: sum on g
        for (size_t mID = 0; mID < numModules; ++mID)
        {
            g[2 * mID] += lambda * 2 * (binDensity[binIdx] - tarDensity) * grad[2 * mID];
            g[2 * mID + 1] += lambda * 2 * (binDensity[binIdx] - tarDensity) * grad[2 * mID + 1];
        }
    }
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{

    double xMaxSum, xMinSum, yMaxSum, yMinSum;

    f = 0.0;
    double *Exp = (double *)calloc(numModules * 4, sizeof(double));

    // LSE
    for (size_t mID = 0; mID < numModules; ++mID)
    {
        Exp[4 * mID] = exp(x[2 * mID] / gamma);
        Exp[4 * mID + 1] = exp(-x[2 * mID] / gamma);
        Exp[4 * mID + 2] = exp(x[2 * mID + 1] / gamma);
        Exp[4 * mID + 3] = exp(-x[2 * mID + 1] / gamma);
    }

    for (size_t nID = 0; nID < _placement.numNets(); ++nID)
    {
        auto Net = _placement.net(nID);
        size_t numPins = Net.numPins();
        xMaxSum = xMinSum = yMaxSum = yMinSum = 0;
        for (size_t pID = 0; pID < numPins; ++pID)
        {
            int mID = Net.pin(pID).moduleId();
            xMaxSum += Exp[4 * mID];
            xMinSum += Exp[4 * mID + 1];
            yMaxSum += Exp[4 * mID + 2];
            yMinSum += Exp[4 * mID + 3];
        }

        f += log(xMaxSum) + log(xMinSum) + log(yMaxSum) + log(yMinSum);
    }
    if (lambda == 0)
        return;

    // density (bell-shape)
    double *binDensity = (double *)calloc(binTotalNum, sizeof(double));

    for (unsigned binIdx = 0; binIdx < binTotalNum; ++binIdx)
    {
        for (size_t mID = 0; mID < numModules; ++mID)
        {
            wrapper::Module mod = _placement.module(mID);
            if (!mod.isFixed())
            {
                binDensity[binIdx] += bellShapeFunc(mID, binIdx, x[2 * mID], x[2 * mID + 1]);
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

double ExampleFunction::bellShapeFunc(unsigned mID, unsigned binIdx, double oldX, double oldY)
{
    auto mod = _placement.module(mID);
    double mW, mH, thetaX, thetaY, dX, dY, ABSdX, ABSdY, aX, bX, aY, bY;
    mW = mod.width(), mH = mod.height();

    aX = 4 / ((binW + mW) * (2 * binW + mW));
    bX = 4 / (binW * (2 * binW + mW));
    aY = 4 / ((binH + mH) * (2 * binH + mH));
    bY = 4 / (binH * (2 * binH + mH));

    unsigned a{binIdx % cutsize}, b{binIdx / cutsize};
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
