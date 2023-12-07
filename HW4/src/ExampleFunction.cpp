#include "ExampleFunction.h"

// minimize 3*x^2 + 2*x*y + 2*y^2 + 7
// double binDensity[BINCUT];
// double grad[BINCUT * 2];
// double xExp[BINCUT * 4];

ExampleFunction::ExampleFunction(wrapper::Placement &placement)
    : _placement(placement)
{
    cout << "Enter ExampleFunction constructor\n";
    boundW = _placement.boundryRight() - _placement.boundryLeft();
    boundH = _placement.boundryTop() - _placement.boundryBottom();
    numModules = _placement.numModules();

    // lambda = 0;
    // eta = 500;
    binCut = 10;

    binTotalNum = binCut * binCut;
    binW = boundW / binCut;
    binH = boundH / binCut;
    binArea = binW * binH;

    // binDensity.resize(binTotalNum);
    // binDensity.assign(binTotalNum, 0.0);
    // grad.resize(2 * numModules);
    // grad.assign(2 * numModules, 0.0);
    // xExp.resize(4 * numModules);
    // xExp.assign(4 * numModules, 0.0);

    // binDensity = new double[binTotalNum]();
    // grad = new double[numModules * 2]();
    // xExp = new double[numModules * 4]();
    avgDensity = 0.0;
    for (unsigned i = 0; i < numModules; ++i)
        avgDensity += _placement.module(i).area();
    avgDensity /= (boundW * boundH);
    // cout << "eta: " << eta <<endl;
    // cout << "binCut: " << binCut <<endl;
    // cout << "binTotalNum: " << binTotalNum <<endl;
    // cout << "binW: " << binW <<endl;
    // cout << "binH: " << binH <<endl;
    // cout << "boundW: " << boundW <<endl;
    // cout << "boundH: " << boundH <<endl;
    // cout << "binArea: " << binArea <<endl;
    // cout << "avgDensity: " << avgDensity <<endl;
    // cout << "numModules: " << numModules <<endl;
    cout << "Done ExampleFunction constructor\n";
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    // f = 3 * x[0] * x[0] + 2 * x[0] * x[1] + 2 * x[1] * x[1] + 7; // objective function
    // g[0] = 6 * x[0] + 2 * x[1];                                  // gradient function of X
    // g[1] = 2 * x[0] + 4 * x[1];                                  // gradient function of Y
    // cout << "Enter: evaluateFG\n";

    fill(g.begin(), g.end(), 0.0);
    f = 0.0;
    // for(int i = 0; i < 2 * numModules; i++) grad[i]=0;
    // for(int i = 0; i < binTotalNum; i++) binDensity[i]=0;
    // fill(grad.begin(), grad.end(), 0.0);
    // fill(binDensity.begin(), binDensity.end(), 0.0);
    // memset(binDensity, 0.0, sizeof(double)*binTotalNum);
    // memset(grad, 0.0, sizeof(double)*numModules*2);

    vector<double> binDensity(binTotalNum, 0.0);
    vector<double> grad(numModules * 2, 0.0);
    // cout << "binDensity.size(): " <<binDensity.size() << endl;
    // cout << "grad.size(): " <<grad.size() << endl;
    // cout << "xExp.size(): " <<xExp.size() << endl;
    // cout << "During evaluateFG... "<< endl;
    fill(grad.begin(), grad.end(), 0.0);
    fill(binDensity.begin(), binDensity.end(), 0.0);
    // fill_n(grad, numModules * 2, 0.0);
    // fill_n(binDensity, binTotalNum, 0.0);

    // wrapper::Module Mod = _placement.module(1);
    double mW{0.0}, mH{0.0}, c{0.0};
    double thetaX{0.0}, thetaY{0.0}, dX{0.0}, dY{0.0}, ABSdX{0.0}, ABSdY{0.0}, aX{0.0}, bX{0.0}, aY{0.0}, bY{0.0};
    for (int a = 0; a < binCut; ++a)
    {
        for (int b = 0; b < binCut; ++b)
        {
            for (unsigned i = 0; i < numModules; ++i)
            {
                mW = _placement.module(i).width();
                mH = _placement.module(i).height(); // seg fault ??
                if (!_placement.module(i).isFixed())
                {
                    c = _placement.module(i).area() / binArea;

                    dX = _placement.module(i).centerX() - (((double)a + 0.5) * binW + _placement.boundryLeft());
                    ABSdX = abs(dX);
                    dY = _placement.module(i).centerY() - (((double)b + 0.5) * binH + _placement.boundryBottom());
                    ABSdY = abs(dY);
                    aX = 4.0 / ((mW + 2 * binW) * (mW + 4 * binW));
                    bX = 4.0 / (binW * (mW + 4 * binW));
                    aY = 4.0 / ((mH + 2 * binH) * (mH + 4 * binH));
                    bY = 4.0 / (binH * (mH + 4 * binH));

                    thetaX = (ABSdX <= mW * 0.5 + binW) ? (1 - aX * ABSdX * ABSdX) : (ABSdX <= mW * 0.5 + binW * 2) ? (bX * pow(ABSdX - 2 * mW - 2 * binW, 2))
                                                                                                                    : 0;
                    thetaY = (ABSdY <= mH * 0.5 + binH) ? (1 - aY * ABSdY * ABSdY) : (ABSdY <= mH * 0.5 + binH * 2) ? (bY * pow(ABSdY - 2 * mH - 2 * binH, 2))
                                                                                                                    : 0;

                    binDensity[a + binCut * b] += c * thetaX * thetaY;

                    if (!_placement.module(i).isFixed())
                    {
                        double signX{(dX >= 0) ? 1.0 : -1.0}, signY{(dY >= 0) ? 1.0 : -1.0};
                        grad[2 * i] += (ABSdX <= mW * 0.5 + binW) ? (-2 * signX * c * aX * ABSdX * thetaY) : (ABSdX <= mW * 0.5 + binW * 2) ? (2 * c * bX * signX * (ABSdX - 2 * binW - 2 * mW) * thetaY) : 0;
                        grad[2 * i + 1] += (ABSdY <= mH * 0.5 + binH) ? (-2 * c * signY * aY * ABSdY * thetaX) : (ABSdY <= mH * 0.5 + binH * 2) ? (2 * c * bY * signY * (ABSdY - 2 * binH - 2 * mH) * thetaX) : 0;
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
    // for(int i = 0; i< _placement.numModules(); ++i){
    //     cout << "origin x: " << x[2*i];
    //     cout << ", origin y: " << x[2*i+1] << " / ";
    //     cout << "updated x: " << x[2*i] - g[2*i];
    //     cout << ", updated y: " << x[2*i+1] - g[2*i+1] << endl;
    // }
    // cout << "Done: evaluateFG\n";
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    // f = 3 * x[0] * x[0] + 2 * x[0] * x[1] + 2 * x[1] * x[1] + 7; // objective function

    // cout << "Enter: evaluateF()\n";
    f = 0.0;
    // memset(binDensity, 0.0, sizeof(double)*binTotalNum);
    // fill(binDensity.begin(), binDensity.end(), 0.0);
    // for(int i = 0; i < BINCUT; i++) binDensity[i]=0;
    // fill_n(binDensity, binTotalNum, 0.0);
    // fill(binDensity.begin(), binDensity.end(), 0.0);
    vector<double> binDensity(binTotalNum, 0.0);

    double mW{0.0}, mH{0.0}, c{0.0};
    double thetaX{0.0}, thetaY{0.0}, dX{0.0}, dY{0.0}, ABSdX{0.0}, ABSdY{0.0}, aX{0.0}, bX{0.0}, aY{0.0}, bY{0.0};
    for (int a = 0; a < binCut; ++a)
    {
        for (int b = 0; b < binCut; ++b)
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
                    bX = 4.0 / (binW * (mW + 4 * binW));
                    aY = 4.0 / ((mH + 2 * binH) * (mH + 4 * binH));
                    bY = 4.0 / (binH * (mH + 4 * binH));

                    thetaX = (ABSdX <= mW * 0.5 + binW) ? (1 - aX * ABSdX * ABSdX) : (ABSdX <= mW * 0.5 + binW * 2) ? (bX * pow(ABSdX - 2 * mW - 2 * binW, 2))
                                                                                                                    : 0;
                    thetaY = (ABSdY <= mH * 0.5 + binH) ? (1 - aY * ABSdY * ABSdY) : (ABSdY <= mH * 0.5 + binH * 2) ? (bY * pow(ABSdY - 2 * mH - 2 * binH, 2))
                                                                                                                    : 0;

                    binDensity[a + binCut * b] += c * thetaX * thetaY;
                }
            }
            f += lambda * pow(binDensity[a + binCut * b] - avgDensity, 2);
        }
    }
    // cout << "Done: evaluateF()\n";
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