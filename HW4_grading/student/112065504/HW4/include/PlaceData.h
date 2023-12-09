#pragma once

class PlaceData
{
public:
    PlaceData()
    {}
    PlaceData(wrapper::Placement &placement, int _binSize, double _gamma)
    {
        binSize = _binSize;
        gamma = _gamma;
        
        double ttlArea = 0;
        for(unsigned mID = 0; mID < placement.numModules(); ++mID)
        {
            ttlArea += placement.module(mID).area();
        }
        maxBinArea = ttlArea / (binSize * binSize);
    }
    int binSize;
    double gamma;

    double maxBinArea{0.0};
};