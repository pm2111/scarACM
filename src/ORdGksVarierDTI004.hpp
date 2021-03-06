#ifndef ORdGKSVARIER_HPP_
#define ORdGKSVARIER_HPP_

#include "HeartGeometryInformation.hpp"

class ORdGksVarierDTI004
{
private:
    HeartGeometryInformation<3>* mpHeartGeomInfo;
    bool mApexbase;
    bool mTransmur;
    bool mInterven;
    int mNormaliser;

public:
    /*
     * Constructor. Takes a heart geometry info object and which gradients we're using.
     */
    ORdGksVarierDTI004(HeartGeometryInformation<3>* pHeartGeomInfo, bool ApexBase, bool Transmural, bool Interven) :
    mpHeartGeomInfo(pHeartGeomInfo),
    mApexbase(ApexBase),
    mTransmur(Transmural),
    mInterven(Interven)
    {
        mNormaliser = (int)mApexbase+(int)mTransmur+(int)mInterven;
    }
    /*
     * Method to return the scaling factor.
     */
    double GetScalingFactor(const Node<3>& rNode);
};

#endif //TT06GKSVARIER_HPP_
