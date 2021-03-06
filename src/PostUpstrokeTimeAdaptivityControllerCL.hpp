#ifndef POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_
#define POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_

#include "AbstractTimeAdaptivityController.hpp"

#include "BidomainProblem.hpp" ///\todo: template over problem?

class PostUpstrokeTimeAdaptivityControllerCL : public AbstractTimeAdaptivityController
{
private:
    double mAllDepolarisedTime;
    bool mHaveSwitchedTimestep;
    double mSmallDt;
    double mLargeDt;
    BidomainProblem<3>& mrProblem; // For tissue, for cells
    double mPacingFrequency; /**< Time period at which to resume small steps. Default is 0 (never reset) */

public:
    /*
     *  Constructor.
     */
    PostUpstrokeTimeAdaptivityControllerCL(double smallDt, double largeDt, BidomainProblem<3>& rProblem, double pacingFrequency=0.0)
        : AbstractTimeAdaptivityController(smallDt, largeDt),
          mAllDepolarisedTime(-1.0), // Nonsense value to indicate not all depolarised
          mHaveSwitchedTimestep(false),
          mSmallDt(smallDt),
          mLargeDt(largeDt),
          mrProblem(rProblem),
          mPacingFrequency(pacingFrequency)
    {
    }
    /*
     *  The method that returns the time step
     */
    double ComputeTimeStep(double currentTime, Vec currentSolution);
};

#endif //POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_
