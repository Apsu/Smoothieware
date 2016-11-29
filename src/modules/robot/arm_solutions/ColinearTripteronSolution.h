#pragma once

#include "libs/Module.h"
#include "libs/Kernel.h"
#include "BaseSolution.h"
#include "libs/nuts_bolts.h"

#include "libs/Config.h"

class ColinearTripteronSolution : public BaseSolution {
    public:
        ColinearTripteronSolution(Config*);
        void cartesian_to_actuator(const float[], ActuatorCoordinates & ) const override;
        void actuator_to_cartesian(const ActuatorCoordinates &, float[] ) const override;

    private:
        float a_a;
        float a_t;
        float a_r;
        float a_x;
        float a_y;
        float b_r;
        float b_x;
        float b_y;
        float g_r;
        float g_x;
        float g_y;
        float d;
};
