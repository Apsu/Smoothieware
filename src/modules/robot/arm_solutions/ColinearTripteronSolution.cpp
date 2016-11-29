#include "ColinearTripteronSolution.h"
#include "ActuatorCoordinates.h"
#include "ConfigValue.h"
#include "checksumm.h"

#include <fastmath.h>

#define arm_angle_checksum            CHECKSUM("arm_angle")
#define alpha_tower_rotation_checksum CHECKSUM("alpha_tower_rotation")
#define beta_tower_rotation_checksum  CHECKSUM("beta_tower_rotation")
#define gamma_tower_rotation_checksum CHECKSUM("gamma_tower_rotation")

ColinearTripteronSolution::ColinearTripteronSolution(Config* config)
{
    // Retrieve arm angle and compute tangent
    a_a = config->value(arm_angle_checksum)->by_default(30.0f)->as_number();
    a_t = tanf(a_a);

    // Retrieve tower rotations
    a_r = config->value(alpha_tower_rotation_checksum)->by_default(0.0f)->as_number();
    b_r = config->value(beta_tower_rotation_checksum)->by_default(120.0f)->as_number();
    g_r = config->value(gamma_tower_rotation_checksum)->by_default(240.0f)->as_number();

    // Compute tower reductions and factor in tangent coefficient
    a_x = sinf(a_r) * a_t; // Can be simplified if alpha rotation is always 0
    a_y = cosf(a_r) * a_t; // Same as above
    b_x = sinf(b_r) * a_t;
    b_y = cosf(b_r) * a_t;
    g_x = sinf(g_r) * a_t;
    g_y = cosf(g_r) * a_t;

    // Compute forward kinematics matrix denominator
    d = a_y*b_x - g_y*b_x - a_x*b_y - a_y*g_x + b_y*g_x + a_x*g_y;
}

void ColinearTripteronSolution::cartesian_to_actuator(const float cartesian_mm[], ActuatorCoordinates &actuator_mm) const {
    /* Inverse Kinematics
     *
     * 3x3 reduction coefficient matrix * 1x3 coordinate matrix
     * [a_x -a_y 1]   [x]
     * [b_x -b_y 1] * [y]
     * [g_x -g_y 1]   [z]
     */

    actuator_mm[ALPHA_STEPPER] = (this->a_x * cartesian_mm[X_AXIS]) - (this->a_y * cartesian_mm[Y_AXIS]) + cartesian_mm[Z_AXIS];
    actuator_mm[BETA_STEPPER ] = (this->b_x * cartesian_mm[X_AXIS]) - (this->b_y * cartesian_mm[Y_AXIS]) + cartesian_mm[Z_AXIS];
    actuator_mm[GAMMA_STEPPER] = (this->g_x * cartesian_mm[X_AXIS]) - (this->g_y * cartesian_mm[Y_AXIS]) + cartesian_mm[Z_AXIS];
}

void ColinearTripteronSolution::actuator_to_cartesian(const ActuatorCoordinates &actuator_mm, float cartesian_mm[]) const {
    /* Forward Kinematics
     *
     * 3x3 inverted reduction coefficient matrix * 1x3 actuator matrix
     * [a_x -a_y 1]      [a]
     * [b_x -b_y 1]^-1 * [b]
     * [g_x -g_y 1]      [g]
     */

    cartesian_mm[X_AXIS] = (actuator_mm[ALPHA_STEPPER]*(g_y-b_y) + actuator_mm[BETA_STEPPER]*(a_y-g_y) + actuator_mm[GAMMA_STEPPER]*(b_y-a_y))/d;
    cartesian_mm[Y_AXIS] = (actuator_mm[ALPHA_STEPPER]*(g_x-b_x) + actuator_mm[BETA_STEPPER]*(a_x-g_x) + actuator_mm[GAMMA_STEPPER]*(b_x-a_x))/d;
    cartesian_mm[Z_AXIS] = (actuator_mm[ALPHA_STEPPER]*(b_y*g_x-b_x*g_y) + actuator_mm[BETA_STEPPER]*(a_x*g_y-a_y*g_x) + actuator_mm[GAMMA_STEPPER]*(a_y*b_x-a_x*b_y))/d;
}
