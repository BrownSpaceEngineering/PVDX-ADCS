#ifndef ADCS_H
#define ADCS_H

#include <stdbool.h>

enum ADCS_MODE{
    DETUMBLING,
    POINTING,
    FAILURE
};

static enum ADCS_MODE adcsMode = DETUMBLING;
static bool finishedInitialDetumbling = false;
static bool initialized = false;
static float magnetometer_t1[3];

void init();
/**
 * At timestep 0, call this initializer instead of the ADCSBody
 */

void adcsBody(float* magnetometer_t0, float* magnetometer_t1, float* gyroscope, bool inSun, float* sun_vector, float* tle, int temperature, int abs_time, int us_since_last_iter, float* output_currents);
/**
 * This is currently a very basic sketch of ADCS code in C, made as a skeleton to build off of in the future.
 * Not everything under the hood has been included, but the function header should be finalized
 * 
 * Magnetometer_t0: Most recent 3 readings
 * Magnetometer_t1: 3 readings from timesteps - 1 ago. 
 * 
 * Gyroscope: Most recent 3 readings
 * Photodiode_inputs: Most recent readings from all 22
 * TLE: Should use a default value if the value wasn't just uplinked. Otherwise, use uplinked values
 * temperature: in celsius, read from SAMD51
 * time: unix time in seconds
 * us_since_last_iter: microseconds since the last iteration
 * outputCurrents: a float pointer that we can dump output currents for 3 magnetorquers inside of
 * 
 * magnetosphere 
 */
bool failureCheck();

bool simulatedSunVector(float* photodiode_inputs, float* sun_vector);

void bdot_detumbling(float* magnetometer_t0, float* magnetometer_t1, float dt, float* moments);

enum ADCS_MODE determine_mode(float* magnetometer, float* gyroscope, bool inSun);

void moments_to_current(float* moments, float* output_currents);
#endif /* ADCS_H */

