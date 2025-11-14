#ifndef POINTING_H
#define POINTING_H

#include <stdbool.h>
static float quaternionState[4];
static float quaternionGoal[4];
static float tleState[6];
static float eciState[6];
static float errorState[6];
static float* kf_covariance[6];

void pointing(float* magnetometer, float* gyroscope, float* sun_vector, bool inSun, float* moments);
void setGoal(float* goal);
void setTLE(float* tle);
bool ukfCollapse();
#endif