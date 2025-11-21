#include "../include/ADCS.h"
#include "../include/UKF.h"
#include <string.h>
#include <stdio.h>
void init(){
    initialized = true;
}
void moments_to_current(float* moments, float* currents){
    memcpy(currents, moments, 3 * sizeof(float));
}
bool failureCheck(){
    return false;
}
void bdot_detumbling(float* magnetometer_t0, float* magnetometer_t1, float dt, float* moments){
    //code copied from bDot.m
    
    float derivative[3];
    float k = 1;
    //may change from finite diffs to something else depending on magnetometer noise
    derivative[0] = (magnetometer_t0[0] - magnetometer_t1[0])/dt;
    derivative[1] = (magnetometer_t0[1] - magnetometer_t1[1])/dt;
    derivative[2] = (magnetometer_t0[2] - magnetometer_t1[2])/dt;

    moments[0] = -k * derivative[0];
    moments[1] = -k * derivative[1];
    moments[2] = -k * derivative[2];
}

enum ADCS_MODE determine_mode(float* magnetometer, float* gyroscope, bool inSun){
    return DETUMBLING;
}

bool simulatedSunVector(float* photodiode_inputs, float* sun_vector){
    //S6 on the design doc
    return false;
}

void adcsBody(float* magnetometer_t0, float* magnetometer_t1, float* gyroscope, bool inSun, float* sun_vector, float* tle, int temperature, int abs_time, int us_since_last_iter, float* output_currents){
    //S4 on the design docs
    if(!initialized){
        init();
        return;
    }
    else if(adcsMode != FAILURE && failureCheck()){
        adcsMode = FAILURE;
        return;
    }
    else if(adcsMode == FAILURE){
        return;
    }
    else{
        adcsMode = determine_mode(magnetometer_t0, gyroscope, inSun);
        if(adcsMode == DETUMBLING){
            float moments[3];
            float dt = (float)us_since_last_iter * 1e-6;
            bdot_detumbling(magnetometer_t0, magnetometer_t1, dt, moments);
            moments_to_current(moments, output_currents);
        }
        else{
            float moments[3];
            pointing(magnetometer_t0, gyroscope, sun_vector, inSun, moments);
            moments_to_current(moments, output_currents);
        }
    }
    
}

int main(){
    float magnetometer_0[3] = {1, 1, 1};
    float magnetometer_1[3] = {-1, -1, -1};
    float dt = 1;
    float moments[3];
    bdot_detumbling(magnetometer_0, magnetometer_1, dt, moments);

    //returns 1
    return (moments[0] == -2 && moments[1] == -2 && moments[2] == -2);
}