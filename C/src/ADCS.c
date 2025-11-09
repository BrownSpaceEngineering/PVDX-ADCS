#include "ADCS.h"
#include "UKF.h"
void init(){
    initialized = true;
}
void adcsBody(float* magnetometer_t0, float* magnetometer_t1, float* gyroscope, float* photodiode_inputs, float* tle, int temp, int time, int us_since_last_iter, float* output_currents){
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
        if(!finishedInitialDetumbling){
            float moments[3];
            bool finished_detumbling = bdot_detumbling(magnetometer_t0, magnetometer_t1, moments);
            if(finished_detumbling){
                adcsMode = POINTING;
            }
            moments_to_current(moments, output_currents);
        }
        else{
            float sun_vector[3];
            bool inSun = simulatedSunVector(photodiode_inputs, sun_vector);
            adcsMode = determine_mode(magnetometer_t0, gyroscope, inSun);
            if(adcsMode == DETUMBLING){
                float moments[3];
                bdot_detumbling(magnetometer_t0, magnetometer_t1, moments);
                moments_to_current(moments, output_currents);
            }
            else{
                float moments[3];
                pointing(magnetometer_t0, gyroscope, sun_vector, inSun, moments);
                moments_to_current(moments, output_currents);
            }
            
        }
    }
    
}