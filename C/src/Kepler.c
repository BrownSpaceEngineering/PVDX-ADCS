#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void propogateOrbitalElements(float semi_major_axis, float eccentricity, float inclination, float ascending_node, float periapsis, float true_anomaly, float time_delta, float *array) {
//Mass of earth is extremely larger than the sattelite, so:

//gravitational  
float u = 3.986004418e14; //[m^3/s^2]
float a = semi_major_axis; //[m]
float e = eccentricity; 
float f = true_anomaly * M_PI / 180.0; //[radians]
float dt = time_delta; //[seconds]


//convert true anomaly to eccentric anamoly
float E =  2 * atan2( sqrt(1 - e) * sin(f/2), sqrt(1 + e) * cos(f/2));

//Get current true anomaly from Kepler's equation
float M = E - e * sin(E);

//define mean motion
float n = sqrt(u/pow(a, 3.0));

//Find mean anomaly at T + dt
float M_new = M + n * dt;

//Use newton's method to iterate until we find the new_true_anomly. 
float E_current = M_new;
float tolerance = 1e-12;
int MAXIMUM_ITERATIONS = 10;

for (int i = 1; i < MAXIMUM_ITERATIONS; ++i) {
    //Calculate the derivative of kepler's equation
    float E_change = ((E_current - e * sin(E_current) - M_new)/(1 - e * cos(E_current)));
  
    //calculate next iteration
    E_current = E_current - E_change;

    //if tolerance is small enough, break out of the loop
    if(fabs(E_change) < tolerance)
        break;
}
array[0] = a;
array[1] = e;
array[2] = inclination;
array[3] = ascending_node;
array[4] = periapsis;

//updated true anomaly from our updated eccentric anomaly
float new_true_anomaly = 2 * atan2( sqrt(1 + e) * sin(E_current/2), sqrt(1 - e) * cos(E_current/2));
array[5] = new_true_anomaly * 180.0 / M_PI; //convert back to degrees

}

void testKepler(float axis, float e0, float inc0, float RAAN0, float argp0, float nu0) {
    FILE *fp = fopen("output.txt", "w");
    for (int j =0; j < 120; ++j) {
        float array[6];
        propogateOrbitalElements(axis, e0, inc0, RAAN0, argp0, nu0, j * 60, array);
        float anomaly = array[5];
        fprintf(fp, "%.5f\n", anomaly);
    }  
}
