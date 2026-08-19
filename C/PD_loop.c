#include "declareFunctions.h"
#include "quat.h"

/* == Helper Functions == */

void cross3(const float a[3], const float b[3], float out[3]){
    /* Computes the cross product of two 3-vectors */
    
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

/* == Conversion Functions == */

void tau_to_m(const float tau[3], const float B[3], float m[3]) {
    /* From the wanted torque and the magnetic field, returns the wanted magnetic moment. */
    float B_cross_tau[3] = {0.0f};
    cross3(B, tau, B_cross_tau);

    float normB = norm(B, 3, L_2);
    float denom = normB * normB;

    m[0] = B_cross_tau[0] / denom;
    m[1] = B_cross_tau[1] / denom;
    m[2] = B_cross_tau[2] / denom;
}

void error_axis_angle(
    const float q_want[4],
    const float q_current[4],
    float e[4])
{
    /* From the wanted quaternion and the current quaternion, computes the error in axis angle form */ 
    float q_inv[4];
    float q_error[4];

    quat_inv(q_current, q_inv);
    quat_multiply(q_want, q_inv, q_error);

    float angle = 2.0 * acos(q_error[0]);

    float s = sin(angle/2.0);

    float axis[3];

    axis[0] = q_error[1] / s;
    axis[1] = q_error[2] / s;
    axis[2] = q_error[3] / s;

    e[0] = axis[0];
    e[1] = axis[1];
    e[2] = axis[2];
    e[3] = angle;
}

/* == Main Loop == */

void PD_loop(
    const float r_e[4],
    const float r_omega[4],
    const float B[3],
    float m[3]
) 
/* From the error (r_e) in axis-angle form, the angular velocity (r_omega) in axis-angle form, and magnetic field (B),
Returns the magnetic moment wanted in that time step. */
{

    float placeholder = 1.0;

    float Kp[3] = {-placeholder, -placeholder, -placeholder};
    float Kd[3] = {-placeholder, -placeholder, -placeholder};

    float tau[3] = {0,0,0};

    for (int i = 0; i < 3; i++) {

        float Pi = Kp[i] * r_e[3] * r_e[i];
        float Di = Kd[i] * r_omega[3] * r_omega[i];

        tau[i] = Pi + Di;
    }

    tau_to_m(tau, B, m);
}

void get_m(
    const float q_want[4],
    const float q_current[4],
    const float r_omega[4],
    const float B[3],
    float m[3]
) {
    /* Obtains the magnetic moment by conducting a PD_loop on the inputs. */
    float e[4];
    error_axis_angle(q_want, q_current, e);
    PD_loop(e, r_omega, B, m);
}
