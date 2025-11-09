import numpy as np
import pyquaternion
from pyquaternion import Quaternion
from simulation_tools import *
from numpy.linalg import cholesky
from scipy.linalg import sqrtm
from filterpy.kalman import UnscentedKalmanFilter, MerweScaledSigmaPoints
import math
from numpy.typing import NDArray
from typing import List, Tuple
from copy import deepcopy
import time
import matplotlib.pyplot as plt

#ALL CREDITS TO https://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf FOR THE EQUATIONS
'''
Overall Idea of a UKF: 
- Like the MEKF, rather than using a 7-element state [quaternion, gyro bias], we use a 6-element state [quaternion error vector, gyro bias] and assume they are both small
- Given a timestep, a measurement, and our last state, we take some informed samples around our estimated state to guess where we are(what our current state is)
'''

'''
Why BDot is probably not great -- 
- As your dt grows smaller and smaller, the noise to signal ratio in the finite-differences 
derivative gets higher and higher
- At the same time, as your dt grows larger and larger, error quaternions grow larger and larger
    - intuitively, a lower dt means more updates means better results
- this dichotemy is worsened by the fact that, technically, there are always 2 solutions to a large error
    - 1. your gyro bias is huge OR 2. your quaternion is off
    - while this dichotemy is always present, it becomes far more exacerbated when both gyro bias and quaternion are 
    used in propagation AND measurement calculations(rather than only quaternion being used in both)
'''

global_Q = np.eye(6)*0.0001
global_R = np.eye(3)*0.001
#global_R = np.eye(6)*0.1
dt = 0.2
timeframe = datetime.now()
kepler_posn = np.array([7e6, math.radians(0.05), math.radians(50), math.radians(311.6218), math.radians(199.2431), math.radians(48.4420)])
xyz_posn = kep_to_cart(kepler_posn)[:3]
magnetosphere_measurements = []
magnetometer_readings = []

gyro_measurement = [0,0,0]

magnetosphere_measurements = []

def rotate(vec, q):
    v = Quaternion(scalar = 0, vector = vec)
    return (q * v * q.inverse).vector
def quaternion_to_rotation(x : Quaternion):
    '''
    Function that converts a quaternion into a rotation vector
    '''
    quat = x.elements
    if quat[0] < 0:
        quat = -quat
    if(quat[0] >= 1):
        return np.zeros(3)
    theta = 2 * math.acos(quat[0])
    if(theta == 0):
        return np.zeros(3)
    vec = theta * quat[1:] / np.sqrt(1-quat[0]**2)
    return vec
def rotation_to_quat(vec):
    '''
    Function that converts a rotation vector into a quaternion
    '''
    if(np.linalg.norm(vec) == 0):
        return Quaternion()
    axis = vec / np.linalg.norm(vec)
    angle = np.linalg.norm(vec)
    return Quaternion(angle = angle, axis = axis)
def get_error_vectors(Y: List[Quaternion], x : Quaternion):
    '''
    Function that gets the error vectors between all quaternions in Y and the quaternion X
    '''
    error_vectors = []
    for quat in Y:
        error_quaternion = quat * x.inverse
        error_vectors.append(quaternion_to_rotation(error_quaternion))
    return error_vectors
def gradient_descent(Y : List[Quaternion], x : Quaternion):
    '''
    Function that performs gradient descent to find the average quaternion, as quaternions are not additive
    '''
    for i in range(100):
        error_vectors = get_error_vectors(Y, x)
        ave = np.mean(error_vectors, axis=0)   # formula 54
        if(np.linalg.norm(ave) < 0.001):
            break
        average_error_quat = rotation_to_quat(ave)
        x = average_error_quat * x    # formula 55
        
    return x, error_vectors


def get_weights(lam, x, alpha, beta):
    '''
    Derived from Van Der Merwe Scaled(Weighted) Sigma Points
    '''
    n = x.shape[0]
    c = .5 / (n + lam)
    covariance_weights = np.full(2*n + 1, c)
    mean_weights = np.full(2*n + 1, c)
    covariance_weights[0] = lam / (n + lam) + (1 - alpha**2 + beta)
    mean_weights[0] = lam/ (n + lam)
    return covariance_weights, mean_weights
def get_sigma_points(lam, x, P):
    '''
    Derived from Van Der Merwe Scaled(Weighted) Sigma Points
    '''
    n = x.shape[0]
    P = (P + P.T) / 2
    P += np.eye(n) * 1e-10
    U =cholesky((lam + n)*P).T
    sigmas = [x]
    for k in range(1, n + 1, 1):#from 1 -> n
        idx = k - 1
        sigmas.append(np.add(x, U[idx]))
    for k in range(n + 1, 2*n + 1, 1):#from n + 1 -> 2n + 1
        idx = k - (n + 1)
        sigmas.append(np.subtract(x, U[idx]))
    return np.array(sigmas)
def calculate_lambda(alpha, x):
    '''
    Derived from Van Der Merwe Scaled(Weighted) Sigma Points
    '''
    n = x.shape[0]
    kappa = 3-n
    return alpha**2 * (n + kappa) - n

def error_sigmas_to_quat_sigmas(error_sigmas, rotation : Quaternion):
    '''
    Converts the 6-element sigma point to a 7-element sigma point through use of the current rotation and attitude error vector
    '''
    quat_sigmas = []
    for sigma in error_sigmas:
        quat_rot = rotation_to_quat(sigma[:3])
        new_quat_rot = rotation * quat_rot
        new_quat_rot = new_quat_rot.normalised
        quat_sigmas.append(np.concatenate([new_quat_rot.elements, sigma[3:]]))
    return np.array(quat_sigmas)

def propagate_quat_sigmas(quat_sigmas):
    '''
    Propagates the 7-element sigma points forward a timestep so that we can guess
    '''
    global gyro_measurement
    new_sigmas = []
    for sigma in quat_sigmas:
        quaternion = Quaternion(sigma[:4])
        w = -sigma[4:7] + gyro_measurement
        change_in_rotation = rotation_to_quat(w * dt)
        new_quaternion = quaternion * change_in_rotation
        
        new_quaternion = new_quaternion.normalised
        new_sigmas.append(np.concatenate([new_quaternion.elements, sigma[4:7]]))
    return np.array(new_sigmas)

def get_measurements(quat_sigmas):
    '''
    Gets the measurements expected of the 7-element sigma points using
    b and bdot
    '''
    measurements = []
    for sigma in quat_sigmas:
        quat = Quaternion(sigma[:4])#rotation from body to ECI
        ECI_to_body = quat.inverse
        
        expected_reading_one = rotate(np.array([1, 0, 0]), ECI_to_body)
        expected_reading_two = rotate(np.array([0, 1, 0]), ECI_to_body)
        
        #expected_reading_two = np.zeros(3) - np.cross(w, magnetometer_readings[-1])
        #measurements.append(np.concatenate((expected_reading_one, expected_reading_two)))
        measurements.append(expected_reading_one)
    return np.array(measurements)

def ensure_positive_definite(P, epsilon=1e-9):
    """Ensure matrix is positive definite"""
    P = (P + P.T) / 2  # Symmetrize
    
    # Check eigenvalues
    eigvals = np.linalg.eigvalsh(P)
    if np.min(eigvals) < epsilon:
        P += np.eye(P.shape[0]) * (epsilon - np.min(eigvals) + epsilon)
    
    return P

def iterate(error_state, rotation : Quaternion, P, obs):
    '''
    Iterates forward given our last 6-element state & covariance, our last rotation, and a measurement
    '''
    global global_Q, global_R
    n = error_state.shape[0]
    alpha = 1e-4
    beta = 2
    P = ensure_positive_definite(P)
    lam = calculate_lambda(alpha, error_state)
    sigmas = get_sigma_points(lam, error_state, P + global_Q)
    quat_sigmas = error_sigmas_to_quat_sigmas(sigmas, rotation)
    propagated_quat_sigmas = propagate_quat_sigmas(quat_sigmas)
    measurements = get_measurements(propagated_quat_sigmas)
    covariance_weights, mean_weights = get_weights(lam, error_state, alpha, beta)

    quaternions_of_propagated_sigmas = []
    for y in propagated_quat_sigmas:
        quaternions_of_propagated_sigmas.append(Quaternion(y[:4]))
    average_quaternion, propagated_error_vectors = gradient_descent(quaternions_of_propagated_sigmas, rotation)
    propagated_errors = np.hstack([propagated_error_vectors, propagated_quat_sigmas[:, 4:]])

    
    mean_error = np.zeros(n)
    for column in range(n):#for each state element
        var_mean = 0
        for row in range(2*n + 1):#each row is the corresponding element for 1 sigma point
            var_mean += propagated_errors[row][column] * mean_weights[row]
        mean_error[column] = var_mean
    mean_measurement = np.zeros(measurements.shape[1])
    for column in range(measurements.shape[1]):#for each measurement
        var_mean = 0
        for row in range(2*n + 1):#each row is the corresponding measurement for 1 sigma point
            var_mean += measurements[row][column] * mean_weights[row]
        mean_measurement[column] = var_mean

    P_hat = np.zeros((n, n))
    
    for i, error in enumerate(propagated_errors):
        err = error - mean_error
        P_hat += covariance_weights[i] * np.outer(err, err)
    
    P_xz = np.zeros((n, measurements.shape[1]))
    for i, error in enumerate(propagated_errors):
        err = error - mean_error
        msmt_err = measurements[i] - mean_measurement
        P_xz += covariance_weights[i] * np.outer(err,msmt_err)
    
    P_zz = np.zeros((measurements.shape[1], measurements.shape[1]))
    for i, msmt in enumerate(measurements):
        msmt_err = msmt - mean_measurement
        P_zz += covariance_weights[i] * np.outer(msmt_err, msmt_err)
    P_vv = P_zz + global_R
    k = P_xz @ np.linalg.inv(P_vv)
    x_hat = mean_error + k@(obs - mean_measurement)

    I_KH = np.eye(n) - k @ np.linalg.inv(P_vv) @ P_xz.T
    P = I_KH @ P_hat @ I_KH.T + k @ global_R @ k.T
    P = ensure_positive_definite(P)

    x_hat_rot = rotation_to_quat(x_hat[:3])
    return x_hat, x_hat_rot * average_quaternion, P


if __name__ == '__main__':
    true_rot = Quaternion(np.random.normal(loc = 0.5, scale = 0.5, size = 4)).normalised
    rot = true_rot.elements + np.random.normal(loc = 0, scale = 0.05, size = 4)
    rot = Quaternion(rot).normalised
    P = np.eye(6) * 0.001
    state = np.zeros(6)

    true_angular_velocity = np.array([0.05, 0.01, 0])#30 degrees per second, which is too fast for us but is good for simulation purposes
    gyro_bias = np.array([0.00005, 0.00005, -0.0001])

    #with a simulated u vector
    rotation_quaternion = Quaternion(scalar = 0, vector = true_angular_velocity)

    magnetometer_readings.append(rotate(np.array([1,0,0]) + np.random.normal(loc = gyro_bias, scale = 0.002, size = 3), true_rot.inverse))
    
    
    diffs = []
    total_bias = np.zeros(3)
    for i in range(15000):
        #with gyro
        gyro_measurement = true_angular_velocity + np.random.normal(loc = gyro_bias, scale = 0.001)
        
        #without gyro
        #gyro_measurement = [0,0,0]
        
        change_in_true_rotation = rotation_to_quat(true_angular_velocity * dt)
        true_rot = true_rot * change_in_true_rotation#from body to reference
        ref_to_body = true_rot.inverse

        given_reading_one = rotate(np.array([1,0,0]) + np.random.normal(loc = 0, scale = 0.01, size = 3), ref_to_body)
        given_reading_one = given_reading_one / np.linalg.norm(given_reading_one)

        given_reading_two = rotate(np.array([0,1,0]) + np.random.normal(loc = 0, scale = 0.01, size = 3), ref_to_body)
        given_reading_two = given_reading_two / np.linalg.norm(given_reading_two)

        magnetometer_readings.append(given_reading_one)

        #print(derivative_reading_one, -np.cross(true_angular_velocity, given_reading_one))


        measurement = given_reading_one
        #measurement = np.concatenate([given_reading_one, given_reading_two])
        diffs.append(np.linalg.norm(np.array(quat_diff(true_rot, rot))))
        if(i % 100 == 0):
            print(true_rot, rot)
            print(quat_diff(true_rot, rot))
        try:
            state, rot, P = iterate(state, rot, P, measurement)
            state[:6] = np.zeros(6)#we need to reset our error vector here, as we already tacked on the error vector to the rotation at the end of the last state! 
            #total_bias += state[3:6]
        except Exception as e:
            print(np.linalg.eigvals(P))
            print(e)
            break
    #print("AVG BIAS: ", total_bias / 10000)
    # plt.title("Pointing Error Over Time Under No Spin -- One Vector")
    # plt.scatter(np.arange(10000) * dt, diffs, s = 3)
    # plt.xlabel("Seconds")
    # plt.ylabel("Magnitude of Degrees(error)")
    # plt.show()
    # print(rot)
    # print(true_rot)
    # print(quat_diff(rot, true_rot))