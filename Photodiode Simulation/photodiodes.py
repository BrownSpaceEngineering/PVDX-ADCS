import numpy as np
import random 
max_voltage = 1.8
def photodiodes(photodiode_vectors, photodiode_measurements, max_voltage):
    #photodiode_vectors is the list of direction vectors. should all be unit vectors
    #photodiode measurements are the raw voltage readings
    #max voltage is the max photodiode voltage reading
    angles = np.array(photodiode_vectors)
    msmts = np.array(photodiode_measurements) / max_voltage
    return (np.linalg.inv(angles.T@angles))@angles.T@msmts
    
    
#usage
#true sun vector
true_vector = np.array([1,2,3])
true_vector = true_vector / np.linalg.norm(true_vector)

#generates some random vectors
vectors = []
msmts = []
for i in range(6):
    new_vec = np.array([random.random() for _ in range(3)])
    new_vec = new_vec / np.linalg.norm(new_vec)
    vectors.append(new_vec)
    proportion = new_vec@true_vector + np.random.normal(loc = 0, scale = 0.05)
    msmts.append(max_voltage * proportion)

measured_vector = photodiodes(vectors, msmts, max_voltage)
measured_vector = measured_vector / np.linalg.norm(measured_vector)
diff = np.clip(true_vector @ measured_vector, 0, 1)

print("Error: " + str(np.rad2deg(np.acos(diff))) + " degrees")