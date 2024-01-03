import numpy as np
    
def func_wing_kinematic_vel(in1, vel_theta1, in3):

       L1 = in3[0]
       L4 = in3[8]
       L6 = in3[13]
       L2a = in3[1]
       L3a = in3[4]
       L5a = in3[9]
       L2b = in3[2]
       L3b = in3[5]
       L5b = in3[10]
       L7b = in3[15]
       L2c = in3[3]
       L3c = in3[6]
       L5c = in3[11]
       L7c = in3[16]
       L5d = in3[12]
       alpha_3 = in3[7]
       theta_1 = in1[0]
       theta_2 = in1[1]
       theta_3 = in1[2]
       theta_5 = in1[3]
       theta_6 = in1[4]
       theta_8 = in1[5]
       theta_9 = in1[6]

       mt1 = np.array([- L1 * np.sin(theta_1),L1 * np.cos(theta_1),- L1 * np.sin(theta_1),L1 * np.cos(theta_1),0.0,0.0,1.0,- L2a * np.sin(theta_2),L2a * np.cos(theta_2),- L2c * np.cos(theta_2) - L2b * np.sin(theta_2),L2b * np.cos(theta_2) - L2c * np.sin(theta_2),0.0,0.0,0.0,0.0,0.0,- L4 * np.sin(theta_3),L4 * np.cos(theta_3),0.0,0.0,0.0,L3c * np.cos(alpha_3 - theta_5) + L3b * np.sin(alpha_3 - theta_5),L3b * np.cos(alpha_3 - theta_5) - L3c * np.sin(alpha_3 - theta_5),0.0,0.0,np.sin(theta_5) * (L3a + L3c * np.sin(alpha_3)) + L3c * np.cos(alpha_3) * np.cos(theta_5),- np.cos(theta_5) * (L3a + L3c * np.sin(alpha_3)) + L3c * np.cos(alpha_3) * np.sin(theta_5),0.0,0.0,0.0,0.0,0.0,L7c * np.cos(theta_6) + L7b * np.sin(theta_6),- L7b * np.cos(theta_6) + L7c * np.sin(theta_6),0.0,0.0,0.0,L5b * np.cos(theta_8) - L5a * np.sin(theta_8),L5a * np.cos(theta_8) + L5b * np.sin(theta_8)])
       mt2 = np.array([- L5d * np.cos(theta_8) - L5c * np.sin(theta_8),L5c * np.cos(theta_8) - L5d * np.sin(theta_8),0.0,0.0,0.0,0.0,0.0,- L6 * np.sin(theta_9),L6 * np.cos(theta_9),0.0])
       A_kv = np.concatenate([mt1, mt2], axis=0).reshape(7,7).T

       h_kv = np.array([[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0]])

       u_kv = np.array([[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[vel_theta1]])
       output_ki_vel = np.concatenate((A_kv,(h_kv).reshape(1,7),u_kv.reshape(1,7)), axis=0)

       return output_ki_vel

if __name__ == '__main__':

    # Linkage length parameters (meter)
    L1  = 7.5/1000
    L2a = 35/1000
    L2b = 5/1000 * np.cos(np.deg2rad(35))
    L2c = 5/1000 * np.sin(np.deg2rad(35))
    L3a = 50/1000
    L3b = 6.71/1000
    L3c = 6/1000
    L4  = 36/1000
    L5a = 11.18/1000
    L5b = 10/1000
    L5c = 78.78/1000
    L5d = -3.89/1000
    L6  = 36/1000
    L7a = 90/1000
    L7b = 35/1000
    L7c = 15/1000

    # Stationary joint positions (x,y) (meter)
    P1  = (np.array([-10, -9.8]) / 1000).reshape(2,1)
    P5  = (22 * np.array([np.cos(np.deg2rad(20)), np.sin(np.deg2rad(20))]) / 1000).reshape(2,1)
    P8  = (22 * np.array([np.sin(np.deg2rad(25)), np.cos(np.deg2rad(25))]) / 1000).reshape(2,1)

    # Wing origin offset about x-axis relative to the body COM
    offset_x  = np.array([-42/1000]).reshape(1,1) # wing com offset in x-axis (meter) -15 old value

    # Other parameters
    alpha_3 = np.deg2rad(-15) # humerus link bend angle

    # For symbolic generated function
    L = np.array([L1, L2a, L2b, L2c, L3a, L3b, L3c, alpha_3, L4, 
                L5a, L5b, L5c, L5d, L6, L7a, L7b, L7c]).reshape(17,1)
    wing_conformation = np.concatenate([L, P1, P5, P8, offset_x], axis=0)
    xk = np.array([-1.5708, 0.8912, 1.4493, 0.1980, -0.2806, 0.4953, -0.8275])
    print(func_wing_kinematic_vel(xk, 4.75 * 2 * np.pi, wing_conformation.flatten()))


