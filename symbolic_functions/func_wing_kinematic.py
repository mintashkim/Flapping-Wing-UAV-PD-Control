import numpy as np
    
def func_wing_kinematic(in1 = None,u_theta1 = None,in3 = None): 

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
    thetad_1 = in1[7]
    thetad_2 = in1[8]
    thetad_3 = in1[9]
    thetad_5 = in1[10]
    thetad_6 = in1[11]
    thetad_8 = in1[12]
    thetad_9 = in1[13]
    
    mt1 = np.array([- L1 * np.sin(theta_1),L1 * np.cos(theta_1),- L1 * np.sin(theta_1),L1 * np.cos(theta_1),0.0,0.0,1.0,- L2a * np.sin(theta_2),L2a * np.cos(theta_2),- L2c * np.cos(theta_2) - L2b * np.sin(theta_2),L2b * np.cos(theta_2) - L2c * np.sin(theta_2),0.0,0.0,0.0,0.0,0.0,- L4 * np.sin(theta_3),L4 * np.cos(theta_3),0.0,0.0,0.0,L3c * np.cos(alpha_3 - theta_5) + L3b * np.sin(alpha_3 - theta_5),L3b * np.cos(alpha_3 - theta_5) - L3c * np.sin(alpha_3 - theta_5),0.0,0.0,np.sin(theta_5) * (L3a + L3c * np.sin(alpha_3)) + L3c * np.cos(alpha_3) * np.cos(theta_5),- np.cos(theta_5) * (L3a + L3c * np.sin(alpha_3)) + L3c * np.cos(alpha_3) * np.sin(theta_5),0.0,0.0,0.0,0.0,0.0,L7c * np.cos(theta_6) + L7b * np.sin(theta_6),- L7b * np.cos(theta_6) + L7c * np.sin(theta_6),0.0,0.0,0.0,L5b * np.cos(theta_8) - L5a * np.sin(theta_8),L5a * np.cos(theta_8) + L5b * np.sin(theta_8)])
    mt2 = np.array([- L5d * np.cos(theta_8) - L5c * np.sin(theta_8),L5c * np.cos(theta_8) - L5d * np.sin(theta_8),0.0,0.0,0.0,0.0,0.0,- L6 * np.sin(theta_9),L6 * np.cos(theta_9),0.0])
    A_k = ((np.concatenate((mt1,mt2),axis=0)).reshape(7,7)).T
    
    # if nargout > 1
    mt3 = np.array([- L1 * thetad_1 ** 2 * np.cos(theta_1) - L2a * thetad_2 ** 2 * np.cos(theta_2) - L3b * thetad_5 ** 2 * np.cos(alpha_3 - theta_5) + L3c * thetad_5 ** 2 * np.sin(alpha_3 - theta_5),- L1 * thetad_1 ** 2 * np.sin(theta_1) - L2a * thetad_2 ** 2 * np.sin(theta_2) + L3c * thetad_5 ** 2 * np.cos(alpha_3 - theta_5) + L3b * thetad_5 ** 2 * np.sin(alpha_3 - theta_5),- L1 * thetad_1 ** 2 * np.cos(theta_1) - L4 * thetad_3 ** 2 * np.cos(theta_3) - L5a * thetad_8 ** 2 * np.cos(theta_8) - L2b * thetad_2 ** 2 * np.cos(theta_2) - L5b * thetad_8 ** 2 * np.sin(theta_8) + L2c * thetad_2 ** 2 * np.sin(theta_2),L5b * thetad_8 ** 2 * np.cos(theta_8) - L2c * thetad_2 ** 2 * np.cos(theta_2) - L1 * thetad_1 ** 2 * np.sin(theta_1) - L4 * thetad_3 ** 2 * np.sin(theta_3) - L5a * thetad_8 ** 2 * np.sin(theta_8) - L2b * thetad_2 ** 2 * np.sin(theta_2)])
    mt4 = np.array([thetad_5 ** 2 * np.cos(theta_5) * (L3a + L3c * np.sin(alpha_3)) - L6 * thetad_9 ** 2 * np.cos(theta_9) + L7b * thetad_6 ** 2 * np.cos(theta_6) - L5c * thetad_8 ** 2 * np.cos(theta_8) - L7c * thetad_6 ** 2 * np.sin(theta_6) + L5d * thetad_8 ** 2 * np.sin(theta_8) - L3c * thetad_5 ** 2 * np.cos(alpha_3) * np.sin(theta_5),thetad_5 ** 2 * np.sin(theta_5) * (L3a + L3c * np.sin(alpha_3)) + L7c * thetad_6 ** 2 * np.cos(theta_6) - L5d * thetad_8 ** 2 * np.cos(theta_8) - L6 * thetad_9 ** 2 * np.sin(theta_9) + L7b * thetad_6 ** 2 * np.sin(theta_6) - L5c * thetad_8 ** 2 * np.sin(theta_8) + L3c * thetad_5 ** 2 * np.cos(alpha_3) * np.cos(theta_5),0.0])
    h_k = (np.concatenate((mt3,mt4),axis=0)).reshape(7,1)
    
    # if nargout > 2
    u_k = np.array([[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[u_theta1]])
    
    #size([A_k;h_k';u_k'])=(9,7)
    return A_k, h_k, u_k
    