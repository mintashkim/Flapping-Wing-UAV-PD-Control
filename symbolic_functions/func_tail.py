import numpy as np
    
def func_tail(in1 = None,in2 = None): 

    Rb_1_1 = in1[13]
    Rb_1_2 = in1[16]
    Rb_1_3 = in1[19]
    Rb_2_1 = in1[14]
    Rb_2_2 = in1[17]
    Rb_2_3 = in1[20]
    Rb_3_1 = in1[15]
    Rb_3_2 = in1[18]
    Rb_3_3 = in1[21]
    dtail_1 = in2[0]
    dtail_2 = in2[1]
    dtail_3 = in2[2]
    w_1 = in1[10]
    w_2 = in1[11]
    w_3 = in1[12]
    xD_1 = in1[7]
    xD_2 = in1[8]
    xD_3 = in1[9]
    x_1 = in1[2]
    x_2 = in1[3]
    x_3 = in1[4]
    pos_aero_tail = np.array([[x_1 + Rb_1_1 * dtail_1 + Rb_1_2 * dtail_2 + Rb_1_3 * dtail_3],[x_2 + Rb_2_1 * dtail_1 + Rb_2_2 * dtail_2 + Rb_2_3 * dtail_3],[x_3 + Rb_3_1 * dtail_1 + Rb_3_2 * dtail_2 + Rb_3_3 * dtail_3]])
    
    #if nargout > 1
    vel_aero_tail = np.array([[xD_1 + dtail_3 * (Rb_1_1 * w_2 - Rb_1_2 * w_1) - dtail_2 * (Rb_1_1 * w_3 - Rb_1_3 * w_1) + dtail_1 * (Rb_1_2 * w_3 - Rb_1_3 * w_2)],[xD_2 + dtail_3 * (Rb_2_1 * w_2 - Rb_2_2 * w_1) - dtail_2 * (Rb_2_1 * w_3 - Rb_2_3 * w_1) + dtail_1 * (Rb_2_2 * w_3 - Rb_2_3 * w_2)],[xD_3 + dtail_3 * (Rb_3_1 * w_2 - Rb_3_2 * w_1) - dtail_2 * (Rb_3_1 * w_3 - Rb_3_3 * w_1) + dtail_1 * (Rb_3_2 * w_3 - Rb_3_3 * w_2)]])
    
    #if nargout > 2
    Ba_tail = np.reshape(np.array([0.0,0.0,1.0,0.0,0.0,- Rb_1_2 * dtail_3 + Rb_1_3 * dtail_2,Rb_1_1 * dtail_3 - Rb_1_3 * dtail_1,- Rb_1_1 * dtail_2 + Rb_1_2 * dtail_1,0.0,0.0,0.0,1.0,0.0,- Rb_2_2 * dtail_3 + Rb_2_3 * dtail_2,Rb_2_1 * dtail_3 - Rb_2_3 * dtail_1,- Rb_2_1 * dtail_2 + Rb_2_2 * dtail_1,0.0,0.0,0.0,0.0,1.0,- Rb_3_2 * dtail_3 + Rb_3_3 * dtail_2,Rb_3_1 * dtail_3 - Rb_3_3 * dtail_1,- Rb_3_1 * dtail_2 + Rb_3_2 * dtail_1]), tuple(np.array([8,3])), order="F")

    return pos_aero_tail, vel_aero_tail, Ba_tail