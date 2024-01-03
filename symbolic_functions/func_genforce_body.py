import numpy as np
    
def func_genforce_body(in1 = None, in2 = None, in3 = None):
    
    Rb_1_1 = in1[13]
    Rb_1_2 = in1[16]
    Rb_1_3 = in1[19]
    Rb_2_1 = in1[14]
    Rb_2_2 = in1[17]
    Rb_2_3 = in1[20]
    Rb_3_1 = in1[15]
    Rb_3_2 = in1[18]
    Rb_3_3 = in1[21]
    dforce_1 = in2[0]
    dforce_2 = in2[1]
    dforce_3 = in2[2]
    ubody_1 = in3[0]
    ubody_2 = in3[1]
    ubody_3 = in3[2]
    output_body_gen_force = np.array([[0.0],[0.0],[ubody_1],[ubody_2],[ubody_3],[- ubody_1 * (Rb_1_2 * dforce_3 - Rb_1_3 * dforce_2) - ubody_2 * (Rb_2_2 * dforce_3 - Rb_2_3 * dforce_2) - ubody_3 * (Rb_3_2 * dforce_3 - Rb_3_3 * dforce_2)],[ubody_1 * (Rb_1_1 * dforce_3 - Rb_1_3 * dforce_1) + ubody_2 * (Rb_2_1 * dforce_3 - Rb_2_3 * dforce_1) + ubody_3 * (Rb_3_1 * dforce_3 - Rb_3_3 * dforce_1)],[- ubody_1 * (Rb_1_1 * dforce_2 - Rb_1_2 * dforce_1) - ubody_2 * (Rb_2_1 * dforce_2 - Rb_2_2 * dforce_1) - ubody_3 * (Rb_3_1 * dforce_2 - Rb_3_2 * dforce_1)]])

    return output_body_gen_force
