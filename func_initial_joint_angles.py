import numpy as np
from scipy.optimize import fsolve
from utility_functions.rotation_matrix import *
from utility_functions import *
from func_create_wing_segment import *

def func_initial_joint_angles(theta_1, wing_conformation):

    x0 = np.zeros((7, 1)) # Joints: [1, 2, 3, 5, 6, 8, 9]
    x0[0] = theta_1

    x0[[1,3]] = fsolve(solve_pos4, [0, 0], args=(theta_1, wing_conformation.flatten()), xtol=1e-8).reshape(2,1)
    x0[[2,5]] = fsolve(solve_pos7, [0, 0], args=(x0[[0,1]], wing_conformation.flatten()), xtol=1e-8).reshape(2,1)
    x0[[4,6]] = fsolve(solve_pos10, [0, 0], args=(x0[[3,5]], wing_conformation.flatten()), xtol=1e-8).reshape(2,1)
    x0[4] = 0.2806

    return x0

def solve_pos4(x_solve, theta_1, wing_conformation):
    # theta_1, wing_conformation = arguments[0], arguments[1]
    L1 = wing_conformation[0]
    L2a = wing_conformation[1]
    L3b = wing_conformation[5]
    L3c = wing_conformation[6]
    P1 = wing_conformation[17:19]
    P5 = wing_conformation[19:21]
    alpha_3 = wing_conformation[7]

    # Output variables
    theta_2 = x_solve[0]
    theta_5 = x_solve[1]

    # Constraint function value
    P4_a = P1 + np.dot(rot_2d(theta_1), np.array([L1, 0])) + np.dot(rot_2d(theta_2), np.array([L2a, 0]))
    P4_b = P5 + np.dot(rot_2d(theta_5 - alpha_3), np.array([-L3b, L3c]))

    f = P4_a - P4_b
    return f

def solve_pos7(x_solve, x_in, wing_conformation):
    # Conformation
    L1 = wing_conformation[0]
    L2b = wing_conformation[2]
    L2c = wing_conformation[3]
    L4 = wing_conformation[8]
    L5a = wing_conformation[9]
    L5b = wing_conformation[10]
    P1 = wing_conformation[17:19]
    P8= wing_conformation[21:23]

    # Input variables
    theta_1 = x_in[0]
    theta_2 = x_in[1]

    # Output variables
    theta_3 = x_solve[0]
    theta_8 = x_solve[1]

    # Constraint function value
    P7_a = P1 + np.dot(rot_2d(theta_1), np.array([L1, 0])) + np.dot(rot_2d(theta_2), np.array([L2b, L2c])) + np.dot(rot_2d(theta_3), np.array([L4, 0]))
    P7_b = P8 + np.dot(rot_2d(theta_8), np.array([-L5a, L5b]))
    f = P7_a - P7_b
    return f

def solve_pos10(x_solve, x_in, wing_conformation):
    # Conformation
    L3a = wing_conformation[4]
    L3c = wing_conformation[6]
    alpha_3 = wing_conformation[7]
    L5c = wing_conformation[11]
    L5d = wing_conformation[12]
    L6 = wing_conformation[13]
    L7b = wing_conformation[15]
    L7c = wing_conformation[16]
    P5 = wing_conformation[19:21]
    P8= wing_conformation[21:23]

    # Input variables
    theta_5 = x_in[0]
    theta_8 = x_in[1]

    # Output variables
    theta_6 = x_solve[0]
    theta_9 = x_solve[1]

    # Constraint function value
    P10_a = P8 + np.dot(rot_2d(theta_8), np.array([L5c, L5d])) + np.dot(rot_2d(theta_9), np.array([L6, 0]))
    P10_b = P5 + np.dot(rot_2d(theta_5), np.array([L3a + L3c * np.sin(alpha_3), L3c * np.cos(alpha_3)]) + np.dot(rot_2d(theta_6), np.array([L7b, L7c])))

    f = P10_a - P10_b

    return f

# if __name__ == '__main__':
#     # Linkage length parameters (meter)
#     L1  = 7.5/1000
#     L2a = 35/1000
#     L2b = 5/1000 * np.cos(np.deg2rad(35))
#     L2c = 5/1000 * np.sin(np.deg2rad(35))
#     L3a = 50/1000
#     L3b = 6.71/1000
#     L3c = 6/1000
#     L4  = 36/1000
#     L5a = 11.18/1000
#     L5b = 10/1000
#     L5c = 78.78/1000
#     L5d = -3.89/1000
#     L6  = 36/1000
#     L7a = 90/1000
#     L7b = 35/1000
#     L7c = 15/1000

#     # Stationary joint positions (x,y) (meter)
#     P1  = (np.array([-10, -9.8]) / 1000).reshape(2,1)
#     P5  = (22 * np.array([np.cos(np.deg2rad(20)), np.sin(np.deg2rad(20))]) / 1000).reshape(2,1)
#     P8  = (22 * np.array([np.sin(np.deg2rad(25)), np.cos(np.deg2rad(25))]) / 1000).reshape(2,1)

#     # Wing origin offset about x-axis relative to the body COM
#     offset_x  = np.array([-42/1000]).reshape(1,1) # wing com offset in x-axis (meter) -15 old value

#     # Other parameters
#     alpha_3 = np.deg2rad(-15) # humerus link bend angle

#     # For symbolic generated function
#     L = np.array([L1, L2a, L2b, L2c, L3a, L3b, L3c, alpha_3, L4, 
#                 L5a, L5b, L5c, L5d, L6, L7a, L7b, L7c]).reshape(17,1)
#     wing_conformation = np.concatenate([L, P1, P5, P8, offset_x], axis=0)

#     print(func_initial_joint_angles(np.deg2rad(-90), wing_conformation.flatten()))