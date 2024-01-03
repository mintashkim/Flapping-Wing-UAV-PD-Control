from sympy import symbols, Matrix, diff, cross, diag, skew, time_derivative, jacobian

# Wing Kinematics
t = symbols('t')

# Wing conformation length variables (constant)
L1, L2a, L2b, L2c, L3a, L3b, L3c, alpha_3 = symbols('L1 L2a L2b L2c L3a L3b L3c alpha_3')
L4, L5a, L5b, L5c, L5d, L6, L7a, L7b, L7c, offset_x = symbols('L4 L5a L5b L5c L5d L6 L7a L7b L7c offset_x')

# Joint positions (constant)
P1 = symbols('P1_1 P1_2')
P5 = symbols('P5_1 P5_2')
P8 = symbols('P8_1 P8_2')

# Joint angles (time varying)
theta_1, thetad_1, thetadd_1 = symbols('theta_1 thetad_1 thetadd_1')
theta_2, thetad_2, thetadd_2 = symbols('theta_2 thetad_2 thetadd_2')
theta_3, thetad_3, thetadd_3 = symbols('theta_3 thetad_3 thetadd_3')
theta_5, thetad_5, thetadd_5 = symbols('theta_5 thetad_5 thetadd_5')
theta_6, thetad_6, thetadd_6 = symbols('theta_6 thetad_6 thetadd_6')
theta_8, thetad_8, thetadd_8 = symbols('theta_8 thetad_8 thetadd_8')
theta_9, thetad_9, thetadd_9 = symbols('theta_9 thetad_9 thetadd_9')

theta = [theta_1, theta_2, theta_3, theta_5, theta_6, theta_8, theta_9]
theta_d = [thetad_1, thetad_2, thetad_3, thetad_5, thetad_6, thetad_8, thetad_9]
theta_dd = [thetadd_1, thetadd_2, thetadd_3, thetadd_5, thetadd_6, thetadd_8, thetadd_9]

# Joint position 4
P4_a = P1 + Matrix([[L1], [0]]) + rot_2d(theta_1) * Matrix([[L2a], [0]])
P4_b = P5 + rot_2d(theta_5 - alpha_3) * Matrix([[-L3b], [L3c]])

# Joint position 7
P7_a = P1 + rot_2d(theta_1) * Matrix([[L1], [0]]) + rot_2d(theta_2) * Matrix([[L2b], [L2c]]) + rot_2d(theta_3) * Matrix([[L4], [0]])
P7_b = P8 + rot_2d(theta_8) * Matrix([[-L5a], [L5b]])

# Joint position 10
P10_a = P8 + rot_2d(theta_8) * Matrix([[L5c], [L5d]]) + rot_2d(theta_9) * Matrix([[L6], [0]])
P10_b = P5 + rot_2d(theta_5) * Matrix([[L3a + L3c * sin(alpha_3)], [L3c * cos(alpha_3)]]) + rot_2d(theta_6) * Matrix([[L7b], [L7c]])

# Temporary function of time declarations
x_theta = symbols('x_theta', real=True)
ft_theta = lambda x, t: x_theta
xx0_theta = symbols('x_theta_1 x_theta_2 x_theta_3 x_theta_4 x_theta_5 x_theta_6 x_theta_7', real=True)
xx_theta = [ft_theta(xx, t) for xx in xx0_theta]
xxd_theta = [diff(xx, t) for xx in xx_theta]

# Constraint equation
pos_constraint = P4_a - P4_b + P7_a - P7_b + P10_a - P10_b + [theta_1]
vel_constraint = time_derivative(pos_constraint, theta + theta_d, theta_dd + xxd_theta, xx_theta, xxd_theta, t)
acc_constraint = time_derivative(vel_constraint, theta + theta_d, theta_dd + xxd_theta, xx_theta, xxd_theta, t)

# vel constraint = A_kv * theta_d + h_kv = u_kv
vel_theta1 = symbols('vel_theta1', real=True)
A_kv = Matrix(jacobian(vel_constraint, theta_d))
h_kv = vel_constraint.subs({theta_d[i]: 0 for i in range(len(theta_d))})
u_kv = Matrix([0, 0, 0, 0, 0, 0, vel_theta1])

# acc_constraint = A_k * theta_dd + h_k = u_k
u_theta1 = symbols('u_theta1', real=True)
A_k = Matrix(jacobian(acc_constraint, theta_dd))
h_k = acc_constraint.subs({theta_dd[i]: 0 for i in range(len(theta_dd))})
u_k = Matrix([0, 0, 0, 0, 0, 0, u_theta1])

# Aerobat Dynamics
mb, mh, mr, g = symbols('mb mh mr g', real=True)  # Mass scalar variables (constant)
Ib = Matrix(symbols('Ib_1 Ib_2 Ib_3', real=True))  # body inertia
Ih = Matrix(symbols('Ih_1 Ih_2 Ih_3', real=True))  # humerus wing segment inertia
Ir = Matrix(symbols('Ir_1 Ir_2 Ir_3', real=True))  # radius wing segment inertia

# Time varying states
pos_body = Matrix(symbols('x_1 x_2 x_3', real=True))  # body CoM linear position
vel_body = Matrix(symbols('xD_1 xD_2 xD_3', real=True))  # body CoM linear velocity
acc_body = Matrix(symbols('xDD_1 xDD_2 xDD_3', real=True))  # body CoM linear acceleration
w_body = Matrix(symbols('w_1 w_2 w_3', real=True))  # body angular velocity
wD_body = Matrix(symbols('wD_1 wD_2 wD_3', real=True))  # body angular acceleration
R_body = Matrix(symbols('Rb_11 Rb_12 Rb_13 Rb_21 Rb_22 Rb_23 Rb_31 Rb_32 Rb_33', real=True)).reshape(3, 3)  # body rotational matrix
RD_body = R_body * skew(w_body)  # R_body rate of change

# wing joint angles [humerus; radius], left only since we assume symmetry
theta_L = Matrix(symbols('theta_WL_1 theta_WL_2', real=True))  # wing angle (left)
thetaD_L = Matrix(symbols('thetaD_WL_1 thetaD_WL_2', real=True))  # wing ang vel (left)
thetaDD_L = Matrix(symbols('thetaDD_WL_1 thetaDD_WL_2', real=True))  # wing ang acc (left)

# Conformation parameters (body frame, relative to body center)
d1L = Matrix([offset_x, P5[0]]) + rot_x(theta_L[0]) * Matrix([0, L3a / 2 + L3c * sin(alpha_3), L3c * cos(alpha_3)])
d2L = d2L + rot_x(theta_L[1]) * Matrix([0, da_rad_L[0], da_rad_L[1], 0])
d1R = Matrix([offset_x, -P5[0], P5[1]]) + rot_x(-theta_L[0]) * Matrix([0, da_hum_R[0], da_hum_R[1], -L3c * sin(alpha_3), L3c * cos(alpha_3)])
d2R = d2R + rot_x(-theta_L[1]) * Matrix([0, da_rad_R[0], da_rad_R[1], 0])

# Temporary function of time declarations
xx = Matrix(symbols('xx_1 xx_2 xx_3 xx_4 xx_5 xx_6 xx_7 xx_8 xx_9 xx_10 xx_11 xx_12 xx_13 xx_14 xx_15 xx_16 xx_17 xx_18 xx_19 xx_20 xx_21 xx_22 xx_23 xx_24 xx_25 xx_26', real=True)).reshape(26, 1)
ft = lambda x, t: symbols('xx', real=True)
xx = [ft(xx[i], t) for i in range(26)]
xxd = [diff(xx[i], t) for i in range(26)]

# wing inertial mass positions
pos_hum_L = pos_body + R_body * d1L
pos_rad_L = pos_body + R_body * d3L
pos_hum_R = pos_body + R_body * d1R
pos_rad_R = pos_body + R_body * d3R

# wing inertial mass velocities
vel_hum_L = time_derivative(pos_hum_L, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
vel_rad_L = time_derivative(pos_rad_L, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
vel_hum_R = time_derivative(pos_hum_R, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
vel_rad_R = time_derivative(pos_rad_R, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)

# wing angular velocities (body frame)
angvel_hum_L = w_body + Matrix([thetaD_L[0], 0, 0])
angvel_rad_L = w_body + Matrix([thetaD_L[1], 0, 0])
angvel_hum_R = w_body + Matrix([-thetaD_L[0], 0, 0])
angvel_rad_R = w_body + Matrix([-thetaD_L[1], 0, 0])

# Energy equations
K_lin = mb / 2 * (vel_body[0] ** 2 + vel_body[1] ** 2 + vel_body[2] ** 2) + \
        mh / 2 * (vel_hum_L[0] ** 2 + vel_hum_L[1] ** 2 + vel_hum_L[2] ** 2) + \
        mr / 2 * (vel_rad_L[0] ** 2 + vel_rad_L[1] ** 2 + vel_rad_L[2] ** 2) + \
        mh / 2 * (vel_hum_R[0] ** 2 + vel_hum_R[1] ** 2 + vel_hum_R[2] ** 2) + \
        mr / 2 * (vel_rad_R[0] ** 2 + vel_rad_R[1] ** 2 + vel_rad_R[2] ** 2)

K_ang = w_body.dot(diag(Ib)).dot(w_body) / 2 + \
        angvel_hum_L.dot(diag(Ih)).dot(angvel_hum_L) / 2 + \
        angvel_rad_L.dot(diag(Ir)).dot(angvel_rad_L) / 2 + \
        angvel_hum_R.dot(diag(Ih)).dot(angvel_hum_R) / 2 + \
        angvel_rad_R.dot(diag(Ir)).dot(angvel_rad_R) / 2

U = mb * [0, 0, g].dot(pos_body) + \
    mh * [0, 0, g].dot(pos_hum_L) + \
    mr * [0, 0, g].dot(pos_rad_L) + \
    mh * [0, 0, g].dot(pos_hum_R) + \
    mr * [0, 0, g].dot(pos_rad_R)

L = K_lin + K_ang - U

# State positions (no body rotation)
q = theta_L + pos_body
qd = time_derivative(q, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
# Standard Euler-Lagrangian
dL_dx = jacobian(L, q)  # del L / del q
dL_dv = jacobian(L, qd)  # del L / del q_dot
dtdt_dLdv = time_derivative(dL_dv, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
eom1 = (dtdt_dLdv - dL_dx).T

# SO(3) dynamics
dL_dw_body = jacobian(L, w_body)  # del L / del w_body
dtdt_dLdw_body = time_derivative(dL_dw_body, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
term1 = dtdt_dLdw_body
term2 = cross(w_body, dL_dw_body)
term3 = Matrix([0, 0, 0])
for i in range(3):
    term3 = term3 + cross(R_body[i, :], jacobian(L, R_body[i, :])).T
eom2 = term1 + term2 + term3.T

# equation of motion
eom = eom1.row_insert(1, eom2)
qd_eom = thetaD_L + thetaDD_L + vel_body + acc_body + wD_body
qdd_eom = time_derivative(qd_eom, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)

# solve for: eom = M*accel + h = generalized forces
M = jacobian(eom, qdd_eom)
h = eom.subs({qdd_eom[i]: 0 for i in range(len(qdd_eom))})

# Aerodynamic Generalized Force
da_hum_L = Matrix(symbols('da_hL_1 da_hL_2', real=True))
da_rad_L = Matrix(symbols('da_rL_1 da_rL_2', real=True))
da_hum_R = Matrix(symbols('da_hR_1 da_hR_2', real=True))
da_rad_R = Matrix(symbols('da_rR_1 da_rR_2', real=True))

da1L = [offset_x, P5[0]] + rot_x(theta_L[0]) * [0, L3a / 2 + L3c * sin(alpha_3), L3c * cos(alpha_3)]
da2L = da2L + rot_x(theta_L[1]) * [0, da_rad_L[0], da_rad_L[1], 0]
da1R = [offset_x, -P5[0], P5[1]] + rot_x(-theta_L[0]) * [0, da_hum_R[0], da_hum_R[1], -L3c * sin(alpha_3), L3c * cos(alpha_3)]
da2R = da2R + rot_x(-theta_L[1]) * [0, da_rad_R[0], da_rad_R[1], 0]

pos_aero_hum_L = pos_body + R_body * da1L
pos_aero_rad_L = pos_body + R_body * da2L
pos_aero_hum_R = pos_body + R_body * da1R
pos_aero_rad_R = pos_body + R_body * da2R

vel_aero_hum_L = time_derivative(pos_aero_hum_L, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
vel_aero_rad_L = time_derivative(pos_aero_rad_L, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
vel_aero_hum_R = time_derivative(pos_aero_hum_R, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
vel_aero_rad_R = time_derivative(pos_aero_rad_R, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)

Ba_hL = jacobian(vel_aero_hum_L, qd_eom).T
Ba_rL = jacobian(vel_aero_rad_L, qd_eom).T
Ba_hR = jacobian(vel_aero_hum_R, qd_eom).T
Ba_rR = jacobian(vel_aero_rad_R, qd_eom).T

# Tail
da_tail = symbols('dtail_1 dtail_2 dtail_3', real=True)
pos_aero_tail = pos_body + R_body * da_tail
vel_aero_tail = time_derivative(pos_aero_tail, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
Ba_tail = jacobian(vel_aero_tail, qd_eom).T

# Thrusters
u_body = Matrix(symbols('ubody_1 ubody_2 ubody_3', real=True))
d_force = Matrix(symbols('dforce_1 dforce_2 dforce_3', real=True))
pos_body_force = pos_body + R_body * d_force
vel_body_force = time_derivative(pos_body_force, theta_L + thetaD_L + pos_body + vel_body + w_body, t, xx, xxd)
Ba_body_force = jacobian(vel_body_force, qd_eom).T
body_gen_force = Ba_body_force.dot(u_body)  # generalized force for forces acting on the body