import numpy as np
from utility_functions.rotation_matrix import *
from symbolic_functions.func_wing_kinematic import *
from symbolic_functions.func_wing_rhRL import *
from symbolic_functions.func_tail import *
from symbolic_functions.func_Mh import *
from symbolic_functions.func_genforce_body import *

class AeroData:
    def __init__(self, p):
        self.ua = np.zeros((8, 1)) # Total generalized forces from aerodynamics
        self.accel = np.zeros((8, 1)) # State accel without body constraint

        # Sectional aerodynamic forces and position (inertial)
        self.pa = np.zeros((3, p.n_blade)) # Position
        self.fa = np.zeros((3, p.n_blade)) # Forces
        self.aoa = np.zeros((1, p.n_blade)) # Angle of attack

    def __str__(self):
        return "ua: {ua}  |  accel: {accel}  |  pa: {pa}  |  fa: {fa}  |  aoa: {}".format(ua=self.ua, accel=self.accel, pa=self.pa, fa=self.fa)
        

def func_eom(xk, xd, xa, u1, u_thruster, u_torque, yaw_damping, p, nargout):
    
    aero_data = AeroData(p)

    w_body = xd[10:13]
    R_body = (xd[13:22].reshape(3, 3))
    
    Ak, hk, uk = func_wing_kinematic(xk, u1, p.wing_conformation.flatten())
    ak = (np.linalg.lstsq(Ak.astype('float64'), (uk.astype('float') - hk.astype('float')), rcond=None))[0]
    
    # Determine fk = [vel ; accel]
    fk = xk * 0
    fk[0:7] = xk[7:14] # KS states velocity
    fk[7:14] = ak.flatten()   # KS states acceleration

    # Wagner's model EOM and aerodynamics
    nWagner = p.n_blade - p.n_blade_tail  # Wagner model length

    # Intitialize output
    fa = xa * 0 # Output dxa / dt
    ua = 0  # Sum of generalized aerodynamics forces
    xa_m = xa.reshape(nWagner, 3) # Reshape xa for easy access
    fa_m = xa_m * 0
    
    # Calculate blade element inertial position, velocities, and Ba(mapping for generalized forces)
    # Ba is the inertial force to generalized force transformation matrix for calculating ua
    pos_s = np.zeros((3, p.n_blade)) # Inertial position
    vel_s = np.zeros((3, p.n_blade)) # Inertial velocity
    Ba_s = np.zeros((8, 3, p.n_blade)) #  ua = Ba * fa, fa = inertial aero force
    vel_s_surf = np.zeros((3, p.n_blade)) # Velocity about the wing axis[front, left, normal]
    e_n = np.zeros((3, p.n_blade)) # Wing's surface normal direction
    aoa = np.zeros((1, p.n_blade)) # Angle of attack(free stream)
    aoa_d = np.zeros((p.n_blade)) # Change in angle of attack due to downwash
    U = np.zeros((p.n_blade)) # Effective air speed
    e_effvel = np.zeros((3, p.n_blade)) # (3,18)
   
    # pos = np.zeros((3,1))
    # vel = np.zeros((3,1))
    # Ba = np.zeros((8,3))
    # Rw = np.zeros((3,3))

    for i in range(p.n_blade):
        # Check which wing segment this blade element index belongs to
        if p.strip_id[i] == 1:
            pos, vel, Ba = func_wing_hL(xd, p.strip_dc[0:2,i], p.wing_conformation.flatten())
            Rw = rot_x(xd[0])
        elif p.strip_id[i] == 2:
            pos, vel, Ba = func_wing_rL(xd, p.strip_dc[0:2,i], p.wing_conformation.flatten())
            Rw = rot_x(xd[1])
        elif p.strip_id[i] == 3:
            pos, vel, Ba = func_wing_hR(xd, p.strip_dc[0:2,i], p.wing_conformation.flatten())
            Rw = rot_x(-xd[0])
        elif p.strip_id[i] == 4:
            pos, vel, Ba = func_wing_rR(xd, p.strip_dc[0:2,i], p.wing_conformation.flatten())
            Rw = rot_x(-xd[1])
        elif p.strip_id[i] == 5:
            pos, vel, Ba = func_tail(xd, p.strip_dc[0:3,i])
            Rw = np.eye(3)
 
        # Wing velocities about wing segment's axis [front, left, normal]
        # Local effective velocity(wing frame)
        vel_surf = (Rw.T) @ (R_body.T) @ (vel - p.airspeed)
        vel_surf_norm = np.linalg.norm(vel_surf, 2) # Air velocity magnitude
        # print(vel_surf_norm)

        # Effective velocity direction(x, z)
        if np.linalg.norm(vel_surf.flatten().take([0,2])) < 1e-6:
            ev = np.zeros(3,1) # Prevent dividing by zero
        else:
            ev = vel_surf / np.linalg.norm(vel_surf, 2) # Unit vector of air flow direction
        alpha_inf = np.arctan2(-vel_surf[2], vel_surf[0]) # Angle of attack
        
        # Record values
        pos_s[:, i] = pos.flatten()
        vel_s[:, i] = vel.flatten() - p.airspeed.flatten()
        vel_s_surf[:, i] = vel_surf.flatten()
        Ba_s[:,:, i] = Ba
        e_n[:, i] = R_body @ Rw @ np.array([0,0,1]).reshape(3,1).flatten()
        aoa[:, i] = alpha_inf
        e_effvel[:, i] = ev.flatten()
        if vel_surf_norm < 1e-4:
            U[i] = 1e-4
        else:
            U[i] = vel_surf_norm

    # Determine An and bn for An * an_dot = bn to solve for fa
    # Follows Boutet formulations

    if(p.flag_use_wagnermodel):

        An = np.zeros((nWagner, nWagner))
        a0 = p.a0
        c0 = p.c0
        n = np.arange(1, nWagner+1)
        for i in range(nWagner):
            An[i,:] = a0 * c0 / U[i] * np.sin(n * p.strip_theta[i])
        
        an = xa_m[:,0]
        bn = np.zeros((nWagner, 1))
       
        for i in range(nWagner):

            # Aero states and effective air speed
            z1 = xa_m[i, 1]
            z2 = xa_m[i, 2]
            # Downwash due to vortex
            wy = -a0 * c0 * U[i] / 4 / p.span_max * (n * np.sin(n * p.strip_theta[i])) / np.sin(p.strip_theta[i]) @ an
            wn = vel_s_surf[2, i]
            w = wn + wy
            # print(w)
            
            # vel_surf = (Rw.T).dot(R_body.T).dot(vel.reshape(3,1) - p.airspeed.reshape(3,1))
            aoa = aoa.reshape(18,1)
            alpha_downwash = np.arctan2(-w, vel_surf[0]) - aoa[i]
            aoa_d[i] = alpha_downwash

            bn[i]= -a0 * c0 / p.strip_c[i] * (np.sin(n * p.strip_theta[i])) @ an + a0 *(w * p.Phi_0 / U[i] + p.phi_a[0] * p.phi_b[0] / (p.strip_c[i] / 2) * z1 + p.phi_a[1] * p.phi_b[1] * z2)
            fa_m[i, 1] = U[i] * wn + p.phi_b[0] / (p.strip_c[i] / 2) * z1
            fa_m[i, 2] = U[i] * wy + p.phi_b[1] / (p.strip_c[i] / 2) * z2
        
        # Calculate aerodynamic states rate of change for an
        fa_m[:,0] = np.linalg.lstsq(An.astype('float32'), bn.astype('float32'), rcond=None)[0].flatten()
        fa = fa_m[:].T.flatten()
    else:
        fa = xa * 0

    # Calculate aerodynamics lift and drag
    pa = np.zeros((3, p.n_blade))
    f_aero = np.zeros((3, p.n_blade))

    for i in range(p.n_blade):
        if p.strip_id[i] == 1 or p.strip_id[i] == 3:
            delta_span = p.span_prox / p.n_blade_prox
        elif p.strip_id[i] == 2 or p.strip_id[i] == 4:
            delta_span = p.span_dist / p.n_blade_dist
        else:
            delta_span = p.tail_width / p.n_blade_tail

            # Rotation matrix from body axis to wing axis
        if p.strip_id[i] == 1: # hL
            Rw = rot_x(xd[0])
        elif p.strip_id[i] == 2: # rL
            Rw = rot_x(xd[1])
        elif p.strip_id[i] == 3: # hR
            Rw = rot_x(-xd[0])
        elif p.strip_id[i] == 4: # rR
            Rw = rot_x(-xd[1])
        elif p.strip_id[i] == 5: # tail
            Rw = np.eye(3)

        # Lift coefficient(wagner for the wing, if enabled)
        if (i < nWagner) and p.flag_use_wagnermodel:
            Gamma = 0.5 * a0 * c0 * U[i] * np.sin(n * p.strip_theta[i]) @ an
            C_lift = -a0 * np.sin(n * p.strip_theta[i]) @ (an + p.strip_c[i] / U[i] * fa_m[:, 0])
        else:
            # Quasi - steady model
            Gamma = 0
            C_lift = lift_coeff(aoa[i], p.aero_model_lift)  # p.aero_model_lift(1,4), G_lift,a number

        # Drag coefficient(quasi - steady)
        aoa = aoa.flatten()
        C_drag = drag_coeff(aoa[i], p.aero_model_drag)
        e_lift = np.cross(e_effvel[:, i].flatten(), np.array([0,1,0])) 
        e_drag = (-e_effvel[:, i]).flatten()  # (3,1)
        lift = p.air_density / 2 * U[i]**2 * C_lift * delta_span * p.strip_c[i] * e_lift
        drag = p.air_density / 2 * U[i]**2 * C_drag * delta_span * p.strip_c[i] * e_drag

        # Combine the drag and lift directions in inertial frame
        f_aero[:,i] = R_body @ Rw @ (drag + lift) # (3,1)
        ua = ua + Ba_s[:, :, i]@(f_aero[:,i].reshape(3,1)) # （8,1)
        pa[:,i] =  pos_s[:,i]


    Md, hd = func_Mh(xd, p.params.flatten(), p.wing_conformation.flatten())
    # Apply joint constraints
    uj = fk[10:12] # Joint acceleration constraint from kinematics

    Jc = np.concatenate((np.eye(2), np.zeros((2,6))),axis=1)  # Constraint jacobian

    ut = hd * 0
    for i in range(p.N_thruster):
        ut = ut + func_genforce_body(xd, p.d_thruster[:,i], u_thruster[:,i])

    for i in range(p.N_thruster):
        ut = ut + ((np.concatenate([np.zeros((3,5)), np.eye(3)], axis=1)).T @ (u_torque[:,i].reshape(3,1))).reshape(8,1)

    ut = ut + (np.concatenate((np.zeros((3,5)), np.eye(3)), axis=1)).T@(yaw_damping.reshape(3,1))
    ut = ut.astype('float64')
    u_a_t_d = np.array([ua + ut - hd]).reshape(8,1)
    
    u_j = uj.reshape(2,1)

    hc = np.concatenate((u_a_t_d, u_j), axis=0)
    hc = hc.astype('float64')
    
    # Optional: add other states constraints
    if p.flag_constraint_dynamics == 1:
        # Constraint lateral movements, roll, and yaw
        Jc = np.concatenate((Jc, np.concatenate((np.zeros((2,5)), np.array([[1,0,0],[0,0,1]])),axis=1)),axis=0)
        im_eye5 = np.concatenate((np.zeros((2, 2)), np.eye(2)), axis=1)
        im_eye6 = np.concatenate((im_eye5, np.zeros((2, 4))), axis=1)
        Jc = np.concatenate((Jc,im_eye6),axis=0)       #(6,8)
        hc = (np.concatenate((hc,np.zeros((4,1))),axis=0)).reshape(10,1)

    elif p.flag_constraint_dynamics == 2:
        Jc=np.concatenate((Jc, np.concatenate((np.zeros((6,2)), np.eye(6)),axis=1)),axis=0)
        hc=(np.concatenate((hc,np.zeros((6,1))),axis=0)).reshape(10,1)

    elif p.flag_constraint_dynamics == 3:
        Jc=np.concatenate((Jc, np.concatenate((np.zeros((3,5)), np.eye(3)),axis=1)),axis=0)
        hc=(np.concatenate((hc,np.zeros((3,1))),axis=0)).reshape(10,1)

    elif p.flag_constraint_dynamics == 4:
        im_eye1=np.concatenate((np.zeros((3,2)), np.eye(3)),axis=1)
        im_eye2=np.concatenate(( im_eye1, np.zeros((3,3))),axis=1)
        Jc = np.concatenate((Jc,im_eye2),axis=0)
        hc=(np.concatenate((hc,np.zeros((3,1))),axis=0)).reshape(10,1)

    elif p.flag_constraint_dynamics == 5:
        im_eye3 = np.concatenate((np.zeros((2,2)), np.eye(2)), axis=1)
        im_eye4 = np.concatenate((im_eye3, np.zeros((2, 4))), axis=1)
        Jc = np.concatenate((Jc, im_eye4), axis=0)
        Jc=np.concatenate((Jc, np.concatenate((np.zeros((3,5)), np.eye(3)),axis=1)),axis=0)
        hc = np.concatenate((hc, np.zeros((5, 1))), axis=0)

    n_jc = len(Jc[:, 0])
    Md_jc = np.concatenate((Md,-Jc.T),axis=1)
    Jc_zeros = np.concatenate((Jc,np.zeros((n_jc, n_jc))),axis=1)
    Mc = np.concatenate((Md_jc, Jc_zeros),axis=0)
    temp = (np.linalg.lstsq(Mc, hc, rcond=None))[0]
    accel = temp[0:8]
    lambda_ = temp[8:]
    R_body_D = R_body.dot(skew(w_body))
    hc=hc.astype('float64')

    # Determine fd
    fd = xd * 0
    fd[0:5] = xd[5:10]
    fd[5:13] = accel.flatten()

    # print(R_body_D[:])
    fd[13:22] = R_body_D.flatten()

    # For plotting
    if(nargout > 3):
        aero_data.ua = ua.flatten() # Total generalized forces
        aero_data.accel = np.linalg.inv(Md)@(ua - hd + Jc[0:2,:].T@lambda_[0:2]).flatten() # Accel without body constraint
        aero_data.pa = pa
        aero_data.fa = f_aero
        aero_data.aoa = aoa
        return fk, fd, fa, aero_data
    else:
        return fk, fd, fa


def strip_states(i, pw, xd, p):

    # Returns strip plunge position and pitch angle, with their velocities.
    q = np.array([0,0]).reshape(2,1)
    dq =np.array([0,0]).reshape(2,1)
    # Calculate inertial position and velocity

    if i <= p.n_blade_prox:
        pos, vel, Ba = func_wing_hL(xd, pw.reshape(2,), p.wing_conformation.reshape(24,1))
    elif i <= p.n_blade_prox + p.n_blade_dist:
        pos, vel, Ba = func_wing_rL(xd, pw.reshape(2,), p.wing_conformation.reshape(24,1))
    elif i <= 2 * p.n_blade_prox + p.n_blade_dist:
        pos, vel, Ba = func_wing_hR(xd, pw.reshape(2,), p.wing_conformation.reshape(24,1))
    elif i <= 2 * p.n_blade_prox + 2 * p.n_blade_dist:
        pos, vel, Ba = func_wing_rR(xd, pw.reshape(2,), p.wing_conformation.reshape(24,1))

    # Calculate body frame wing position and velocity
    Rb = (xd[13:21].reshape(3, 3)).T
    pos_body = Rb.T*(pos - xd[2:4])
    vel_body = Rb.T*(vel - xd[7:9])

    # pos_body = pos
    # vel_body = vel

    # Calculate pitch angle(angle of attack?)
    Rb_T = Rb.T
    pitch = np.asin(-Rb_T[0, 2])

    # Calculate plunge, pitch, and their rates
    q = np.array([pos_body[2],pitch]).reshape(2,1)
    dq = np.array([vel_body[2],xd[11]]).reshape(2,1)

    return q, dq, Ba

def func_wing_surface_pos(id, s, pos_wing, x_pitch):
    # Inputs:
    # id = wing id, id = 1: LH, 2: LR, 3: RH, 4: RR
    # s = unitless span position
    # pos_wing = wing vertices, in wing frame
    # x_pitch = pitch axis position in x - axis, in wing frame
    # outputs
    # pw = [x; y] or [front;left] position on the wing, in wing local frame
    # c = strip chord length
    # xe = distance from mid-chord point to the pitch axis along the wing

    pos_wing_hL = pos_wing[:,0:5]
    pos_wing_rL = pos_wing[:,5:11]
    pos_wing_hR = pos_wing[:,11:16]
    pos_wing_rR = pos_wing[:,16:22]

    if id == 1:
        # Left humerus wing
        id_f = np.array([2, 1]).reshape(1,2)# Front side index
        id_r = np.array([3, 4, 0]).reshape(1,3) # Rear side index
        pos_temp = pos_wing_hL
    elif id == 2:
        # Left radius wing
        id_f = np.array([0, 5, 4]).reshape(1,3)# Front side index
        id_r = np.array([1, 2, 3, 4]).reshape(1,4) # Rear side index
        pos_temp = pos_wing_rL
    elif id == 3:
        # Right humerus wing
        id_f = np.array([1, 2]).reshape(1,2) # Front side index
        id_r = np.array([0, 4, 3]).reshape(1,3) # Rear side index
        pos_temp = pos_wing_hR
    elif id == 4:
        # Right radius wing
        id_f = np.array([4, 5, 0]).reshape(1,3) # Front side index
        id_r = np.array([4, 3, 2, 1]).reshape(1,4) # Rear side index
        pos_temp = pos_wing_rR

    # Determine chord length at the specific wingspan
    span = pos_temp[1, id_f[-1]-1] - pos_temp[1, id_f[0]-1]
    pos_span = span * s
    xf, xr = interp_wingpos(pos_temp[:, id_f], pos_temp[:, id_r], pos_span)

    # Outputs
    pw = np.array([x_pitch,pos_span]).reshape(2,1) # Position of the pitch axis and span on the wing local frame
    c = xf - xr # Strip chord length
    xe = (xf + xr) / 2 - x_pitch # Pitch axis relative to the chord center

    return pw, c, xe

def interp_wingpos(xf, xr, y):
    # id order: left humerus-radius, right humerus-radius
    # outputs front and back side positions in front side
    i = 2
    j = 2
    while (y > xf[1,i-1] and i < len(xf[1,:])):
        i = i+1
    while (y > xr[1,j-1] and j < len(xr[1,:])):
        j = j+1
    qf = (xf[0,i-1] - xf[0, i-1-1]) * (y - xf[1,i-1-1]) / (xf[1,i-1] - xf[1,i-1-1]) + xf[0,i-1-1]
    qr = (xr[0,j-1] - xr[0, j-1-1]) * (y - xr[1,j-1-1]) / (xr[1,j-1] - xr[1,j-1-1]) + xr[0,j-1-1]

    return qf, qr

def lift_coeff(a, params):
    c = params[0]+ params[1] * np.sin(np.radians(params[2] * np.degrees(a) + params[3]))
    return c

def drag_coeff(a, params):
    c = params[0] + params[1] * np.cos(np.radians(params[2] * np.degrees(a) + params[3]))
    return c