import numpy as np
from utility_functions.rotation_matrix import *
from utility_functions.runge_kutta import march_rk4
from func_initial_joint_angles import *
from symbolic_functions.func_wing_kinematic_vel import *

class Data:
    def __init__(self, time, Ns, p):
        # Data recording struct (simulation states)
        self.t = time # time
        self.xk = np.zeros((14, Ns)) # KS states (joint angles and velocity)
        self.xd = np.zeros((22, Ns)) # Dynamical states (body and wing surfaces)
        self.xa = np.zeros((3*(p.n_blade - p.n_blade_tail), Ns)) # Wagner model states
        self.angles = np.zeros((3, Ns)) # Euler roll, pitch, yaw angles
        self.u1 = np.zeros((1, Ns)) # Crank / motor accel eration
        self.x = np.zeros((2, Ns))
        self.y = np.zeros((2, Ns))

        # Data recording struct (aerodynamic forces data)
        self.ua = np.zeros((8, Ns)) # Total generalized aerodynamic force
        self.ft = np.zeros((3*p.N_thruster, Ns)) # Thruster forces (inertial)
        self.tt = np.zeros((3*p.N_thruster, Ns))
        self.pt = np.zeros((3*p.N_thruster, Ns)) # Thruster forces position (inertial)
        self.accel = np.zeros((8, Ns)) # Dynamic state accelleration due to aerodynamic forces
        self.pa = np.zeros((3, p.n_blade, Ns)) # Blade elements aerodynamic force applied position (inertial)
        self.fa = np.zeros((3, p.n_blade, Ns)) # Blade elements forces (inertial)
        self.loadcell = np.zeros((6, Ns))
        self.aoa = np.zeros((p.n_blade, Ns)) # Angle of attack

        # Others
        self.i_end = Ns; # Simulation end step count, in case when simulation terminated early
        

def func_simulation(simulation_parameter, k_opt):
# function [cost, data] = func_simulation(p_base, k_opt)
# Flappy numerical simulation

# Inputs: 
# simulation = simulation base parameters. Will be edited using k_opt if used in an optimization
# k_opt = optimization varible valkes. usage depends on the optimization

# Output: 
# cost = cost value for optimization
# data = simulation data, to be used for plotting

    # Optimization setup (if used in an optimization)

    # Setup the optimization, if used
    cost = 0

    # Edit some of the parameters to be optimized
    p = simulation_parameter

    # Initialize simulation

    time = np.arange(0, p.t_end, p.dt) # Time stamp
    Ns = time.size # Simulation step count

    data = Data(time, Ns, p)

    # Equation of motion states and inputs
    xk = np.zeros(14) # Kinematic states (KS linkage angles)
    xd = np.zeros(22) # Dynamic states (body and wing states)
    xa = np.zeros(3 * (p.n_blade - p.n_blade_tail)) # Aerodynamics states (for Wagner model)
    uk = np.zeros(1) # Kinematic inputs (KS crank acceleration)

    # Initial rotation matrix
    Rb0 = rot_y(p.pitch_init) @ rot_x(p.roll_init) # Initial pitch and roll angle
    xd[13:22] = Rb0.flatten()

    # Initial states (kinematic), joint order: 1,2,3, ,5,6, ,8,9
    # Joint angles
    xk[0:7] = func_initial_joint_angles(p.theta1_init, p.wing_conformation).flatten()
    # Joint velocities
    kinematic_vel_out = func_wing_kinematic_vel(xk[0:7], p.flapping_freq * 2 * np.pi, p.wing_conformation.flatten())
    A_kv = kinematic_vel_out[0:7,:]
    h_kv = kinematic_vel_out[7,:]
    u_kv = kinematic_vel_out[8,:]
    xk[7:14] = np.linalg.solve(A_kv, -h_kv + u_kv)

    # Wing initial states (dynamics), follows the same values as KS states
    xd[0:2] = xk[3:5]      
    xd[5:7] = xk[10:12]

    # body initial velocity
    xd[7:10] = p.body_init_vel.flatten()

    # Numerical Simulation

    f_thruster = np.zeros((3,4)) # inertial force
    
    pitch_error = 0
    roll_error = 0
    pitch_error_i = 0
    roll_error_i = 0

    vel_error = 0
    vel_error_i = 0

    pitch_error_d = 0
    roll_error_d = 0

    wc = np.exp(-p.wc_lp * 2 * np.pi * p.dt)

    for i in range(Ns):
        # Display progress, 10% increments.
        # Comment out if you don't want to clutter the screen.
        if np.mod(i,(Ns - 1)/10) == 0:
            print("Progress: {percent}%".format(percent=i/(Ns - 1)*100))

        # Measurements
        
        # Calculate Euler angles from body rotation matrix
        R_body = xd[13:22].reshape((3,3)).T # body to inertial
        R_body_T = R_body.T
        pitch = np.arcsin(-R_body_T[0,2])
        yaw = np.arctan2(R_body_T[0,1], R_body_T[0,0])
        roll = np.arctan2(R_body_T[1,2], R_body_T[2,2])
        print("Pitch: {pitch}  |  Yaw: {yaw}  |  Roll: {roll}".format(pitch=pitch, yaw=yaw, roll=roll))
        
        # Prevent roll-over from -pi to pi
        if(roll > np.pi * 0.99):
            pitch = -np.pi - pitch

        # Controller

        # Forward Velocity
        vel_forward = xd[7]
        
        # Roll and pitch error
        pitch_error_old = pitch_error
        roll_error_old = roll_error
        pitch_error = p.pitch_ref - pitch
        roll_error = p.roll_ref - roll
        vel_error = p.vel_ref - vel_forward

        # Derivative, use low pass to prevent huge spikes
        pitch_rate = (pitch_error - pitch_error_old)/p.dt
        roll_rate = (roll_error - roll_error_old)/p.dt
        pitch_error_d = wc * pitch_error_d + (1-wc) * pitch_rate
        roll_error_d = wc * roll_error_d + (1-wc) * roll_rate
        
        # Integral
        pitch_error_i = pitch_error_i + pitch_error * p.dt
        roll_error_i = roll_error_i + roll_error * p.dt
        vel_error_i = vel_error_i + vel_error * p.dt
            
        # Desired u
        u_pitch = p.pitch_dir * (p.KP_P * pitch_error + p.KI_P * pitch_error_i + p.KD_P * pitch_error_d)
        u_roll = p.roll_dir * (p.KP_R * roll_error + p.KI_R * roll_error_i + p.KD_R * roll_error_d)
        u_vel = p.KP_V * vel_error + p.KI_V * vel_error_i
        
        
        # Pitch output forces
        x1 = u_vel + u_pitch 
        x2 = u_vel - u_pitch
        x1 = np.clip(x1, 0.0001, p.thruster_max_force)
        x2 = np.clip(x2, 0.0001, p.thruster_max_force)
        
        # Roll output forces
        y1 = u_roll
        y2 = -u_roll
        y1 = np.clip(y1, 0.0001, p.thruster_max_force)
        y2 = np.clip(y2, 0.0001, p.thruster_max_force)
        

        # Flapping speed tracking controller (adjust crank acceleration)
        if(time[i] > p.wait_time):
            u1 = p.kd * (p.flapping_freq * 2 * np.pi - xk[7])  # Crank acceleration
        else:
            u1 = 0
        
        # Pre-multiply by R_body for inertial frame forces
        f_thruster = np.zeros((3,4)) # inertial force (inertial)
        f1 = np.matmul(R_body, np.array([x1, 0, 0]).T).reshape(3,1) # Top propeller, thrust direction x positive.(x1>0)
        f2 = np.matmul(R_body, np.array([x2, 0, 0]).T).reshape(3,1) # Bottom propeller, thrust direction x positive.(x2>0)
        f3 = np.matmul(R_body, np.array([0, -y1, 0]).T).reshape(3,1) # Left propeller, thrust direction y negative.(y1<0)
        f4 = np.matmul(R_body, np.array([0, y2, 0]).T).reshape(3,1) # Right propeller, thrust direction y positive.(y2>0)
        f_thruster = np.concatenate([f1, f2, f3, f4], axis=1) # 3x4
        t_thruster = np.zeros((3,4)) # thruster torque (body)
        # Rotational damping in yaw.
        yaw_rate = xd[12]
        kd = 0.15
        yaw_damper = -kd * yaw_rate
        yaw_damping = np.array([0, 0, yaw_damper]).reshape(3,1)

        # Record state and control data
        
        data.xk[:,i] = xk
        data.xd[:,i] = xd
        data.xa[:,i] = xa 
        data.u1[:,i] = u1
        
        data.x[:,i] = np.array([x1, x2])
        data.y[:,i] = np.array([y1, y2])

        data.ft[:,i] = f_thruster[:].flatten()
        data.tt[:,i] = t_thruster[:].flatten()

        data.angles[:,i] = np.array([roll, pitch, yaw])
    

        # Simulation time march (RK4)
        
        # States of the next time step
        # print(xk) # 4, 11, 13 initial value error
        # print(xd) # 4, 11, 13 initial value error

        xk, xd, xa, aero_data = march_rk4(xk, xd, xa, u1, f_thruster, t_thruster, yaw_damping, p)

        # Check if states become unstable and terminate simulation
        if(np.sum(np.isnan(xk)) or np.sum(np.isnan(xd)) or np.sum(np.isnan(xa))):
            print("UNSTABLE, TERMINATING EARLY")
            data.i_end = i
            break # If unstable, end early

        pitch_rate_current = xd[11]

        # Record aerodynamics data (from the RK4 calculations)
        data.ua[:,i] = aero_data.ua
        data.accel[:,i] = aero_data.accel
        data.pa[:,:,i] = aero_data.pa
        data.fa[:,:,i] = aero_data.fa
        data.aoa[:,i] = aero_data.aoa

    return cost, data