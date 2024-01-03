from func_create_wing_segment import *
from func_simulation import *
from plot_data import *
from plot_animation import *
import numpy as np

class Simulation_Parameter:
    def __init__(self):
        # Setup

        # Flags during parameter setup
        self.flag_use_optimized_model = 0 # 0 = dickinson's model, 1 = load cell optimization model
        self.flag_reverse_wing = 0 # 0 = aerobat's forward sweep, 1 = reverse sweep

        # Flags during simulation
        self.flag_use_wagnermodel = 1 # 1 = wagner model, 0 = quasi-steady model

        # Tail stabilizer
        # NOTE: Tail uses quasi-steady aerodynamics
        self.flag_use_tail_stabilizer = 1 # 1 if to add a tail stabilizer

        # Constraint dynamic states (prevent change in states, assuming zero IC)
        # 0 = no constraint
        # 1 = constraint lateral movements, roll, and yaw
        # 2 = fully constraint position and attitude
        # 3 = constraint orientation only
        # 4 = constraint position only
        # 5 = constraint all but z position
        self.flag_constraint_dynamics = 0

        # PI controller velocity tracking
        self.vel_ref = 3
        self.KP_V = 2 
        self.KI_V = 2

        # PID controller pitch tracking
        self.pitch_ref = np.deg2rad(-20)
        self.pitch_dir = 1
        self.KP_P = 8
        self.KD_P = 0.1
        self.KI_P = 1

        # PID controller roll tracking
        self.roll_ref = np.deg2rad(0)
        self.roll_dir = 1
        self.KP_R = 2
        self.KD_R = 0.1
        self.KI_R = 0.1

        self.thruster_max_force = 0.5
        self.wc_lp = 20 # low pass cutoff frequency hz


        # Simulation parameters (general)

        self.pitch_init = np.deg2rad(-16) # initial pitch angle. positive is pitch down
        self.roll_init = np.deg2rad(-5) # initial roll angle, positive is roll right(back side view)
        self.theta1_init = np.deg2rad(-90) # initial crank gear angle
        self.body_init_vel = np.array([2.08, 0, 0]).reshape(3, 1) # initial body velocity (m/s)

        # Simulation time parameters
        self.dt = 0.00002 # simulation time step (s). smaller is more stable but slower
        self.t_end = 10 # simulation end time (s)
        self.wait_time = 0.01 # wait time before start applying flapping speed tracking controller

        # Airflow settings
        self.airspeed = np.array([0, 0, 0]).reshape(3, 1) # freestream air speed (inertial frame)
        self.air_density = 1 # density of air (kg/m^3)

        # Flapping and crank controller gains
        self.kd = 400 # crank speed tracking controller gain
        self.flapping_freq = 4.75 # flapping frequency in Hz

        self.flag_use_tail_stabilizer = True

        if(self.flag_use_tail_stabilizer):
            self.n_blade_tail = 2 # tail stabilizer blade elements
        else:
            self.n_blade_tail = 0

        # Blade element numbers (for lifting line theory, assume constant blade elements)
        self.n_blade_prox = 2 # number of blade elements per side (proximal)
        self.n_blade_dist = 6 # number of blade elements per side (distal)

        # Total blade elements
        self.n_blade = 2 * (self.n_blade_prox + self.n_blade_dist) + self.n_blade_tail

        # Length vectors to thrusters in body frame
        self.N_thruster = 4
        self.d_thruster = np.zeros((3, self.N_thruster)) # [d1,d2,d3] = body frame displacement from body com to force pos, DI = [DX;DY;DZ]
    
        # d1 = np.array([-0.15, 0, 0.015])
        # d2 = np.array([-0.13, 0, -0.03])
        d1 = np.array([-0.15, 0, 0.02]).reshape(3,1)
        d2 = np.array([-0.13,  0, -0.03]).reshape(3,1)
        d3 = np.array([0, 0.01, 0.07]).reshape(3,1)
        d4 = np.array([0, -0.01, 0.07]).reshape(3,1)
        self.d_thruster = np.concatenate([d1, d2, d3, d4], axis=1)


        # KS Dimension parameters

        # NOTE: taken from solidworks

        # Linkage length parameters (meter)
        self.L1  = 7.5/1000
        self.L2a = 35/1000
        self.L2b = 5/1000 * np.cos(np.deg2rad(35))
        self.L2c = 5/1000 * np.sin(np.deg2rad(35))
        self.L3a = 50/1000
        self.L3b = 6.71/1000
        self.L3c = 6/1000
        self.L4  = 36/1000
        self.L5a = 11.18/1000
        self.L5b = 10/1000
        self.L5c = 78.78/1000
        self.L5d = -3.89/1000
        self.L6  = 36/1000
        self.L7a = 90/1000
        self.L7b = 35/1000
        self.L7c = 15/1000

        # Stationary joint positions (x,y) (meter)
        self.P1  = (np.array([-10, -9.8]) / 1000).reshape(2,1)
        self.P5  = (22 * np.array([np.cos(np.deg2rad(20)), np.sin(np.deg2rad(20))]) / 1000).reshape(2,1)
        self.P8  = (22 * np.array([np.sin(np.deg2rad(25)), np.cos(np.deg2rad(25))]) / 1000).reshape(2,1)

        # Wing origin offset about x-axis relative to the body COM
        self.offset_x  = np.array([-42/1000]).reshape(1,1) # wing com offset in x-axis (meter) -15 old value

        # Other parameters
        self.alpha_3 = np.deg2rad(-15) # humerus link bend angle

        # For symbolic generated function
        L = np.array([self.L1, self.L2a, self.L2b, self.L2c, self.L3a, self.L3b, self.L3c, self.alpha_3, self.L4, 
                    self.L5a, self.L5b, self.L5c, self.L5d, self.L6, self.L7a, self.L7b, self.L7c]).reshape(17,1)
        self.wing_conformation = np.concatenate([L, self.P1, self.P5, self.P8, self.offset_x], axis=0)

        # Body and aerodynamic blade elements dimensions 

        # Body dimensions (meter)
        self.body_length = 40/1000  
        self.body_radius = 25/1000

        # Wing segment span (meter)
        self.span_prox = 50/1000
        self.span_dist = 132/1000

        # Proximal/inner wing segment front and rear x and y positions
        self.wing_xf_prox = np.array([63.5, 76.5]) / 1000 # front
        self.wing_xr_prox = np.array([-95.6, -111]) / 1000 # rear
        self.wing_y_prox  = np.array([0, 50]) / 1000 # y axis positions

        # Distal/outer wing segment front and rear x and y positions
        self.wing_xf_dist = np.array([76.5, 80.7, 118, 101]) / 1000 # front
        self.wing_xr_dist = np.array([-111, -115, 17.05, 101]) / 1000 # rear
        self.wing_y_dist  = np.array([0, 8.5, 84, 132]) / 1000 # y axis positions

        # Tail stabilizer dimensions (meter)
        self.tail_center = np.array([-0.12, 0, 0.01]) # tail center position
        self.tail_width = 0.050 # tail width
        self.tail_chord = 0.030 # tail chord

        # Reverse the wing sweep
        if(self.flag_reverse_wing):
        # Switch positions between wing's front and rear edge positions
            temp1 = self.wing_xf_prox
            temp2 = self.wing_xr_prox
            self.wing_xf_prox = -temp2
            self.wing_xr_prox = -temp1
            
            temp1 = self.wing_xf_dist
            temp2 = self.wing_xr_dist
            self.wing_xf_dist = -temp2
            self.wing_xr_dist = -temp1
        
        # Maximum wing span
        self.span_max = 2 * (self.body_radius + self.span_dist + self.span_prox)

        # create the wing segment data for calculating aerodynamic forces
        self.strip_id, self.strip_dc, self.strip_c, self.strip_theta = func_create_wing_segment(self);


        # Aerodynamic and mass parameters  
        
        # Wagner function coefficients, a*exp(b*t_bar)
        self.phi_a = np.array([-0.165, -0.335])   # for wagner lift model
        self.phi_b = np.array([-0.045, -0.3])
        self.Phi_0 = 1 + sum(self.phi_a)
        self.a0 = 2 * np.pi  # slope of the lift coefficient. approximately 2*pi for thin airfoil
        self.c0 = 160/1000  # base chord along wing's axis of symmetry (i'm using proximal wing's base chord)

        # Dynamic parameters 

        # Total mass around 24 grams
        self.mb = 60/1000 # 20g
        self.mh = 2/1000 # 1.25/1000; # times 2
        self.mr = 4/1000 # (2.25 + 1.5)/1000; # times 2
        self.m_total = self.mb + 2 * (self.mh + self.mr)

        # Body inertia tensor, assume solid cylinder
        self.Ib = self.mb * np.array([self.body_radius**2 / 2,
                                      (3 * self.body_radius**2 + self.body_length**2) / 12,
                                      (3 * self.body_radius**2 + self.body_length**2) / 12]).reshape(3,1)

        # Humerus wing assume thin plate, 70mm chord, 50mm span, 0.5 mm thick
        self.Ih = self.mh/12 * np.array([(0.05**2 + 0.0005**2), (0.07**2 + 0.0005**2), (0.05**2 + 0.07**2)]).reshape(3,1)

        # Radius wing assume thin plate, 70mm chord, 130mm span, 0.5 mm thick
        self.Ir = self.mr/12 * np.array([(0.13**2 + 0.0005**2), (0.07**2 + 0.0005**2), (0.13**2 + 0.07**2)]).reshape(3,1)

        # Other
        self.g = 9.8

        self.params = np.concatenate([np.array([self.mb, self.mh, self.mr, self.g]).reshape(4,1), self.Ib, self.Ih, self.Ir], axis=0)
        self.chord_length = 70*10^-3

        # Aerodynamic modeling (lift and drag coefficient model)
        if(self.flag_use_optimized_model):
            # From "aero_opt_04-15-2021_exp2_v2.mat" in the simulation_data folder
            self.aero_model_lift = np.array([-0.4802, 6.5892, 2.4974, -6.3210])
            self.aero_model_drag = np.array([2.2481, -1.7594, 2.4843, -9.3407])
        else:
            # Dickinson model
            self.aero_model_lift = np.array([0.225, 1.58, 2.13, -7.2])
            self.aero_model_drag = np.array([1.92, -1.55, 2.04, -9.82])
    
if __name__ == '__main__':
    # Simulation
    print('Start Simulation')
    simulation_parameter = Simulation_Parameter()

    # This is the function to run the simulation. 
    k_opt = 0
    cost, data = func_simulation(simulation_parameter, k_opt)

    # NOTE: Can be modified to be used in an optimization function by using the 
    # "k_opt" and "cost" variables. Works as long as the cost value is the first
    # output of the function.
    # Example: k_opt = fmincon(  @(k) func_simulation(p, k), ... );   
    # This main script can be included in the function for optimization as well.

    print("Simulation complete!")

    # Post process data

    # Sum the total force acting on the body from all blade elements
    data.force_sum = np.zeros((3,len(data.t)))
    for i in range(len(data.t)):
        data.force_sum[:,i] = sum(data.fa[:,:,i],2)

    # Aerodynamic lift and drag (convert from newton to grams)
    mean_lift = np.mean(data.force_sum[2,:])
    mean_drag = np.mean(data.force_sum[0,:])
    print('Mean lift = {mean_lift} grams'.format(mean_lift = mean_lift/9.8*1000))
    print('Mean drag = {mean_drag} grams'.format(mean_drag = mean_drag/9.8*1000))

    # Plotting

    # save('simulation_data/dataset_cdc_test.mat', 'data', 'p');

    # plot_data()
    # plot_animation()

# TODO
# 1. Array index -> find :
# 2. Data class
# 3. save
# 4. np.zeros() dimension