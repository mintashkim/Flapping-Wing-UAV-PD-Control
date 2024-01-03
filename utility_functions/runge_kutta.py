from func_eom import *

# 4th order Runge-Kutta numerical integration
# aero_data contains the aerodynamic forces and applied positions of the first
# time integration subcomponent.

def march_rk4(xk, xd, xa, u1, f_thruster, t_thruster, yaw_damping, p):
    
    fk1, fd1, fa1, aero_data = func_eom(xk, xd, xa, u1, f_thruster, t_thruster, yaw_damping, p, nargout=4)
    fk2, fd2, fa2 = func_eom(xk + fk1*p.dt/2, xd + fd1*p.dt/2, xa + fa1*p.dt/2, u1, f_thruster, t_thruster, yaw_damping, p, nargout=3)
    fk3, fd3, fa3 = func_eom(xk + fk2*p.dt/2, xd + fd2*p.dt/2, xa + fa2*p.dt/2, u1, f_thruster, t_thruster, yaw_damping, p, nargout=3)
    fk4, fd4, fa4 = func_eom(xk + fk3*p.dt, xd + fd3*p.dt, xa + fa3*p.dt, u1, f_thruster, t_thruster, yaw_damping, p, nargout=3)

    xk_next = xk + (fk1/6 + fk2/3 + fk3/3 + fk4/6) * p.dt
    xd_next = xd + (fd1/6 + fd2/3 + fd3/3 + fd4/6) * p.dt
    xa_next = xa + (fa1/6 + fa2/3 + fa3/3 + fa4/6) * p.dt

    return xk_next, xd_next, xa_next, aero_data