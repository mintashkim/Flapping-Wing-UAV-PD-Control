import numpy as np
from scipy.interpolate import interp1d

def func_create_wing_segment(p):
    # Inputs:
    # p = parameters containing strip data and wing conformation
    # Outputs:
    # id = wing segment id, {hL, rL, hR, rR} = {1, 2, 3, 4}
    # dc = distance vector from wing segment origin to quarter chord force position.
    # c = segment chord length
    # theta = acos(2 * y / fullspan),for lifting line theory change of variable
    # blade element id, also order of segments: {1, ..., 4} = {hL, rL, hR, rR}

    id = np.concatenate((np.ones((1, p.n_blade_prox)) * 1,np.ones((1, p.n_blade_dist)) * 2,np.ones((1, p.n_blade_prox)) * 3,np.ones((1, p.n_blade_dist)) * 4,np.ones((1, p.n_blade_tail)) * 5),axis=1).reshape(18,) #tail
    # Left wing segments positions
    # NOTE: assume symmetric wing.
    # Left wing center of pressure positions(local frame)
    yL_prox = np.arange(1 / 2 / p.n_blade_prox, 1, 1 / p.n_blade_prox) * p.span_prox
    xfL_prox = np.array(interp1d(p.wing_y_prox, p.wing_xf_prox)(yL_prox)).reshape(2, )
    xrL_prox = np.array(interp1d(p.wing_y_prox, p.wing_xr_prox)(yL_prox)).reshape(2, )

    yL_dist = np.arange(1 / 2 / p.n_blade_dist, 1, 1 / p.n_blade_dist) * p.span_dist
    xfL_dist = np.array(interp1d(p.wing_y_dist, p.wing_xf_dist)(yL_dist)).reshape(6, )
    xrL_dist = np.array(interp1d(p.wing_y_dist, p.wing_xr_dist)(yL_dist)).reshape(6, )
    
    dc_x_prox = (xfL_prox + xrL_prox) / 2 + (xfL_prox - xrL_prox) / 4
    dc_x_dist = (xfL_dist + xrL_dist) / 2 + (xfL_dist - xrL_dist) / 4
    dc_x = np.concatenate((dc_x_prox, dc_x_dist, dc_x_prox, dc_x_dist), axis=0)
    dc_y = np.concatenate((yL_prox, yL_dist, -yL_prox, -yL_dist), axis=0)

    dc_x_n = dc_x.reshape(1, 16)
    dc_y_n = dc_y.reshape(1, 16)
    dc_wing_1 = np.concatenate((dc_x_n, dc_y_n), axis=0)
    dc_wing = np.concatenate((dc_wing_1, dc_x_n * 0), axis=0)

    x_tail = p.tail_center[0] + p.tail_chord * 0.25 * np.ones((1, p.n_blade_tail))
    y_tail = p.tail_center[1]+ (np.arange(1 / 2 / p.n_blade_tail,1, 1 / p.n_blade_tail) - 0.5) * p.tail_width
    z_tail = p.tail_center[2] * np.ones((1, p.n_blade_tail))

    x_tail_n = x_tail.reshape(1,2)
    y_tail_n = y_tail.reshape(1,2)
    z_tail_n = z_tail.reshape(1,2)

    dc_tail = np.concatenate((x_tail_n, y_tail_n, z_tail_n),axis=0)
    dc=np.concatenate((dc_wing,dc_tail),axis=1).reshape(3,18)

    # Calculate element chord length
    c_prox = xfL_prox - xrL_prox
    c_dist = xfL_dist - xrL_dist
    c_prox_n=c_prox.reshape(1,2)
    c_dist_n=c_dist.reshape(1,6)

    c_wing=np.concatenate((c_prox_n, c_dist_n, c_prox_n, c_dist_n),axis=1)

    c_tail = np.ones((1, p.n_blade_tail)) * p.tail_chord
    c_wing_n=c_wing.reshape(1,16)
    c_tail_n=c_tail.reshape(1,2)

    c=np.concatenate((c_wing_n,c_tail_n),axis=1).reshape(18,)

    # Calculate theta
    yL_prox_n=yL_prox.reshape(1,2)
    yL_dist_n=yL_dist.reshape(1,6)
    yL = p.body_radius + np.concatenate((yL_prox_n,yL_dist_n + p.span_prox),axis=1)
    yL_n=yL.reshape(1,8)
    y = np.concatenate((yL_n, -yL_n),axis=1).reshape(16,)

    theta = np.arccos(y * 2 / p.span_max) # Not including tail, since it 's quasi-steady

    return id, dc, c, theta