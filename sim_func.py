import numpy as np
import sdeint
import process_observed_data as proc
import importlib
importlib.reload(proc)

def u_QuasiPeriodic(tspan, omega, b):
    N = len(tspan)
    u = b*( np.sin(np.sqrt(2)*omega*tspan) + np.cos(np.sqrt(3)*omega*tspan)).reshape(1, N)
    return u

def G_whitenoise(x,t, sigma):
    num_node = len(x)
    return np.diag(np.ones(num_node)* sigma)

def u_OU(tspan, k, sigma):
    N = len(tspan)
    def OU_det(y, t):
        return - k*y
    def OU_stc(y, t):
        return G_whitenoise(y,t, sigma)
    u = (sdeint.itoint(f=OU_det, G=OU_stc, y0=np.array([0.0]), tspan=tspan)).reshape(1, N)
    return u

def u_PhaseShift(tspan1, tspan2, tspan3, omega1, omega2, sigma):
    def det1(y,t):
        return np.array([0.0])
    def det2(y,t):
        return np.array([omega2-omega1])
    def stc(y, t):
        return G_whitenoise(y, t, sigma)
    u1 = (sdeint.itoint(f=det1, G=stc, y0=np.array([0.0]), tspan=tspan1))
    u2 = (sdeint.itoint(f=det2, G=stc, y0=u1[-1], tspan=tspan2))
    u3 = (sdeint.itoint(f=det1, G=stc, y0=u2[-1], tspan=tspan3))
    u = np.append(np.append(u1[:-1], u2[:-1]), u3)
    u = u.reshape(1,len(u))
    return u

def phase_mod_QuasiPeriodic(tspan, omega, b):
    N= len(tspan)
    phase_linear = (omega * tspan).reshape(1,N)
    u = u_QuasiPeriodic(tspan, omega, b)
    phase = (phase_linear + u).reshape(1,N)
    return phase

def phase_mod_OU(tspan, omega, k, sigma):
    N=len(tspan)
    phase_linear = (omega * tspan) .reshape(1,N)
    u = u_OU(tspan, k, sigma)
    phase = (phase_linear + u ).reshape(1,N)
    return phase

def phase_mod_PhaseShift(tspan1, tspan2, tspan3, omega1, omega2, sigma):
    tspan = np.append(np.append(tspan1[:-1], tspan2[:-1]), tspan3)
    N=len(tspan)
    phase_linear = (omega1 * tspan) .reshape(1,N)
    u = u_PhaseShift(tspan1, tspan2, tspan3, omega1, omega2, sigma)
    phase = (phase_linear + u ).reshape(1,N)
    return phase

def x_QuasiPeriodic(tspan, omega, b, max_cut_len=None):
    if max_cut_len==None: 
        tau = tspan[1]-tspan[0] 
        max_cut_len=2*round(2*np.pi / omega / tau ) 
    phase = phase_mod_QuasiPeriodic(tspan, omega, b)
    x_tmp = np.cos(phase)
    x, indx_s, indx_e = proc.data_preprocess(x_tmp, max_cut_len, return_indx=True)
    phase = phase[:, indx_s:indx_e]
    return phase, x


def x_OU(tspan, omega, k, sigma, max_cut_len=None):
    if max_cut_len==None: 
        tau = tspan[1]-tspan[0] 
        max_cut_len=2*round(2*np.pi / omega / tau ) 
    phase = phase_mod_OU(tspan, omega, k, sigma)
    x_tmp = np.cos(phase)
    x, indx_s, indx_e = proc.data_preprocess(x_tmp, max_cut_len, return_indx=True)
    phase = phase[:, indx_s:indx_e]
    return phase, x

def x_ampmod_QuasiPeriodic(tspan, omega, b, r, nu, max_cut_len=None):
    if max_cut_len==None: 
        tau = tspan[1]-tspan[0] 
        max_cut_len=2*round(2*np.pi / omega / tau ) 
    phase = phase_mod_QuasiPeriodic(tspan, omega, b)
    amp = 1 + r*np.cos(nu * tspan)
    x_tmp = amp * np.cos(phase)
    x, indx_s, indx_e = proc.data_preprocess(x_tmp, max_cut_len, return_indx=True)
    phase = phase[:, indx_s:indx_e]
    return phase, x

def x_PhaseShift(tspan1, tspan2, tspan3, omega1, omega2, sigma, max_cut_len=None):
    tspan = np.append(np.append(tspan1, tspan2), tspan3)
    N=len(tspan)
    if max_cut_len==None: 
        tau = tspan[1]-tspan[0] 
        max_cut_len=2*round(2*np.pi / omega1 / tau ) 
    phase = phase_mod_PhaseShift(tspan1, tspan2, tspan3, omega1, omega2, sigma)
    x_tmp = np.cos(phase)
    x, indx_s, indx_e = proc.data_preprocess(x_tmp, max_cut_len, return_indx=True)
    phase = phase[:, indx_s:indx_e]
    return phase, x