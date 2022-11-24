import numpy as np
from scipy import signal


def phase_HT(data):
    """
    Reconstruct the phase from observed oscillatory signal.

    Parameters
    ----------
        data: ndarray with shape (1, N)
            Observed oscillatory signal.
    
    Returns
    -------
        phase_h: ndarray with shape (1, N)
            Reconstructed phase signal.
    """
    analytic_signal = signal.hilbert(data)

    phase_h = np.angle(analytic_signal[:, 0]).reshape((1, 1)) + np.append(np.zeros((1,1)), np.cumsum(np.angle(analytic_signal.T[1:] / analytic_signal.T[:-1]), axis=0).T, axis=1) 
    phase_h = phase_h - 2*np.pi*round(phase_h[0][0] / (2*np.pi))
    return phase_h


def reconst_mod_fft(ft_mod, m):
    """
    Reconstruct the true phase-modulation from phase-modulation reconstructed via the HT.

    Parameters
    ----------
        ft_mod: ndarray with shape (1,N)
            Fourier transform of phase-modulation. 
        m: int
            Index coresponding to the effective frenqency of the oscillatory signal.
    
    Returns
    -------
        ft_mod_new: ndarray with shape (1,N)
            Reconstructed Fourier transform of phase-modulation.
    """
    N = np.size(ft_mod, axis=1)
    ft_mod_new = (np.zeros(N) + 1j*np.zeros(N)).reshape(1, N)
    Nh = int(N/2)
    for n in np.arange(2*m)+Nh-2*m+1:
        ft_mod_new[:, n] = 2*ft_mod[:, n]
    for n in np.arange(Nh - 3*m)[::-1] + m+1:
        ft_mod_new[:, n] = 2*ft_mod[:, n] + ft_mod_new[:, n+2*m]
    ft_mod_new[:, m] = 2*np.real(ft_mod[:, m]+ 0.5*ft_mod_new[:, 3*m]) + 1j*np.imag(ft_mod[:, m]+ 0.5*ft_mod_new[:, 3*m])
    for n in np.arange(m)[::-1]:
        ft_mod_new[:, n] = ft_mod[:, n] + 0.5 * np.conjugate(ft_mod_new[:, 2*m - n]) + 0.5*ft_mod_new[:, n+2*m]
    for n in np.arange(Nh+1, N):
        ft_mod_new[:, n] = np.conjugate(ft_mod_new[:, N-n])
    
    return ft_mod_new

def phase_postprocess(phase, return_outliers=False):
    """
    Detect outliers in the phase difference and replace the phase at that point with linear interpolation.

    Parameters
    ----------
        phase: ndarray with shape (1, N)
            Phase signal.
        return_outliers: bool, optional
            Whether or not return the list of detected outliers.

    Returns
    ----------
        phase_new: ndarray with shape (1, N)
            Postprocessed phase signal.
        outlier_list: ndarray
            List of outliers.

    """
    K = np.size(phase, axis=0) 
    N = np.size(phase, axis=1) #Length of the data
    phase_new = np.array([])
    outlier_list = np.array([])
    for k in range(K):
        diff = phase[k][1:] - phase[k][:-1]
        diff_mad = np.sum(np.abs(diff - np.median(diff)))/  (N-1)
        outlier = np.abs(diff - np.median(diff)) > 3*diff_mad
        outlier_list = np.append(outlier_list.reshape(k-1, N-1), outlier.reshape(1, N-1), axis=0)
        phase_new_k = np.zeros(N)
        n=0
        while (n<N-1): 
            if outlier[n] ==True:
                J=1
                while(n+J < N-2 and outlier[n+J] ==True): 
                    J+=1 
                for j in range(J):
                    if (n==0 or n+J-1==N-1): phase_new_k[n+j] = phase[k][n+j]
                    else:
                        phase_new_k[n+j] = (phase[k][n-1]* (J-j)+ phase[k][n+J]*(j+1)) / (J+1)
                n += J 
            else: 
                phase_new_k[n] = phase[k][n]
                n+=1
        phase_new_k[N-1] = phase[k][N-1]
        phase_new_k = phase_new_k.reshape(1, N)
        phase_new = np.append(phase_new.reshape(len(phase_new), N), phase_new_k, axis=0)
    if return_outliers==True: return phase_new, outlier_list
    else: return phase_new

def phase_reconst(data, tau, max_cut_len=None):
    """
    Extended Hilbert transform
    
    Parameters
    ----------
        data:
            ndarray with size of (1, N)
            Oscillatory signal.
            Should be preprocessed to smoothly connect from the end point to the start point, which guarantees the accuracy of Fourier tansform.
        tau: float
            Time step size. 
        max_cut_len: int
            Maximum length of cutting. 

    Returns
    -------
        phase_new: ndarray with size of (1, N)
        Phase 
    """

    phase_h = phase_HT(data)

    #Cut the data so that the initial point of phase_h becomes zero.
    if max_cut_len==None:  
        freq_tmp = np.mean(phase_h[0][1:] - phase_h[0][:-1])/tau/(2*np.pi)
        max_cut_len=2*round(2*np.pi / freq_tmp / tau ) 
    initindx = np.argmin( np.abs( (phase_h[:, :max_cut_len] + np.pi) % (2*np.pi) -np.pi))
    phase_h = phase_h[:, initindx:]
    phase_h = phase_h - 2*np.pi*round(phase_h[0][0]/(2*np.pi)) # set initial phase at zero

    N = np.size(phase_h, axis=1)
    tspan = np.linspace(0, N*tau-tau, N)

    efffreq = np.mean(phase_h[0][1:] - phase_h[0][:-1])/tau/(2*np.pi) #Effective frequency
    m = np.rint(efffreq * (N*tau)).astype(np.int64) 
    
    efffreq_trend = 2*np.pi*efffreq* tspan.reshape(1, N) #Linear trend of the phase 

    phase_mod = phase_h - phase_h[0][0] - efffreq_trend #Phase-modulation
    ft_mod = np.fft.fft(phase_mod)
    ft_mod_new = reconst_mod_fft(ft_mod, m)

    phase_mod_new = (np.fft.ifft(ft_mod_new)).reshape(1,N)

    phase_new = np.real((phase_h[0][0]).reshape(1,1) + efffreq_trend + phase_mod_new)
    phase_new = phase_new - 2*np.pi*round(phase_new[0][0]/(2*np.pi))
    phase_new = phase_postprocess(phase_new)

    return phase_new, phase_h
        
