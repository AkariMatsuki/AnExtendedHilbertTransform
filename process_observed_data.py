import numpy as np

def data_preprocess(data, max_cut_len, return_indx=False):
    """
    Preprocessing the observed signal. Cut the data so that the initial and final phase becomes zero (mod 2pi).

    Parameters
    ----------
        data: ndarray with shape(1,N)
            An oscillatory signal. 
        max_cut_len: int
            Maximum length of cutting. 
    
    Returns
    -------
        output_data: ndarray with (1, N)
            Preprocessed data.
        indx_s: int
            Index of the original data corresponding to the initial point of the preprocessed data.
        indx_e: int
            Index of the original data corresponding to the final point of the preprocessed data.
    """
    N = np.size(data, axis=1)
    indx_s = np.argmax(data[0][:max_cut_len]) 
    indx_e = N-max_cut_len + np.argmax(data[0][-max_cut_len:])
    N_cut = indx_e-indx_s 
    output_data = (data[:, indx_s:indx_e]).reshape(1,N_cut)
    if return_indx==True: return output_data, indx_s, indx_e
    else: return output_data
