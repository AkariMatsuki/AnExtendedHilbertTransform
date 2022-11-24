import numpy as np
import Extended_Hilbert_transform as ext_ht
import sys

file_name = sys.argv[1]
tau = float(sys.argv[2])

x = np.loadtxt(file_name)
x = x.reshape(1, len(x))

phase_exth, phase_h = ext_ht.phase_reconst(x, tau)

np.savetxt("phase_exth.txt", phase_exth)
np.savetxt("phase_h.txt", phase_h)
