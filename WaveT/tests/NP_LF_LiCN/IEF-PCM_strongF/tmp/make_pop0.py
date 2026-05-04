import numpy as np
a=np.loadtxt("c_t_1.dat", usecols=(2,3))
b=a[:,0]+1j*a[:,1]
pop2=np.sqrt(b*np.conj(b))
pop1=b*np.conj(b)
np.savetxt("pop0_sqrt.dat", np.real(pop2))
np.savetxt("pop0.dat", np.real(pop1))

