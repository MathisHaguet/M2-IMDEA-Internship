""" 
    This python code, calculate the loss due to air inside the tube. 
    Internship @ UdeS - Crash 
    Date   : 06/03/2023 to 28/07/2023
    Author : HAGUET Mathis 
    Supervisors : ROBIN Olivier and MELON Manuel
    
"""
#------------------------------------------------------------------------------
# Packages needed to run this code

import matplotlib.pyplot as plt 
import numpy as np
# from tqdm import trange

#%%
class Air:
    c0 = 343.318219 # sound celerity [m.s-1]
    rho = 1.204786  # air density [kg.m-3]
    mu = 1.6e-5     # ideal gaz viscosity [kg.m-1.s-1] 
    kappa = 0.02225  # thermal conductivity [W.m-1.K-1]
    Cp = 1006       # spcific heat @ constant pressure [J.kg-1.K-1]
    T0 = 293.15     # temperature [K]
    P0 = 101325     # atmospheric pressure [Pa]
    Cv = Cp - (P0/(rho*T0)) # specific heat per unit mass
    nu = mu/rho  
    nup = kappa/(rho*Cv)
    gamma = Cp/Cv 
    
def LossRecTube(f,Width,Height):
    w = 2 * np.pi * f # pulsation [rad.s-1]
    Nit = 200          # number of iteration
    a = Width/2        
    b = Height/2
    
    temp_rho = 0
    temp_k   = 0 
    
    RHO = 0 
    K   = 0 
    
    # SUM 1 (m)
    for m in range(Nit):
        alpha = (m + 0.5)*(np.pi/a)
        # print(m)
        # SUM 2 (n)

        for n in range(Nit):
            beta  = (n + 0.5)*(np.pi/b)
            # print(n)
            temp_rho += 1/(alpha**2 * beta**2 * ( alpha**2 + beta**2 + (1j*w/Air.nu)))
            temp_k   += 1/(alpha**2 * beta**2 * (alpha**2 + beta**2 + (1j*w*Air.gamma/Air.nup)))

    RHO = Air.rho * ((Air.nu * a**2 * b**2) / (4 * 1j *w)) * 1/temp_rho
    B = (4*1j*w *(Air.gamma-1)) / (Air.nup * a**2 * b**2)
    K   = 1/Air.P0 * ( 1 - B * temp_k)
    
    # c_loss = np.sqrt(1/K/RHO)
    c_loss = np.sqrt(1/(K*RHO))

    
    
    return RHO , K , c_loss


#%%

f = np.linspace(20, 5000, 2000)
RHO , K , c_loss = LossRecTube(f, 0.1, 0.1)


#%%
fig, ax = plt.subplots(ncols=3 , facecolor='white' , figsize=(12,4))

ax[0].semilogx(f,abs(c_loss) , 'k' , linewidth="2")
ax[0].set_ylabel(r"$|c(\omega)|$ [m/s]", fontsize=20)
ax[0].set_xlabel("frequency [Hz]", fontsize=20)
# ax[0,0].set_ylim((330 , 344))
ax[0].grid()

ax[1].semilogx(f,abs(RHO), 'k' , linewidth="2")
ax[1].set_ylabel(r"$|\rho(\omega)|$ [kg/m$^3$]", fontsize=20)
ax[1].set_xlabel("frequency [Hz]", fontsize=20)
ax[1].grid()

ax[2].semilogx(f,abs(K), 'k' , linewidth="2")
ax[2].set_ylabel(r"$|K(\omega)|$ [N/m]", fontsize=20)
ax[2].set_xlabel("frequency [Hz]" , fontsize=20)
ax[2].grid()



plt.tight_layout()
plt.show()

print(abs(c_loss[-1]))


