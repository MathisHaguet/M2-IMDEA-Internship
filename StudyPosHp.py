""" 
    This python code, is used to show how the poisition of the source in the tube 
    impact the resonances. 
    Internship @ UdeS - Crash 
    Date   : 06/03/2023 to 28/07/2023
    Author : HAGUET Mathis 
    Supervisors : ROBIN Olivier and MELON Manuel
    
"""
#------------------------------------------------------------------------------
# Packages needed to run this code

import numpy as np 
import matplotlib.pyplot as plt
from ClassFunction import *
import os
# from tqdm import trange
path = os.getcwd()

#%%

class Radiation_HpShift():
    def __init__(self,f,LS,Duct,U,r):
        
        self.f = f
        self.w = 2*np.pi*f
        self.k = self.w/Air.c0
        self.LS = LS
        self.r = r
        self.Duct = Duct
        self.U = U
        s = 1j*self.w
        
        self.Qd = ((self.U*self.LS.Sd)/(self.LS.Bl*self.LS.Qec)) * ((s/self.LS.Wc) / ( (s)**2/(self.LS.Wc**2) + (s/(self.LS.Qtc*self.LS.Wc)) + 1 ))
        # self.ZRT = s * 0.61 * ((Air.rho)/(np.pi*self.Duct.Radius))
        self.ZRT = self.Duct.Z2 
        # self.ZLF = ((Air.rho*Air.c0*(self.k*self.Duct.Radius)**2)/(2*np.pi*self.Duct.Radius**2)) + 1j*((8*Air.rho*Air.c0*self.k*self.Duct.Radius)/(3*np.pi**2*self.Duct.Radius**2))
        # self.ZTf = s * ((Air.rho*self.Duct.Length)/(np.pi**2*self.Duct.Radius)) + self.ZRT
        
        Cm_tot = Parallel(self.LS.Cms , (self.LS.Cab / self.LS.Sd**2)) 
        Mm_tot = self.LS.Mms 
        
        self.Ze_tot = self.LS.Re +  s * self.LS.Le + (( s * Cm_tot * self.LS.Bl**2 ) \
                                / (( s**2 * Mm_tot * Cm_tot) + (s * self.LS.Rms * Cm_tot + 1)))
        
        self.Ze_coil = self.LS.Re + s*self.LS.Le

        self.P = (self.LS.Bl * self.U) / (self.Ze_coil * self.LS.Sd)
        
        self.Zae = (self.LS.Bl**2 / self.Ze_coil ) / self.LS.Sd**2             # Electrical Imp coil -> Acoutic Domain
        self.Zas = self.LS.Zms/(self.LS.Sd**2)                                 # Mechanical Imp  -> Acoutic Domain
        self.Za0 = self.Duct.Imp                                               # Duct Imp -> Acoustic Domain
        self.Zab = 1/(s*self.LS.Cab)                                           # Imp at Back -> Acoustic Domain
        self.Zaf = np.real(self.ZRT)                                           # Imp at the Front -> Acoustic Domain
        
        self.Q0 = self.P / (self.Zae  + self.Zas + self.Zab + self.Za0 + self.ZRT)
        self.Qv = self.Q0 / (np.cos(self.k*(self.Duct.Length )) + 1j*np.sin(self.k*(self.Duct.Length ))/self.Duct.Zc) * ADMITANCE 
        
        
        self.Prear = 1j*self.k*Air.rho*Air.c0 * self.Qd * np.exp(-1j*self.k*(self.r)) / (4*np.pi*(self.r))
        
        self.Pvent = 1j*self.k*Air.rho*Air.c0 * self.Qv * np.exp(-1j*self.k*(self.r)) / (4*np.pi*(self.r))
        # self.Ptot =  1j*self.k*Air.rho*Air.c0 * (self.Qd * np.exp(-1j*self.k*(self.r))/2*np.pi*(self.r) - self.Qv * np.exp(-1j*self.k*self.r)/2*np.pi*self.r )



#%%

path = os.getcwd()
Data = np.loadtxt(path + "/DATA/AKABAK SIMU/REC_PosCompa.txt")
f = Data[:,0]
# f = DATA[:,0]
# f  = np.logspace(np.log10(20),np.log10(20e3),DSP.fs)

TB_W4_1658SB = Loudspeaker(f ,  61.1 , 4.2 , 0.21e-3 , 7.08 , 12.97e-3, 10.4, 5e-2 , 1 , BackVolume=8)
P_SDS_135F25CP02_04 = Loudspeaker(f,57.7 , 3.3, 1.22e-3, 7.39, 20.6e-3, 3.89, PouceEnMetre(5.5)/2, 1, BackVolume=9)

# Delta_D = P_SDS_135F25CP02_04.rd * 1
Delta_D = TB_W4_1658SB.rd * 1

# Delta_D = -0.02
Shift = 1/3
LenTotDuct = 1.4 + Delta_D

Driver_pos = LenTotDuct*Shift
LenRightDuct = LenTotDuct - Driver_pos #+ Delta_D
LenLeftDuct = Driver_pos 

REC_Tot = RectangularTube(f, LenTotDuct, 0.1, 0.1, 0.01, DuctType="Open" )
REC_Left = RectangularTube(f, LenLeftDuct, 0.1, 0.1, 0.01, DuctType="Closed" )
REC_Right = RectangularTube(f, LenRightDuct, 0.1, 0.1, 0.01, DuctType="Open" )

TUBES = [REC_Left,REC_Right]
Duct = Connection(f, *TUBES)

TEST = Parallel(REC_Left.Imp,REC_Right.Imp )

global ADMITANCE
ADMITANCE = (REC_Left.Imp)/(REC_Right.Imp + REC_Left.Imp)

U = 1
r = 1 

RAD = Radiation_HpShift(f, TB_W4_1658SB , Duct , U , r )
PYTHON = Dict(RAD.f, RAD.Pvent, label=r"Analytic model : Speaker @  L/3",color="k")
RAD_TOT = Radiation(f, TB_W4_1658SB , REC_Tot , U , r)
TOT = Dict(RAD_TOT.f, RAD_TOT.Pvent, label="Analytic model : Speaker @ 0",color="k",linestyle="--")

path = os.getcwd()
Data = np.loadtxt(path + "/DATA/AKABAK SIMU/REC_PosCompa.txt")
AKABAK = Dict(Data[:,0],Data[:,7] + 1j*Data[:,8], label="AKABAK : Speaker @ L/3",color="r")

PlotFRF(TOT,PYTHON,AKABAK  , Title="" , linewidth=3,loc="upper left")

# PlotFRF(AKABAK_DUCT,TOT , Title="" , linewidth=3,loc="upper left")

# plt.figure(facecolor='white')

# # plt.semilogx(f , 20*np.log10(abs(RAD.Za0 )/2e-5) , label="Python 1/3 Duct")
# # plt.semilogx(f , 20*np.log10(abs(Duct.Imp)/2e-5) , label="Python End")
# # plt.semilogx(f , 20*np.log10(abs(TEST)/2e-5) , label="Python 1/3 Test")
# plt.semilogx(f , 20*np.log10(abs(RAD.Pvent)/2e-5) , label="Python 1/3" , linewidth=3)
# plt.semilogx(f , 20*np.log10(abs(RAD_TOT.Pvent)/2e-5) , label="Python end", linewidth=3)
# plt.semilogx(f , 20*np.log10(abs(Data[:,7] + 1j*Data[:,8])/2e-5) , label="Akabak 1/3", linewidth=3)
# plt.grid(which="both")
# plt.xlim([20,1000])
# plt.legend(loc="upper left")

#%%
TUBES = []

Length = 1.4 
NbrSlice = 200

D = 0.1*3/4
W = 0.1
H = 0.1
Width = VoigtSection(H, D, W, Length, NbrSlice)


LengthSlice = Length/NbrSlice
REC_Tot = CylindricalTube(f, Length , 0.05 , 0.02 , DuctType="Open") #RectangularTube(f, Length, 0.1, 0.1, 0.01, DuctType="Open" )

for i in range(NbrSlice):
    TUBES.append(RectangularTube(f, LengthSlice, 0.1 , Width[i] , 0.02 , DuctType="Open"))
    
    # TUBES.append(CylindricalTube(f, LengthSlice, Width[i]/2 , 0.02,DuctType="Open"))

plt.figure(facecolor='white')
plt.plot(Width , "o-")
plt.ylabel("Width (m)")
plt.xlabel("Tube Number")
plt.grid(which="both")
plt.show()

CONNECT = Connection(f, *TUBES)
ImpTest = CONNECT.T[0,1,:] / CONNECT.T[1,1,:]

Diff = abs(20*np.log10(np.abs(CONNECT.ConnectImp)/2e-5)-20*np.log10(np.abs(REC_Tot.Imp)/2e-5))
plt.figure(facecolor='white')
plt.semilogx(f,20*np.log10(np.abs(CONNECT.ConnectImp)/2e-5) , label="Connected Tubes")
plt.semilogx(f,20*np.log10(np.abs(REC_Tot.Imp)/2e-5) , label="1 Tube")
# plt.semilogx(f,Diff)
plt.legend()
plt.grid(which="both")
plt.xlim([20,1000])
plt.show()

#%%
# plt.close('all')

RAD_TOT = Radiation(f, TB_W4_1658SB , REC_Tot , U=1 , r=1)
RAD_CONNECT = Radiation(f, TB_W4_1658SB , CONNECT , U=1 , r=1)

ind = 7
DATA = np.loadtxt("C:/Users/Acer/Documents/STAGE - UdeS/CODES/AKABAK SIMU/COMPA SLOPE VOIGT/RecSlope_HpEnd.txt")
REC_TAPERED_AKABAK = Dict(DATA[:,0], DATA[:,ind] + 1j*DATA[:,ind+1] , "AKABAK Tapered")

# RAD_TOT = Radiation(f, TB_W4_1658SB , REC_Tot , U , r)
Connect = Dict(RAD_CONNECT.f, RAD_CONNECT.Pvent , "Py Tapered : {} element(s)".format(NbrSlice))
One = Dict(RAD_TOT.f, RAD_TOT.Pvent , "Py Straight : 1 element")

PlotFRF(One,Connect,REC_TAPERED_AKABAK, Title="Connect function conparison")

# plt.figure()
# for i in trange(NbrSlice):
#     plt.semilogx(f,20*np.log10(np.abs(TUBES[i].Imp)))
# plt.grid(which="both")
# plt.xlim([20,1000])
# plt.show()


#%%

def fbg(l,S,V):
    f = Air.c0/(2*np.pi) * np.sqrt((S)/(l*V))
    return f

#%%


    