import numpy as np
import matplotlib.pyplot as plt
from ClassFunction import * 

#%%
# f  = np.logspace(np.log10(20),np.log10(5e3),40000)

DATA_CYL = np.loadtxt("/Users/mathishaguet/Documents/STAGE - UdeS/CODES/CYL_REAL_IMAG_SPECIFIC_IMP.txt") # REAL & IMAG
DATA_REC = np.loadtxt("/Users/mathishaguet/Documents/STAGE - UdeS/CODES/REC_REAL_IMAG_SPECIFIC_IMP.txt") # REAL & IMAG
DATA_CYL_THIRD = np.loadtxt("/Users/mathishaguet/Documents/STAGE - UdeS/CODES/CYL_REAL_IMAG_SPECIFIC_IMP_HP_THIRD.txt") # REAL & IMAG

f = DATA_CYL[:,0]
#==================
Fr = 50
WidthHole = 0.05
MaxWidth  = 0.15
MaxHeight = 0.15
MaxDepth  = 0.20
#==================

NbrBend , LengthBend , WidthBend , HeightBend  = BendCalculation(QWL(Fr), MaxWidth, MaxHeight, MaxDepth)

print("Resonance Frequency : {} Hz".format(Fr))
print("Bend Number : {}".format(NbrBend)) 
print("Length Bend : {} m".format(LengthBend))
print("Width Bend  : {} m".format(WidthBend))
print("Height Bend : {} m ".format(HeightBend))

#%% OPEN/CLOSED
Fr = 50
Length = QWL(Fr)
Width = 0.1
Height = 0.1

REC = RectangularTube(f ,1.7, 0.1, 0.1, 0.01 ,DuctType="Open" )
CYL = CylindricalTube(f, 1.7, 0.1, 0.01,DuctType="Open")

fig , ax = plt.subplots(nrows=2)

ax[0].semilogx(f,20*np.log10(np.abs(REC.Imp)) ,"k", label="REC")
ax[0].semilogx(f,20*np.log10(np.abs(CYL.Imp)) ,"b", label="CYL")
ax[0].set_ylabel(r'$\mathbb{R}(Z_{in})$ [dB]' ,fontsize=13)

ax[1].semilogx(f,np.angle(REC.Imp) ,"k", label="REC")
ax[1].semilogx(f,np.angle(CYL.Imp) ,"b", label="CYL")
ax[1].set_ylabel(r'$\mathbb{Im}(Z_{in})$ [dB]' ,fontsize=13)

for i in range(2):
    ax[i].grid(which="both")
    ax[i].legend()
    ax[i].set_xlabel('Frequency [Hz]',fontsize=13)
    
plt.tight_layout()
plt.show()

#%%

Tube = RectangularTube(f,1.7, 0.1, 0.1, 0.01 ,DuctType="Open" )
Z1 = Tube.Impedance(f)
Z2 = Tube.RadImp(f)
Ztest = (1j*2*np.pi*f*Air.rho*1.715)/0.01
plt.figure()
plt.semilogx(f,10*np.log10(np.abs(Z1)) , label="Z1")
plt.semilogx(f,10*np.log10(np.abs(Z2)) , label="Z2")
plt.semilogx(f,10*np.log10(np.abs(Ztest)) , label="Z2")

plt.grid(which="both")
plt.legend()
plt.show()

#%%
TUBES = []
ind = []

for i in range(NbrBend):
    TUBES.append(RectangularTube(f,LengthBend-WidthHole/2, WidthBend , HeightBend , 0.001 , DuctType="Open"))

if NbrBend > 1:
    for i in range(NbrBend):
        ind.append(i)
        if i %2 ==0:
            ind.remove(i)
        # print(ind)
    for i in ind:
        # print(i)
        _ , _ , LBp = BendConnection(TUBES[i-1], TUBES[i+1], WidthHole)
        TUBES.insert(i , RectangularTube(f,LBp, WidthHole , HeightBend , 0.001 , DuctType="Open"))
else:
    TUBES = TUBES

# Tube1 = RectangularTube(1/2, 0.1, 0.15, 0.01)
# Tube2 = RectangularTube(1.7/2, 0.1, 0.1, 0.01)

Test = Connection(f, *TUBES)
ImpTest = Test.Impedance(f)

# %%
# plt.close('all')

plt.figure()
plt.semilogx(f,10*np.log10(np.abs(ImpTest)) , label="Connected Tubes")
plt.legend()
plt.grid(which="both")
plt.show()

plt.figure()
for i in range(2):
    plt.semilogx(f,10*np.log10(np.abs(TUBES[i].Impedance(f))) , label="Tubes {}".format(i+1))
plt.grid(which="both")
plt.legend()
plt.show()

#%%
plt.close("all")


AKABAK_REC = DATA_REC[:,1] + 1j*DATA_REC[:,2]
AKABAK_CYL = DATA_CYL[:,1] + 1j*DATA_CYL[:,2]
AKABAK_CYL_THIRD = DATA_CYL_THIRD[:,1] + 1j*DATA_CYL_THIRD[:,2]

REC = RectangularTube(f, 1.715, 0.1, 0.1, 0.01, DuctType="Open" )
CYL = CylindricalTube(f, 1.715, 0.1, 0.01, DuctType="Open")

# OpenDuct = RectangularTube(f,1.715*2/3, 0.1, 0.1, 0.01 , DuctType="Open")
# ClosedDuct = RectangularTube(f,1.715/3, 0.1, 0.1, 0.01 , DuctType="Closed")
# ZT = Parallel(OpenDuct.Impedance(f), ClosedDuct.Impedance(f))

OpenDuct = CylindricalTube(f, 1.715*2/3, 0.1, 0.01,DuctType="Open")
ClosedDuct = CylindricalTube(f, 1.715/3, 0.1, 0.01,DuctType="Closed")
ZT = Parallel( ClosedDuct.Impedance(f),OpenDuct.Impedance(f))

#==============================================================================

fig , ax = plt.subplots(nrows=2 , sharex=True)

ax[0].semilogx(f,20*np.log10(np.abs(ZT)) ,"k" , label = "HP OFFSET 1/3")
ax[0].semilogx(f,20*np.log10(np.abs(CYL.Imp)) , label = "HP WO OFFSET")
ax[0].semilogx(f,20*np.log10(np.abs(AKABAK_CYL_THIRD)/CYL.Section_In) , label = "HP OFFSET 1/3 AKABAK")
ax[0].set_ylabel(r'$|Z_{in}|$ [dB]' ,fontsize=13)


ax[1].semilogx(f,np.angle(ZT) ,"k" , label = "HP OFFSET 1/3")
ax[1].semilogx(f,np.angle(REC.Imp) , label = "HP WO OFFSET")
ax[1].semilogx(f,np.angle(AKABAK_CYL_THIRD) , label = "HP OFFSET 1/3 AKABAK")
ax[1].set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r"$-\pi$" ,r"$-\frac{\pi}{2}$", r"$0$" ,r"$\frac{\pi}{2}$",r"$\pi$" ])
ax[1].set_ylabel(r'$\mathbb{\phi}(Z_{in})$ [rad]' ,fontsize=13)

for i in range(2):
    ax[i].grid(which="both")
    ax[i].legend()
    ax[i].set_xlabel('Frequency [Hz]',fontsize=13)
    
ax[0].set_title("Rectangular Duct")
plt.tight_layout()
plt.show()

#==============================================================================

fig , ax = plt.subplots(nrows=2 , sharex=True)
ax[0].semilogx(f,20*np.log10(np.abs(REC.Imp)) ,"k", label="PYTHON")
ax[0].semilogx(f,20*np.log10(np.abs(AKABAK_REC)/REC.Section_In) ,"b", label="AKABAK")
ax[0].set_ylabel(r'$|Z_{in}|$ [dB]' ,fontsize=13)

ax[1].semilogx(f,np.angle(REC.Imp) ,"k", label="PYTHON")
ax[1].semilogx(f,np.angle(AKABAK_REC),"b", label="AKABAK")
ax[1].set_ylabel(r'$\mathbb{\phi}(Z_{in})$ [rad]' ,fontsize=13)
ax[1].set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r"$-\pi$" ,r"$-\frac{\pi}{2}$", r"$0$" ,r"$\frac{\pi}{2}$",r"$\pi$" ])

for i in range(2):
    ax[i].grid(which="both")
    ax[i].legend()
    ax[i].set_xlabel('Frequency [Hz]',fontsize=13)
    
ax[0].set_title("Rectangular Duct")
plt.tight_layout()
plt.show()

#==============================================================================

fig , ax = plt.subplots(nrows=2 , sharex=True)
ax[0].semilogx(f,20*np.log10(np.abs(CYL.Imp)) ,"k", label="PYTHON")
ax[0].semilogx(f,20*np.log10(np.abs(AKABAK_CYL)/CYL.Section_In) ,"b", label="AKABAK")
ax[0].set_ylabel(r'$|Z_{in}|$ [dB]' ,fontsize=13)

ax[1].semilogx(f,np.angle(CYL.Imp) ,"k", label="PYTHON")
ax[1].semilogx(f,np.angle(AKABAK_CYL),"b", label="AKABAK")
ax[1].set_ylabel(r'$\mathbb{\phi}(Z_{in})$ [rad]' ,fontsize=13)
ax[1].set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r"$-\pi$" ,r"$-\frac{\pi}{2}$", r"$0$" ,r"$\frac{\pi}{2}$",r"$\pi$" ])
for i in range(2):
    ax[i].grid(which="both")
    ax[i].legend()
    ax[i].set_xlabel('Frequency [Hz]',fontsize=13)
    
ax[0].set_title("Cylindrical Duct")
plt.tight_layout()
plt.show()


#%%
Height = 0.1
Width = 0.1
Length = 1.715
D = Width/2 * 1/3
NbrPoints = 10

S_Voight = VoightSection(Height, D, Width, Length, NbrPoints)

plt.figure()
plt.plot(S_Voight)
plt.ylabel(r"Section [m$^2$]")
plt.xlabel("Division Number of the Pipe")
plt.grid()

TUBES = []
for s in S_Voight:
    TUBES.append(RectangularTube(f,Length/NbrPoints, np.sqrt(s), np.sqrt(s), 0.01 , Section_In=s , DuctType="Open"))

Voight_Test = Connection(f, *TUBES)

#%%

DATA = np.loadtxt("//Users/mathishaguet/Documents/STAGE - UdeS/CODES/VOIGHT_1_Third.txt")

plt.close("all")
plt.figure()
# plt.semilogx(DATA[:,0],DATA[:,2] , label="AKABAK")
plt.semilogx(f,20*np.log10(np.abs(Voight_Test.Impedance(f))) , label="Connected Tubes")
# plt.semilogx(f,20*np.log10(np.abs(StraightOpenDuct.Impedance(f))) , label="Straight Open Duct")
plt.legend()
plt.grid(which="both")
plt.show()

# plt.figure()
# for i in [0,10,20,30,40]:
#     plt.semilogx(f,10*np.log10(np.abs(TUBES[i].Impedance(f))) , label="Tubes {}".format(i+1))
# plt.grid(which="both")
# plt.legend()
# plt.show()