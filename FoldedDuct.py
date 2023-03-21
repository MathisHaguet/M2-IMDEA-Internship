import numpy as np
import matplotlib.pyplot as plt
from ClassFunction import * 

#%%
f  = np.logspace(np.log10(20),np.log10(5e3),40000)

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

#%%

TubeREC = RectangularTube(f ,1.7, 0.1, 0.1, 0.01 ,DuctType="Open" )
TubeCYL = CylindricalTube(f ,1.7, 0.05, 0.01 ,DuctType="Open")

ZREC = TubeREC.Impedance(f)
ZCYL = TubeCYL.Impedance(f)

plt.figure()
plt.semilogx(f,10*np.log10(np.abs(ZREC)) , label="Z REC")
plt.semilogx(f,10*np.log10(np.abs(ZCYL)) , label="Z CYL")

plt.grid(which="both")
plt.legend()
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

DATA = np.loadtxt("/Users/mathishaguet/Documents/STAGE - UdeS/CODES/HP_OFFSET.txt")


StraightOpenDuct = RectangularTube(f,1.7, 0.1, 0.1, 0.01 , DuctType="Open")
OpenDuct = RectangularTube(f,1.7*2/3, 0.1, 0.1, 0.01 , DuctType="Open")
ClosedDuct = RectangularTube(f,1.7/3, 0.1, 0.1, 0.01 , DuctType="Closed")
ZT = Parallel(OpenDuct.Impedance(f), ClosedDuct.Impedance(f))

plt.figure()
# plt.semilogx(DATA[:,0],DATA[:,1] , label="AKABAK")
# plt.semilogx(DATA[:,0],DATA[:,2] , label="AKABAK")

plt.semilogx(f,10*np.log10(np.abs(StraightOpenDuct.Impedance(f))) , label="Straight Open Duct")
# plt.semilogx(f,10*np.log10(np.abs(OpenDuct.Impedance(f))) , label="Open Duct")
# plt.semilogx(f,10*np.log10(np.abs(ClosedDuct.Impedance(f))) , label="Closed Duct")
plt.semilogx(f,10*np.log10(np.abs(ZT)) ,"k", label="Parallel Duct")

plt.legend()
plt.grid(which="both")
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