import numpy as np
import matplotlib.pyplot as plt

#%% 


class Air:
    c0 = 343
    rho = 1.204

def PouceEnMetre(Pouce):
    return Pouce * 0.0254

def QWL(Fr):
    return Air.c0/(4*Fr)


def BendCalculation(LengthTot , MaxWidth , MaxHeight , MaxDepth):
    Bend = 0
    i = 1
    while LengthTot/i > MaxDepth:
        i += 1
    Bend = i

    if Bend %2 != 0:
         Bend += 1
            
    if LengthTot < MaxDepth:
         LengthBend = LengthTot
         Bend = 1
               
        
    LengthBend = LengthTot/Bend
    WidthBend = MaxWidth/Bend
    HeightBend = MaxHeight/Bend
    
    return Bend , LengthBend , WidthBend , HeightBend


class RectangularTube:
    def __init__(self,Length,Height,Width,Thickness,DuctType="open"):
        self.Length = Length
        self.Height = Height
        self.Width  = Width
        self.Thickness = Thickness
        self.DuctType = DuctType
        self.Section_In  = self.Height * self.Width
        self.Section_Out = (self.Height + 2*self.Thickness)+(self.Width + 2*self.Thickness)
        self.Zc = Air.rho*Air.c0 / self.Section_In
        self.T = self.Transfermatrix(f)
        
        
    def Transfermatrix(self,f):
        K = 2*np.pi*f/Air.c0
        T =  np.array([[ np.cos(K*self.Length) , 1j*self.Zc*np.sin(K*self.Length) ] ,
                       [ 1j*np.sin(K*self.Length)/self.Zc , np.cos(K*self.Length) ]])
        return T
            
    # def Impedance(self,f):
    #     K = 2*np.pi*f/Air.c0
    #     Z = Air.rho*Air.c0/(np.pi*self.Section_Out**2)
    #     ZT = Z/1j*np.tan(K*self.Length)
    #     return ZT
            
    def Impedance(self,f):
        if self.DuctType == "Open":
            Imp = self.T[0,1,:] / self.T[1,1,:]
        if self.DuctType == "Closed":
            Imp = self.T[0,0,:] / self.T[1,0,:]
        return Imp

def BendConnection(T1,T2,WidthHole):
    L1p = T1.Length - WidthHole/2
    L2p = T2.Length - WidthHole/2
    LBp = T1.Width/2 + T2.Width/2
    return L1p , L2p , LBp

def Parallel(A,B):
    return (A*B)/(A+B)

class Connection:
    def __init__(self,f,*Tubes):
        self.T = self.connect(f,*Tubes)
    
    def connect(self,f,*Tubes):
        T_tot = np.zeros((2,2,len(f)))
        for k in range(len(f)):
            T = np.array([[1,0],[0,1]],dtype=complex)
            for i in range(len(Tubes)):
                T = T @ Tubes[i].Transfermatrix(f[k])
            T_tot[:,:,k] = np.abs(T)
        return T_tot
    
    def Impedance(self,f):
        return self.T[0,1,:] / self.T[1,1,:]
        
    
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

TUBES = []
ind = []

for i in range(NbrBend):
    TUBES.append(RectangularTube(LengthBend-WidthHole/2, WidthBend , HeightBend , 0.001 , DuctType="Open"))

if NbrBend > 1:
    for i in range(NbrBend):
        ind.append(i)
        if i %2 ==0:
            ind.remove(i)
        # print(ind)
    for i in ind:
        # print(i)
        _ , _ , LBp = BendConnection(TUBES[i-1], TUBES[i+1], WidthHole)
        TUBES.insert(i , RectangularTube(LBp, WidthHole , HeightBend , 0.001 , DuctType="Open"))
else:
    TUBES = TUBES

# Tube1 = RectangularTube(1/2, 0.1, 0.15, 0.01)
# Tube2 = RectangularTube(1.7/2, 0.1, 0.1, 0.01)

Test = Connection(f, *TUBES)
ImpTest = Test.Impedance(f)

# Tube = RectangularTube(1.7, 0.1, 0.1, 0.01)
# Imp = Tube.Impedance(f)

# %%
plt.close('all')

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


StraightOpenDuct = RectangularTube(1.7, 0.1, 0.1, 0.01 , DuctType="Open")
OpenDuct = RectangularTube(1.7*2/3, 0.1, 0.1, 0.01 , DuctType="Open")
ClosedDuct = RectangularTube(1.7/3, 0.1, 0.1, 0.01 , DuctType="Closed")
ZT = Parallel(OpenDuct.Impedance(f), ClosedDuct.Impedance(f))

plt.figure()
plt.semilogx(DATA[:,0],DATA[:,1] , label="AKABAK")
plt.semilogx(DATA[:,0],DATA[:,2] , label="AKABAK")

plt.semilogx(f,20*np.log10(np.abs(StraightOpenDuct.Impedance(f))) , label="Straight Open Duct")
# plt.semilogx(f,10*np.log10(np.abs(OpenDuct.Impedance(f))) , label="Open Duct")
# plt.semilogx(f,10*np.log10(np.abs(ClosedDuct.Impedance(f))) , label="Closed Duct")
plt.semilogx(f,20*np.log10(np.abs(ZT)) ,"k", label="Parallel Duct")

plt.legend()
plt.grid(which="both")
plt.show()

