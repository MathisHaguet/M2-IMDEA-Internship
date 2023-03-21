import numpy as np

class Air:
    c0 = 343
    rho = 1.204

def PouceEnMetre(Pouce):
    return Pouce * 0.0254

def QWL(Fr):
    return Air.c0/(4*Fr)

def Fc(Length):
    return Air.c0/(2*Length)


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
    def __init__(self,f,Length,Height,Width,Thickness,DuctType="open"):
        self.Length = Length
        self.Height = Height
        self.Width  = Width
        self.Thickness = Thickness
        self.Section_In = self.Width * self.Height
        self.DuctType = DuctType
        self.Section_Out = (self.Height + 2*self.Thickness)+(self.Width + 2*self.Thickness)
        self.f = f
        self.Zc = Air.rho*Air.c0 / self.Section_In
        
        self.T   = self.Transfermatrix(f)
        self.Imp = self.Impedance(f)

        # if self.Section_In != None:
        #     self.Section_In = Section_In
        # else:
        #     self.Section_In  = self.Height * self.Width
        
        
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
        Imp = 0
        if self.DuctType == "Open":
            Imp = self.T[0,1,:] / self.T[1,1,:]
        if self.DuctType == "Closed":
            Imp = self.T[0,0,:] / self.T[1,0,:]
        
        return Imp
    
class CylindricalTube:
    def __init__(self,f,Length,Radius,Thickness,DuctType="open"):
        self.Length = Length
        self.Radius  = Radius
        self.Thickness = Thickness
        self.DuctType = DuctType
        self.f = f
        self.Section_In  = np.pi*Radius**2
        self.Zc = Air.rho*Air.c0 / self.Section_In
        self.T   = self.Transfermatrix(f)
        self.Z2 = self.Radimp(f)
        self.Imp = self.Impedance(f)
        
        # self.Section_Out = (self.Height + 2*self.Thickness)+(self.Width + 2*self.Thickness)
        
        # if self.Section_In != None:
        #     self.Section_In = Section_In
        # else:
        #     self.Section_In  = np.pi*Radius**2
        
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

    def Radimp(self,f):
        Z2 = 0
        w = 2*np.pi*f
        k = w/Air.c0
        n1 = 0.167
        d1 = 1.393
        d2 = 0.457
        # n1 = 0.182
        # d1 = 1.825
        # d2 = 0.649
        # Beta = 1/2
        # eta = 0.6133
        # a1 = 0.800
        # a2 = 0.266
        # a3 = 0.0263
        # b1 = 0.0599
        # b2 = 0.238
        # b3 = 0.0153
        # b4 = 0.00150
        
        Beta = 1
        eta = 0.8216
        a1 = 0.730
        a2 = 0.372
        a3 = 0.0231
        b1 = 0.244
        b2 = 0.723
        b3 = 0.0198
        b4 = 0.00366
        
        # Z2 = self.Zc * (k*self.Radius)**2/2
        
        # Z2 = ((n1-d1)*1j*k*self.Radius + d2*(1j*k*self.Radius)**2)/ \
        #       (2 - (d1+n1)*1j*k*self.Radius + d2*(1j*k*self.Radius)**2)

        R = (1 + a1*(k*self.Radius)**2) / \
            (1 + (Beta+a1)*(k*self.Radius)**2 + a2*(k*self.Radius)**4 + a3*(k*self.Radius)**6)
        
        L = (self.Radius*eta * (1 + b1*(k*self.Radius)**2)) / \
            (1 + b2*(k*self.Radius)**2 + b3*(k*self.Radius)**4 + b4*(k*self.Radius)**6)
           
        Z2 = -1j*np.tan(k*L - 1j*1/2*np.log(R))
        
        return Z2


    def Impedance(self,f):
        Imp = 0
        if self.DuctType == "Open":
            Imp = (self.T[0,0,:]*self.Z2 + self.T[0,1,:] )/ (self.T[1,1,:] + self.T[1,0,:]*self.Z2)

        if self.DuctType == "Closed":
            Imp = self.T[0,0,:] / self.T[1,0,:]
        return Imp
    
def VoightSection(H,D,W,L,NbrPoints):
    y = np.linspace(0, L , NbrPoints)
    S = (-D*H / L) * y + (W*H)
    return S 

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

    def RadImp(self,f):
        return self.T[0,1,:]
    