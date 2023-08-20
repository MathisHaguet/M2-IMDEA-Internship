import numpy as np
import matplotlib.pyplot as plt
""" 
    This python code regroupe all the classes and functions used for
    the protyping of the labyrinth speaker. 
    Internship @ UdeS - Crash 
    Date   : 06/03/2023 to 28/07/2023
    Author : HAGUET Mathis 
    Supervisors : ROBIN Olivier and MELON Manuel
    
"""
#------------------------------------------------------------------------------
# Packages needed to run this code

from scipy.special import jv
import scipy.signal as signal
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider, Button
from pandas import read_excel
import csv
import pandas as pd
import matplotlib as mpl
# from cmcrameri import cm
# from tqdm import trange

#------------------------------------------------------------------------------
""" GENERAL CLASS & FUNCTION """
#------------------------------------------------------------------------------
class DSP:
    fs = int(48e3)
#------------------------------------------------------------------------------  
# Classe : fixed parameters in air
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
    
#------------------------------------------------------------------------------
# Function : Calculation of the loss inside a rectangular TUbe
def LossRecTube(f,Width,Height):
    Nit = int(len(f)/1000)
    w = 2 * np.pi * f # pulsation [rad.s-1]
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
    c_loss = np.sqrt(1/K/RHO)
    return RHO , c_loss

#------------------------------------------------------------------------------
# Function : inch to meter
def PouceEnMetre(Pouce):
    return Pouce * 0.0254
#------------------------------------------------------------------------------
# Function : Bessel function 
def J1(k,a,theta):
    if theta == 0:
        J1 = jv(0,k*a*np.sin(0)) + jv(2,k*a*np.sin(0))
    else :
        J1  = 2*jv(1,k*a*np.sin(theta))/(k*a*np.sin(theta))
    return J1
#------------------------------------------------------------------------------
""" TUBE CLASS & FUNCTION """
#------------------------------------------------------------------------------
# Function : calculate the length of a quater wave length tube for a given frequency
def QWL(Fr):
    return Air.c0/(4*Fr)
#------------------------------------------------------------------------------
# Function : calculate the first cut off frequency in a tube
def Fc(Length):
    return Air.c0/(2*Length)
#------------------------------------------------------------------------------
# Function : calculate the global dimenssions of a box depending on the number of bends
#            while keeping a constant volume
def Bendcalculation(L,W,H,N):
    Ltot = L/N
    Wtot = W
    Htot = H*N
    print("Bend Number : {:.2f}".format(N))
    print("Total Length : {:.2f} m".format(Ltot))
    print("Total Width : {:.2f} m ".format(Wtot))
    print("Total Hight : {:.2f} m".format(Htot))
    print("Volume init : {:.4f} m^3".format(L*W*H))
    print("Volume final : {:.4f} m^3".format(Ltot*Htot*Wtot))
    return N , Ltot , Wtot , Htot
#------------------------------------------------------------------------------

# bend optimisation pour des dimensions contraintes !
# def BendCalculation(LengthTot , MaxWidth , MaxHight , Bend):
    # opti contrainte longueure 
    # Bend = 0
    # i = 1
    # while LengthTot/i > MaxDepth:
    #     i += 1
    # Bend = i

    # if Bend %2 != 0: usefull to return in front with the hp
    #      Bend += 1
            
    # if LengthTot < MaxDepth:
    #      LengthBend = LengthTot
    #      Bend = 1
        
    # LengthBend = LengthTot/Bend
    # WidthBend = MaxWidth/Bend
    # HightBend = MaxHight/Bend
    
    # print("Bend Number :" , Bend)
    # print("Bend Length :" , LengthBend)
    # print("Bend Width :" , WidthBend)
    # print("Bend Hight :" , HightBend)
    
    # return Bend , LengthBend , WidthBend , HightBend
#------------------------------------------------------------------------------

# Function : calculate the effective length/resonance frequency of a plied tube 
def Leff(Length, Width , NbBend):
    Leff = Length - (4-np.pi)*(Width/2)*(NbBend-1)
    print("Effective Length : " ,Leff)
    Fres = Air.c0 / (4*(Length - (4-np.pi)*Width/2*(NbBend-1)))
    print("Freq res : " , Fres)
    return Leff , Fres

#------------------------------------------------------------------------------

# Class : create a rectangulare tube object
class RectangularTube:
    def __init__(self,f,Length,Height,Width,Thickness,DuctType="open"):
        self.Length = Length
        self.Height = Height
        self.Width  = Width
        self.f = f

      # tacking into acount losses into the model 
      
        # self.rho_loss , self.c_loss = LossRecTube(self.f, self.Width, self.Height)       
        # print(len(self.rho_loss))
        # print(len(self.c_loss))
        # print(len(self.f))
        
        
        self.rho_loss = Air.rho
        self.c_loss = Air.c0
        
        self.Thickness = Thickness
        self.Reff = self.Width/(np.sqrt(np.pi))
        self.Section_In = np.pi*self.Reff**2
        self.DuctType = DuctType
        self.Section_Out = (self.Height + 2*self.Thickness)+(self.Width + 2*self.Thickness)
        
        self.Zc = self.rho_loss*self.c_loss/ self.Section_In
        
        self.T   = self.Transfermatrix(f)
        self.Z2  = self.Radimp(f)
        self.Imp = self.Impedance(f)

        # if self.Section_In != None:
        #     self.Section_In = Section_In
        # else:
        #     self.Section_In  = self.Height * self.Width
        
    # Function : calculate the tansfermatrix of the tube
    def Transfermatrix(self,f):
        K = 2*np.pi*f/self.c_loss
        T =  np.array([[ np.cos(K*self.Length) , 1j*self.Zc*np.sin(K*self.Length) ] ,
                       [ 1j*np.sin(K*self.Length)/self.Zc , np.cos(K*self.Length) ]])
        return T
    
    # Function : calculate the radiation impedance of the tube
    def Radimp(self,f):
        w = 2*np.pi*f
        k = w/self.c_loss
        zc = self.rho_loss*self.c_loss / self.Section_In
        """ Simple Aproximation ( from Guide of the Guides - JPD )"""
        eta = 0.597
        # eta = 0.811
        R = 1 - (k*self.Reff*eta)**2 / 2
        # Z2 =  1j*np.tan(k*self.Length - 1j*1/2*np.log(R))
        Z2 = zc*(1-R)/(1+R)
        return Z2
    
    # Function : calculate the acoustic impedance of the tube
    def Impedance(self,f):
        Imp = 0
        if self.DuctType == "Open":
            Imp = self.T[0,1,:] / self.T[1,1,:]
        elif self.DuctType == "Closed":
            Imp = self.T[0,0,:] / self.T[1,0,:]
        return Imp
    
#------------------------------------------------------------------------------
# Class : create a cylindrical tube object
class CylindricalTube:
    def __init__(self,f,Length,Radius,Thickness,DuctType="open"):
        self.Length = Length
        self.Radius  = Radius
        self.Thickness = Thickness
        self.DuctType = DuctType
        self.f = f
        self.Section_In  = np.pi*self.Radius**2
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
        w = 2*np.pi*f
        k = w/Air.c0
        
        """Unflanged case"""
        Beta = 1/2
        eta = 0.6133
        n1 = 0.167
        d1 = 1.393
        d2 = 0.457
        
        R = -(1-n1*1j*k*self.Radius) / (1-d1*1j*k*self.Radius - d2 * k**2 * self.Radius**2) 
        Z2 = self.Zc*(1+R)/(1-R) 
        # a1 = 0.800
        # a2 = 0.266
        # a3 = 0.0263
        # b1 = 0.0599
        # b2 = 0.238
        # b3 = -0.0153
        # b4 = 0.00150
        
        """FLanged case"""
        # Beta = 1
        # eta = 0.8216
        # a1 = 0.730
        # a2 = 0.372
        # a3 = 0.0231
        # b1 = 0.244
        # b2 = 0.723
        # b3 = -0.0198
        # b4 = 0.00366
        
        # R = (1 + a1*(k*self.Radius)**2) / \
        #     (1 + (Beta+a1)*(k*self.Radius)**2)
        # + a2*(k*self.Radius)**4 + a3*(k*self.Radius)**6
        
        
        # R = -(1-n1*1j*k*self.Radius) / (1-d1*1j*k*self.Radius - d2 * k**2 * self.Radius**2) 
        
        # L = self.Radius * eta * (1 + b1*(k*self.Radius)**2) / \
        #     (1 + b2*(k*self.Radius)**2)
        # + b3*(k*self.Radius)**4 + b4*(k*self.Radius)**6 
            
        # Z2 = self.Zc * 1j*np.tan(k*L - 1j*0.5*np.log(R)) 

        """ Simple Aproximation ( from Guide of 2 Guides - JPD )"""
        # eta = 0.6133
        # # # eta = 0.811
        # R = 1 - (k * self.Radius * eta )**2 / 2
        # Z2 = self.Zc * (1-R)/(1+R) 
        
        return Z2
    
    def Impedance(self,f):
        Imp = 0
        if self.DuctType == "Open":
            Imp = self.T[0,1,:] / self.T[1,1,:]
        if self.DuctType == "Closed":
            Imp = self.T[0,0,:] / self.T[1,0,:]
        return Imp
#------------------------------------------------------------------------------
# Function : return a vector containing the width values of a voigt pipe 
def VoigtSection(H,D,W,L,NbrPoints):
    y = np.linspace(0, L , NbrPoints)
    Section = (-D*H / L) * y + (W*H)
    Width = (D/L)*y + W
    return Width
#------------------------------------------------------------------------------
# Function : calculate the equivalent length of 2 tubes connected together 
#            based on the J.Panzer paper "Acoustics - Ducts with Bend".
def BendConnection(T1,T2,WidthHole):
    L1p = T1.Length - WidthHole/2
    L2p = T2.Length - WidthHole/2
    LBp = T1.Width/2 + T2.Width/2
    return L1p , L2p , LBp

#------------------------------------------------------------------------------
#Function : put 2 object in parallel from electric point vue 
def Parallel(A,B):
    return (A*B)/(A+B)

#------------------------------------------------------------------------------
# Function : Calculate frequency of the mode of a room 
def RoomMode(Lx,Ly,Lz,l,m,n):
    F = 343/2 * np.sqrt((l/Lx)**2 + (m/Ly)**2 + (n/Lz)**2)
    print(" f = {:.1f} Hz".format(F))
    return F

#------------------------------------------------------------------------------
# Function : create a monopole object
def Monopole(xs, ys, x, y, w, k, q , tau):
    rs = np.sqrt((xs - x)**2 + (ys - y)**2)
    p_mono = 1j * Air.rho * k * Air.c0 * q * np.exp(-1j * k * rs) / (2 * np.pi * rs) * np.exp(-1j*w*tau)
    return p_mono
#------------------------------------------------------------------------------
#Class : connect tubes together using their transfer matrix
class Connection:
    def __init__(self,f,*Tubes):
        self.Length  = Tubes[-1].Length
        self.Zc = Tubes[-1].Zc
        self.LengthTot  = self.Tube_Info(*Tubes)
        self.T = self.connect(f,*Tubes)
        self.Z2 = self.Radimp(f , *Tubes)
        self.Imp = self.Impedance(f , *Tubes)
        self.ConnectImp = self.T[0,1,:] / self.T[1,1,:]
        self.Tubes = Tubes
        self.DuctType="Connection"
        
    def Tube_Info(self,*Tubes):
        TotLenght = 0
        for Tube in Tubes:
            TotLenght += Tube.Length
        return TotLenght
        
    def connect(self,f,*Tubes):
        T_tot = np.zeros((2,2,len(f)),dtype=complex)
        for k in range(len(f)):
            T = np.array([[1,0],[0,1]],dtype=complex)
            for i in range(len(Tubes)):
                T = T @ Tubes[i].Transfermatrix(f[k])
            T_tot[:,:,k] = np.abs(T)
        return T_tot

    def Radimp(self,f,*Tubes):
        if type(Tubes[-1]) == CylindricalTube :
            w = 2*np.pi*f
            k = w/Air.c0
            
            """Unflanged case"""
            Beta = 1/2
            eta = 0.6133
            n1 = 0.167
            d1 = 1.393
            d2 = 0.457
            R = -(1-n1*1j*k*Tubes[-1].Radius) / (1-d1*1j*k*Tubes[-1].Radius - d2 * k**2 * Tubes[-1].Radius**2) 
            Z2 = self.Zc*(1+R)/(1-R) 
            
        else:
            w = 2*np.pi*f
            k = w/Air.c0
            
            """ Simple Aproximation ( from Guide of 2 Guides - JPD )"""
            eta = 0.597
            # eta = 0.811
            R = 1 - (k*Tubes[-1].Reff*eta)**2 / 2
            # Z2 =  -1j*np.tan(k*L - 1j*1/2*np.log(R)))
            Z2 = self.Zc*(1-R)/(1+R)
            
        return Z2
    

    # For now this Impedance Function is working only ofr the case of 2 ducts
   
    def Impedance(self,f,*Tubes):
        Imp = np.zeros((len(Tubes),len(f)), dtype='complex')
        for i in range(len(Tubes)):
            if Tubes[i].DuctType == "Open":
                Imp[i,:] = Tubes[i].T[0,1,:] / Tubes[i].T[1,1,:]

            if Tubes[i].DuctType == "Closed":
                Imp[i,:] = Tubes[i].T[0,0,:] / Tubes[i].T[1,0,:]
        
        Imp = Parallel(Imp[1,:], Imp[0,:])  ## To do general case with Parallel
        return Imp

#------------------------------------------------------------------------------
#Class : create a loudspeaker object
class Loudspeaker():
    def __init__(self ,f, Fr , Re , Le , Bl , Mms , Qms , rd  , r , BackVolume):
        self.Fr = Fr
        self.Wr = self.Fr * 2 * np.pi
        self.Re = Re
        self.Le = Le
        self.Bl = Bl
        self.Mms = Mms
        self.Qms = Qms
        self.rd = rd
        self.Sd = rd**2 * np.pi
        self.Zc = Air.rho*Air.c0/self.Sd
        self.Cms = 1 / (self.Wr**2 * self.Mms)
        self.Kms = 1 / self.Cms
        self.Rms = 1 / (self.Wr * self.Cms * self.Qms)
        self.f = f
        self.Qes = self.Re / (self.Wr*self.Cms*(self.Bl)**2)
        self.Qts = (self.Qms*self.Qes)/(self.Qms+self.Qes)
        self.Vas = Air.rho * Air.c0**2 * self.Cms * self.Sd**2
        self.BackVolume = BackVolume/1000

        if self.BackVolume == False:
            self.alpha = 0
            self.Vb = 0
        else:
            self.Vb = self.BackVolume
            self.alpha = self.Vas / self.Vb
            
        self.Fc = np.sqrt(1+self.alpha) * self.Fr
        self.Wc = self.Fc * 2 * np.pi
        self.Qec = np.sqrt(1+self.alpha) * self.Qes
        self.Qmc = np.sqrt(1+self.alpha) * self.Qms
        self.Qtc = Parallel(self.Qmc, self.Qec)        
        
        self.Cab = self.BackVolume / (Air.rho * Air.c0**2)
        
        self.Ze = self.Ze(f)
        self.Zms = self.Zms(f)
        self.r = r
        self.Zrad = self.Zrad(f)

        self.x = 0  # loudspeaker x position (baffle x-position) [m]
        self.y = r  # loudspeaker y position in the baffle [m]
            
        # calculate the digital filter coefficients
        self.b, self.a = self.digital_filter_coefficients()
    

    def Ze(self,f):
        s = 1j*2*np.pi*self.f
        return self.Re + s * self.Le + ( ( s * self.Cms*self.Bl**2 ) / ( s**2 * self.Mms * self.Cms + s * self.Rms * self.Cms + 1))

    def Zms(self,f):
        s = 1j*2*np.pi*self.f
        return self.Rms + s*self.Mms + 1/(s*self.Cms)
        
    def Zrad(self,f):
        k = 2*np.pi*self.f / Air.c0
        s = 1j*2*np.pi*self.f
        q = ((1*self.Sd)/(self.Bl*self.Qec)) * ((s/self.Wc) / ( (s)**2/(self.Wc**2) + (s/(self.Qtc*self.Wc)) +1 ))

        Pinfb_piston = 1j*k*Air.rho*Air.c0 *q* np.exp(-1j*k*self.r)/(2*np.pi*self.r) * J1(k,self.rd*2,0)
        Pdrive_piston = (1+np.cos(0))/2 * J1(k,self.rd*2,0) * Pinfb_piston
        Pe_piston = 2*jv(1,k*self.rd*2)/(k*self.rd*2)*(-np.cos(0)/2) * np.exp( -1j * k * self.rd) * jv(0,k*self.rd*np.sin(0)) * Pinfb_piston
        Pt_piston = (Pdrive_piston + Pe_piston ) / Pinfb_piston
        
        Z = self.Zc*k*self.Sd*np.exp(-1j*k*self.r) / 2*np.pi*self.r
        return Z
    
    def analog_filter_coefficients(self):
        """ analog filter simulatig transfer function of a loudspeaker
            acceleration vs. voltage 
            H(s) = A(s)/U(s)
        Returns:
            B, A: coefficients of analog filter
        """
        
        A = np.zeros(4)
        A[0] = self.Le*self.Mms/self.Bl
        A[1] = (self.Re*self.Mms+self.Le*self.Rms)/self.Bl
        A[2] = (self.Re*self.Rms+self.Le*self.Kms)/self.Bl+self.Bl
        A[3] = self.Re*self.Kms/self.Bl

        B = np.array([1, 0, 0])
        return B, A

    def digital_filter_coefficients(self):
        """ digital filter simulatig transfer function of a loudspeaker
            acceleration vs. voltage 
            H(z) = A(z)/U(z)
        Returns:
            b, a: coefficients of digital filter        
        """
        B, A = self.analog_filter_coefficients()
        b, a = signal.bilinear(B, A, DSP.fs)
        return b, a
    
    def ThieleAndSmallParam(self):
        print("###########################")
        print("Cms = {:.0f} um/N ".format(self.Cms*1e6))
        print("Mms = {:.2f} g ".format(self.Mms*1e3))
        print("Fr  = {:.2f} Hz".format(self.Fr))
        print("Qms = {:.2f}".format(self.Qms))
        print("Qes = {:.3f} ".format(self.Qes))
        print("Qts = {:.2f}".format(self.Qts))
        print("Rms = {:.2f} Kg/s".format(self.Rms))
        print("Bl  = {:.2f} N/A ".format(self.Bl))
        print("Vas = {:.2f} L".format(self.Vas*1e3))
        print("Le  = {:.2f} mH".format(self.Le*1e3))
        print("Re  = {} ohm".format(self.Re))
        print("###########################")

#------------------------------------------------------------------------------
# Class : create a linstening point object 
class ListeningPoint:
    """ class for storing the coordinates r, theta of a point in 2D Space
        and creating automatically coordinates x and y
    """

    def __init__(self, r, theta):
        self.r = r  # distance from (x=0, y=0) [m]
        self.theta = theta  # angle [deg]

        self.x = r * np.cos(np.pi/180 * theta)
        self.y = r * np.sin(np.pi/180 * theta)

#------------------------------------------------------------------------------
class Propagation:
    """ class simulation a propagation of an acoustic wave between two points 
        (listening point and speaker) using an FIR filter
    """

    def __init__(self, speaker, point):
        self.point = point  # object of class ListeningPoint
        self.speaker = speaker  # object of class Loudspeaker

        self.a = 1  # filter coefficients a
        self.b = self.filter_coefficients()  # filter coefficients b

    def get_distance(self):
        """ calculate the distance r between listening point and speaker
        """
        return np.sqrt((self.speaker.x-self.point.x)**2 + (self.speaker.y-self.point.y)**2)

    def filter_coefficients(self):
        """ FIR filter design simulating a propagation of an acoustic wave between two points """
        
        r = self.get_distance()  # distance between listening point and speaker
        N = int(DSP.fs*r/Air.c0)  # [samples]
        b = np.zeros(N+1)
        b[N] = Air.rho*self.speaker.Sd / (2*np.pi*r)
        return b
#------------------------------------------------------------------------------
# Class : calculate the acoustic radiation of a loudspeaker into a tube into air
class Radiation():
    def __init__(self,f,LS,Duct,U,r):
        
        self.f = f
        self.w = 2*np.pi*f
        self.k = self.w/Air.c0
        self.LS = LS
        self.r = r
        self.Duct = Duct
        self.U = U
        s = 1j*self.w
        
        if self.Duct.DuctType == "Connection":
            Length = self.Duct.LengthTot
            Imp = self.Duct.ConnectImp
            print("Ici")
        else:
            Imp = self.Duct.Imp
            Length = self.Duct.Length
        # print(Length)
        
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
        self.ZTL = Imp                                                         # Duct Imp -> Acoustic Domain
        self.Zab = 1/(s*self.LS.Cab)                                           # Imp at Back -> Acoustic Domain
        self.Zaf = self.ZRT                                                    # Imp at the Front -> Acoustic Domain
        
        self.Q0 = self.P / (self.Zae  + self.Zas + self.Zab + self.ZTL + self.Zaf)
        self.Qv = self.Q0 / (np.cos(self.k*Length) + 1j*np.sin(self.k*Length)/self.Duct.Zc) 
                
        self.Prear = 1j*self.k*Air.rho*Air.c0 * self.Qd * np.exp(-1j*self.k*(self.r)) / (4*np.pi*(self.r))
        self.Pvent = 1j*self.k*Air.rho*Air.c0 * self.Qv * np.exp(-1j*self.k*(self.r)) / (4*np.pi*(self.r))
        
        # self.Ptot =  1j*self.k*Air.rho*Air.c0 * (self.Qd * np.exp(-1j*self.k*(self.r))/2*np.pi*(self.r) - self.Qv * np.exp(-1j*self.k*self.r)/2*np.pi*self.r )

#------------------------------------------------------------------------------

def FRF(N, *filters):
    """ claculate FRF as a multiplication of individual FRFs provided in filters list

    Args:
        N: number of points to plot the FRF

    Returns:
        F, H: frequency axis (in [Hz]), FRF of the filters in series 
    """
    H = 1
    for filter in filters:
        W, Htemp = signal.freqz(filter.b, filter.a, worN=N)
        H *= Htemp

    F = W/(2*np.pi) * DSP.fs

    return F, H


""" Data post traitment and plot figures generation"""
#------------------------------------------------------------------------------

def Plot_free_field(x, y, xS, yS, p, text_xlabel, text_ylabel, text_title):
    fig, ax = plt.subplots(figsize=(10,10))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    p_dBSPL = 20 * np.log10(np.abs(p) / 2e-5 / np.sqrt(2))
    amp_plot = ax.pcolormesh(x, y, p_dBSPL, cmap=cm.batlow,vmin=0,vmax=140)
    cbar = fig.colorbar(amp_plot, cax=cax, orientation='vertical')
    cbar.set_label('Module Pression (dB SPL)', rotation=270, labelpad=20)
    ax.set_title(text_title, y=1.04)
    ax.set_xlabel(text_xlabel)
    ax.set_ylabel(text_ylabel)
    ax.set_aspect('equal', 'box')
#------------------------------------------------------------------------------

def MeshGrid(Lx,Ly,d0):
    Nx = int(Lx / d0)  # nombre de points du maillage selon Ox
    Ny = int(Ly / d0)  # nombre de points du maillage selon Oy
    x_vec = np.linspace(0,Lx, Nx)
    y_vec = np.linspace(-Ly/2,Ly/2, Ny)
    [x, y] = np.meshgrid(x_vec, y_vec)
    return x,y

#------------------------------------------------------------------------------
# Function : replace or change a caratere ina  vector
def ReplaceFloat(Vec,C1,C2):
    for i in range(len(Vec)):
        temp = Vec[i].replace(C1,C2)
        Vec[i] = float(temp)
    return  Vec  
     
#------------------------------------------------------------------------------
# Function : create a dictionary of the FRF data
#           usefull to formalize the data and compare simulation/measure
def Dict(f,FRF,label=None,color=None,linestyle="-",PowerSpec=None):
    if PowerSpec is True:
        P = 2*np.angle(FRF)
    else:
        P = np.angle(FRF)
    return {"f_axis":f,"MFRF":np.abs(FRF),"PFRF":P,"label":label, "color":color,"linestyle":linestyle}

#------------------------------------------------------------------------------
# Function : load measurement data from a CSV file and create a dictionnary of them
def LoadData(path , label=None ,color=None,linestyle="-" , PowerSpec=None):
    data = pd.read_csv(path, sep=';',skiprows=19).to_numpy() #,usecols=[1,2,3])) 
    f = ReplaceFloat(data[:,0],',','.').astype(float)
    Real = ReplaceFloat(data[:,1],',','.').astype(complex)
    Imag = ReplaceFloat(data[:,2],',','.').astype(complex)
    FRF = Real + 1j*Imag
    if PowerSpec is True:
        FRF = np.sqrt(FRF)
    else:
        FRF = FRF
    return Dict(f, FRF, label,color,linestyle,PowerSpec=PowerSpec)

#------------------------------------------------------------------------------
# Function : Display the FRF
def PlotFRF(*FRF , Title="" ,SAR=False , linewidth=3 , loc="best" , Ylabel=r"$P_{vent}$"):
    fig , ax = plt.subplots(nrows=2 , sharex=True , facecolor='white')
    
    MAX_dB = 0
    MIN_dB = 60
    
    for H in FRF:
        if MAX_dB <= max(20 * np.log10(H['MFRF']/2e-5)) + 5:
            MAX_dB =  max(20 * np.log10(H['MFRF']/2e-5)) + 5
        else:
            MAX_dB = MAX_dB
        if MIN_dB >= min(20 * np.log10(H['MFRF']/2e-5)) - 5: 
            MIN_dB = min(20 * np.log10(H['MFRF']/2e-5)) - 5
        else:
            MIN_dB = MIN_dB
        
        # print(MIN_dB)
        ax[0].semilogx( H['f_axis'], 20 * np.log10(H['MFRF']/2e-5), label=H["label"] , linewidth=linewidth , color=H["color"],linestyle=H["linestyle"])
        ax[1].semilogx(H['f_axis'], H["PFRF"]  , label=H["label"] , linewidth=linewidth , color=H["color"],linestyle=H["linestyle"])     
        
    ax[0].set_ylim([MIN_dB , MAX_dB])
    ax[1].set_ylim([-np.pi-0.6,np.pi+0.6])
    
    if SAR is True:
        ax[0].fill_between(H['f_axis'],MAX_dB,MIN_dB,where=H['f_axis']<=57, facecolor='grey', alpha=.3 , label="SAR")
        ax[1].fill_between(H['f_axis'], np.pi+0.6,-np.pi-0.6,where=H['f_axis']<=57, facecolor='grey', alpha=.3, label="SAR")

    ax[0].set_xlim([H['f_axis'][0],H['f_axis'][-1]])
    ax[0].legend(loc=loc,fontsize=15)
    ax[0].set_ylabel('|{}| [dBSPL]'.format(Ylabel) ,fontsize=20)

    ax[1].set_ylabel(r'$\mathbb{\phi}$'+'({}) [rad]'.format(Ylabel) ,fontsize=20)
    ax[1].set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r"$-\pi$" ,r"$-\frac{\pi}{2}$", r"$0$" ,r"$\frac{\pi}{2}$",r"$\pi$" ],fontsize=15)
    ax[1].set_xlabel('Frequency [Hz]',fontsize=20)

    for i in range(2):
        ax[i].grid(which="both")

    ax[0].set_title("{}".format(Title))
    
    plt.tight_layout()
    plt.show()

#------------------------------------------------------------------------------
# Function : return the indice of freq value in a f axis vector
def indfreq(f,ind):
    i = np.argwhere(f==f.flat[np.abs(f - ind).argmin()])[0][0]
    return i
    
#------------------------------------------------------------------------------

def Fplot(*SIG ,  Title=""  , linewidth=3 , loc="best"):
    fig , ax = plt.subplots(nrows=2 , sharex=True , facecolor='white')

    MAX_dB = 0
    
    for H in SIG:    
        ax[0].semilogx( H['f_axis'], 20 * np.log10(H['MFRF']/2e-5) , label=H["label"] , linewidth=linewidth)
        ax[1].semilogx(H['f_axis'],H["PFRF"]  , label=H["label"] , linewidth=3)     
        
    ax[0].legend(loc=loc,fontsize=15)
    ax[0].set_ylabel(r'$|P_l|$ [dBSPL]' ,fontsize=20)

    ax[1].set_ylabel(r'$\mathbb{\phi}(P_{vent})$ [rad]' ,fontsize=20)
    ax[1].set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r"$-\pi$" ,r"$-\frac{\pi}{2}$", r"$0$" ,r"$\frac{\pi}{2}$",r"$\pi$" ])
    
    for i in range(2):
        ax[i].grid(which="both")
        ax[i].set_xlabel('Frequency [Hz]',fontsize=20)
    
    ax[0].set_title("{}".format(Title))
    
    plt.tight_layout()
    plt.show()
 
#------------------------------------------------------------------------------
# Function : return a color fade between to chosen color
def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)
  
#------------------------------------------------------------------------------
# Function : display the directivity 
def PlotDirect(folder, deg_val , freq_val , dB=None):
    D = np.zeros((len(deg_val),len(freq_val)))
    THETA = np.zeros((len(deg_val),len(freq_val)))
    fig , ax = plt.subplots(subplot_kw={"projection":"polar"}, facecolor='white')
    ax.set_thetamax(180)
    
    for f in range(len(freq_val)):
        for i in range(len(deg_val)):
            path = folder.format(deg_val[i])
            data_temp = LoadData(path , "°{}".format(deg_val[i]) , PowerSpec=True)
            indf = indfreq(data_temp["f_axis"], freq_val[f])
            D[i,f]=data_temp["MFRF"][indf]
            THETA[i,f] = deg_val[i]*np.pi/180
            
        # LINEAR
        if dB is None:
            ax.plot(THETA[:,f], D[:,f]/max(D[:,f])  , 
                    color=colorFader("orange", "blue", f/len(freq_val)),
                    linewidth=2 , label = "{:.0f} Hz".format(freq_val[f]))
            # ax.set_title("Directivity patern (linear normalized)")

        else:   
            ax.plot(THETA[:,f],20*np.log10(D[:,f]/2e-5)  , 
                color=colorFader("orange", "blue", f/len(freq_val)),
                linewidth=2 , label = "{:.0f} Hz".format(freq_val[f]))
            # ax.set_title("Directivity patern [dB SPL]")

    
    plt.legend(loc="center")
    plt.tight_layout()
    plt.show()