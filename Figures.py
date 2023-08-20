""" 
    This code include most of the post traitment of the simulated/measured data 
    and the figures used in the repport. 
    Internship @ UdeS - Crash 
    Date   : 06/03/2023 to 28/07/2023
    Author : HAGUET Mathis 
    Supervisors : ROBIN Olivier and MELON Manuel
    
"""
#------------------------------------------------------------------------------
# Packages needed to run this code

import numpy as np
from pandas import read_excel
import csv
import matplotlib.pyplot as plt
import os
from ClassFunction import *
from scipy.io import wavfile
import pandas as pd
import matplotlib.style
import matplotlib as mpl

mpl.style.use('classic')


#%%----------------------------------------------------------------------------

f  = np.logspace(np.log10(20),np.log10(20e3),200)
point = ListeningPoint(r=1, theta=0)

#%%----------------------------------------------------------------------------

path = "DATA/MEASUREMENTS/MESURES SSA/MESURES MATHIS/SS_Hp_Alone.csv"
HpAlone = LoadData(path , "Measure" , color=None)

# WOOFER USED FOR TESTS PROTOTYPES
TB_W4_1658SB = Loudspeaker(HpAlone["f_axis"] ,  61.1 , 4.2 , 0.21e-3 , 7.08 , 12.97e-3, 10.4, 5e-2 , 1 , BackVolume=8)

# WOOFER FINAL PROTOTYPE
P_SDS_135F25CP02_04 = Loudspeaker(f,63 , 3.3, 1.22e-3, 7.39, 20.6e-3, 3.89, PouceEnMetre(5.5)/2, 1, BackVolume=False)
PropaP = Propagation(P_SDS_135F25CP02_04, point)
_, HP = FRF(2*np.pi*f/DSP.fs, P_SDS_135F25CP02_04 , PropaP)
RAD_P = Dict(f, HP, "Woofer")

# TWEETER FINAL PROTOTYPE
P_DX25TG59_04 = Loudspeaker(f, 640, 2.93, 0.02e-3, 2.41 , 0.40e-3, 2.46, 0.0323/2 , 1, BackVolume=False)
PropaP_T = Propagation(P_DX25TG59_04, point)
_, HP_T = FRF(2*np.pi*f/DSP.fs, P_DX25TG59_04 , PropaP_T)
RAD_P_T = Dict(f, HP_T, "Tweeter")

#%%----------------------------------------------------------------------------
""" HP ALONE """
Fplot(RAD_P , RAD_P_T , Title="")

#%%----------------------------------------------------------------------------
""" HP ALONE + REAR VOLUME """

path = "DATA/MEASUREMENTS/MESURES SSA/MESURES MATHIS/SS_Hp_Alone.csv"
HpAlone = LoadData(path , "Measure" , color=None)

CYL = CylindricalTube(HpAlone["f_axis"], 1.4, 0.1/2, 0.01, DuctType="Open")
CYL_RAD = Radiation(HpAlone["f_axis"], TB_W4_1658SB , CYL , 1 , 1)
HpRad_AlonePython = Dict(CYL_RAD.f, CYL_RAD.Prear , "Model" , color="k")

PlotFRF(HpAlone,HpRad_AlonePython,SAR=True,Ylabel=r"$P_{rad}$")

#%%----------------------------------------------------------------------------
""" HP END + CYL TUBE  """

path = "DATA/MEASUREMENTS/MESURES SSA/MESURES MATHIS/SS_CYL_Octaves.csv"
HpCyl = LoadData(path , "Measure")

path = os.getcwd()
DATA = np.loadtxt("DATA/AKABAK SIMU/CYL_LEM_vent@1mBem.txt") 
CYL_AKABAK = Dict(DATA[:,0], DATA[:,1] + 1j*DATA[:,2] , "AKABAK",color="r")

TB_W4_1658SB = Loudspeaker(HpCyl["f_axis"] ,  61.1 , 4.2 , 0.21e-3 , 7.08 , 12.97e-3, 10.4, 5e-2 , 1 , BackVolume=8)

# CYL = CylindricalTube(HpCyl["f_axis"], 1.43, 0.1/2, 0.01, DuctType="Open")
CYL = RectangularTube(HpCyl["f_axis"], 1.43, 0.1, 0.1, 0.02 , DuctType="Open")

CYL_RAD = Radiation(HpCyl["f_axis"], TB_W4_1658SB , CYL , 1 , 1)
CYL_PYTHON = Dict(CYL_RAD.f, CYL_RAD.Pvent , "Analytic model" , color="k")

# Cylindrical Straight Tube
PlotFRF(HpCyl,CYL_PYTHON,CYL_AKABAK,Title="",linewidth=3,SAR=True)

#%%----------------------------------------------------------------------------
""" HP L/3 + REC Tapered TUBE  """



path = os.getcwd()
DATA = np.loadtxt(path + "/DATA/AKABAK SIMU/COMPA SLOPE VOIGT/REC_Tapered1Quart_Hp1Third.txt") 
REC_TAPERED_AKABAK = Dict(DATA[:,0], DATA[:,1] + 1j*DATA[:,2] , "AKABAK" , color = "red")

path = "DATA/MEASUREMENTS/MESURES SSA/MESURES MATHIS/SS_REC_Tapered_Octaves.csv"
HpRecTapered = LoadData(path , "Measure")

TB_W4_1658SB = Loudspeaker(HpRecTapered["f_axis"] ,  61.1 , 4.2 , 0.21e-3 , 7.08 , 12.97e-3, 10.4, 5e-2 , 1 , BackVolume=8)

# Rectangular Tapered Tube & Speaker in position 2
PlotFRF(REC_TAPERED_AKABAK, HpRecTapered  ,SAR=True , loc="lower left")

#%%----------------------------------------------------------------------------

""" HP END + REC Tapered TUBE """
path = "DATA/AKABAK SIMU/COMPA SLOPE VOIGT/RecSlope_HpEnd.txt"
DATA = np.loadtxt(path)
DICT_TEMP = []
Labels = ["H/4","","H/2","","2H/3","","3H/4","","H"]
c1='blue' 
c2='orange' 
n = 9

for i in range(1,10):
    if i%2 != 0:
        real = i
        imag = real+1
        # print(real-1, imag)
    DICT_TEMP.append(Dict(DATA[:,0] , DATA[:,real] + 1j*DATA[:,imag] , label=Labels[real-1],color=colorFader(c1, c2, real/n)))

PlotFRF(DICT_TEMP[0],DICT_TEMP[2],DICT_TEMP[4],DICT_TEMP[6],DICT_TEMP[8], Title="" , linewidth=3)

#%%----------------------------------------------------------------------------

""" HP END + LABYRINTHE TUBE  """

path = "DATA/MEASUREMENTS/MESURES SSA/MESURES MATHIS/SS_LABYRINTHE.csv"
HpLAB = LoadData(path , "Measure")

path = os.getcwd()
DATA = np.loadtxt(path + "/DATA/AKABAK SIMU/FOLDED_N5.txt") 
FOLDED_AKABAK_REC = Dict(DATA[:,0], DATA[:,1] + 1j*DATA[:,2] , "Rectangular" ,color="r")

TB_W4_1658SB = Loudspeaker(HpLAB["f_axis"] ,  61.1 , 4.2 , 0.21e-3 , 7.08 , 12.97e-3, 10.4, 5e-2 , 1 , BackVolume=8)
CYL = CylindricalTube(HpLAB["f_axis"], 1.43, 0.1/2, 0.01, DuctType="Open")
CYL_RAD = Radiation(HpLAB["f_axis"], TB_W4_1658SB , CYL , 1 , 1)
CYL_PYTHON = Dict(CYL_RAD.f, CYL_RAD.Pvent , "PYTHON")

PlotFRF(HpLAB,FOLDED_AKABAK_REC , Title = "" , linewidth=3 , SAR=True)

#%%----------------------------------------------------------------------------

""" FOLDED 5 BENDS CYL  """

path = os.getcwd()
DATA = np.loadtxt(path + "/DATA/AKABAK SIMU/FOLDED_5Bends_Cyl.txt")
DATA_SIDE =  np.loadtxt(path + "/DATA/AKABAK SIMU/FOLDED_5Bends_Side_Cyl.txt")

FOLDED_AKABAK_CYL = Dict(DATA[:,0], DATA[:,1] + 1j*DATA[:,2] , "Cylindrical" , color="k")
FOLDED_SIDE_AKABAK = Dict(DATA_SIDE[:,0], DATA_SIDE[:,1] + 1j*DATA_SIDE[:,2] , "Folded Side : Leff = 1.04832 m")
TB_W4_1658SB = Loudspeaker(f,  61.1 , 4.2 , 0.21e-3 , 7.08 , 12.97e-3, 10.4, 5e-2 , 1 , BackVolume=8)


PlotFRF(FOLDED_AKABAK_CYL,FOLDED_AKABAK_REC)

#%%----------------------------------------------------------------------------

""" FOLDED 5 BENDS REC  """

path = os.getcwd()
DATA = np.loadtxt(path + "/DATA/AKABAK SIMU/FOLDED_N5.txt")
DATA_SIDE =  np.loadtxt(path + "/DATA/AKABAK SIMU/FOLDED_5Bends_Side_Rec.txt")

FOLDED_AKABAK = Dict(DATA[:,0], DATA[:,1] + 1j*DATA[:,2] , "AKABAK")
FOLDED_SIDE_AKABAK = Dict(DATA_SIDE[:,0], DATA_SIDE[:,1] + 1j*DATA_SIDE[:,2] , "Folded Side : Leff = 1.04832 m \n L = 1220 m" )

path = "DATA/MEASUREMENTS/MESURES SSA/MESURES MATHIS/SS_LABYRINTHE.csv"
HpLAB = LoadData(path , "Measure")

PlotFRF(HpLAB , FOLDED_AKABAK , SAR=True)

#%%----------------------------------------------------------------------------

""" PROTOTYPE FINAL V3  """
path = os.getcwd()
DATA = np.loadtxt(path + "/DATA/AKABAK SIMU/Pf_Woofer.txt",dtype=str)
DATA = np.char.replace(DATA, ",", ".")
DATA = DATA.astype(float)
f = DATA[:,0]

AKABAK_DUCT = Dict(f, DATA[:,1] + 1j*DATA[:,2], "AKABAK" , color="r", linestyle="-" )

DATA = np.loadtxt(path + "/DATA/AKABAK SIMU/Pf_Woofer_Tweeter.txt",dtype=str)
DATA = np.char.replace(DATA, ",", ".")
DATA = DATA.astype(float)
AKABAK_ALL = Dict(f , DATA[:,1] + 1j*DATA[:,2] , "AKABAK" , color="r" , linestyle="-")


""" Peerless by Tymphany """
#------------------------------------------------------------------------------
P_SDS_135F25CP02_04 = Loudspeaker(f,57.7 , 3.3, 1.22e-3, 7.39, 20.6e-3, 3.89, PouceEnMetre(5.5)/2, 1, BackVolume=19)
PropaP = Propagation(P_SDS_135F25CP02_04, point)
_, HP = FRF(2*np.pi*f/DSP.fs, P_SDS_135F25CP02_04 , PropaP)
#------------------------------------------------------------------------------
P_DX25TG59_04 = Loudspeaker(f, 640, 3.78, 0.02e-3, 2.41 , 0.40e-3, 2.46, 0.0323/2 , 1, BackVolume=False)
PropaP_T = Propagation(P_DX25TG59_04, point)
_, HP_T = FRF(2*np.pi*f/DSP.fs, P_DX25TG59_04 , PropaP_T)

RAD_P = Dict(f, HP, "Woofer" , color="k" , linestyle="-")
RAD_P_T = Dict(f, HP_T, "Tweeter" , color="k" ,linestyle="--")
#------------------------------------------------------------------------------


Fplot(RAD_P , RAD_P_T)

#%%




# P_SDS_135F25CP02_04 = Loudspeaker(f,57.7 , 3.3, 1.22e-3, 7.39, 20.6e-3, 3.89, PouceEnMetre(5.5)/2, 1, BackVolume=9)
# REC = RectangularTube(f, 1.466 , 0.1, 0.1, 0.018)
# RAD = Radiation(f, P_SDS_135F25CP02_04 , REC , 1, 1)
# RAD_PROTOFINAL = Dict(f, RAD.Pvent , label="Python : Straight Duct", color="k")

""" MEASUREMENTS """

""" MEASUREMENTS """
path = "DATA/MESURES PF/PF_WOOFER_ALONE.csv"
Pf_Woofer_Alone = LoadData(path , "Measure : Wood" , color="b")

path = "DATA/MESURES PF/PF_WOOFER_TWEETER_WO_CROSSOVER.csv"
Pf_All_nogain = LoadData(path , "Measure : Wood" , color="b",linestyle="-",PowerSpec=False)

path = "DATA/MESURES PF/PF_WOOFER_TWEETER_WO_CROSSOVER_Gainadjust.csv"
Pf_All = LoadData(path , "Measure : Wood ", color="b")

path = "DATA/MESURES PF/PF_EQ1_WODelay.csv"
Pf_EQ1 = LoadData(path , "Measure : Woofer + Tweeter + EQ1", color="b")

path = "DATA/MESURES PF/PF_WOOFER_ALONE_WO_Foam.csv"
Pf_WO_FOAM = LoadData(path , "Measure : Woofer", color="k")

path = "DATA/MESURES PF/PF_WOOFER_ALONE_Foam1.csv"
Pf_FOAM1 = LoadData(path , "Measure : Wood + foam  (Fibrous)", color="orange")

path = "DATA/MESURES PF/PF_WOOFER_ALONE_Foam2.csv"
Pf_FOAM2 = LoadData(path , "Measure : Wood + foam (Porous)", color="green")

path = "DATA/MESURES PF/PF_WOOFER_FOAM_EQ.csv"
Pf_WT_F_EQ = LoadData(path , "Measure : Woofer + foam + EQ" , PowerSpec=True, color="b")

path = "DATA/MESURES PF/PF_WOOFER_TWEETER_EQ_FOAM_VF1.csv"
Pf_VF1 = LoadData(path , "Measure : Woode + EQ" , PowerSpec=True, color="orange")

"""  PLOT """

PlotFRF(AKABAK_DUCT, Pf_Woofer_Alone ,linewidth=3, loc="lower left")

PlotFRF(AKABAK_ALL,Pf_All_nogain , loc="lower left" , Ylabel=r"$P_{rad}$")

PlotFRF(AKABAK_ALL,Pf_All , loc="lower left" , Ylabel=r"$P_{rad}$")

PlotFRF(Pf_All, Pf_EQ1 , loc="lower left",Ylabel=r"$P_{rad}$")

PlotFRF(Pf_Woofer_Alone,Pf_FOAM1,Pf_FOAM2 , loc="lower left",Ylabel=r"$P_{vent}$")

PlotFRF(AKABAK_DUCT, Pf_WT_F_EQ , loc="lower left",Ylabel=r"$P_{rad}$")

PlotFRF(RAD_P,RAD_P_T ,AKABAK_ALL,Pf_VF1, loc="lower left",Ylabel=r"$P_{rad}$")

PlotFRF(RAD_P,RAD_P_T ,Pf_All,Pf_VF1, loc="lower left",Ylabel=r"$P_{rad}$")



#%%

""" Measurements PLA"""

path = "DATA/MESURES PF/PLA/PF_PLA_WOOFER.csv"
Pf_Woofer_PLA = LoadData(path , "Measure : PLA" , color="c",linestyle="-")

path = "DATA/MESURES PF/PLA/PF_PLA_WOOFER_TWEETER_PHASE_ASSSIGNED.csv"
Pf_All_PLA = LoadData(path , "Measure : PLA" , color="c",linestyle="-",PowerSpec=True)

path = "DATA/MESURES PF/PLA/PLA_WOOFER_TWEETER_FOAM_EQ.csv"
Pf_VF_PLA = LoadData(path , "Measure : PLA" , color="c",linestyle="-",PowerSpec=True)


PlotFRF(Pf_All_PLA,Pf_VF_PLA)

PlotFRF(AKABAK_DUCT,Pf_Woofer_Alone, Pf_Woofer_PLA ,linewidth=3, loc="lower left")
PlotFRF(AKABAK_ALL,Pf_All,Pf_All_PLA , loc="lower left" , Ylabel=r"$P_{rad}$")
PlotFRF(RAD_P,RAD_P_T ,Pf_All,Pf_VF1,Pf_All_PLA, loc="lower left",Ylabel=r"$P_{rad}$")

#%%
""" COMPARISON """

PlotFRF(AKABAK_DUCT, Pf_Woofer_PLA , Pf_Woofer_Alone , linewidth=3, loc="lower left", SAR=True)
PlotFRF(AKABAK_ALL ,Pf_All_PLA,Pf_All, loc="lower right" , Ylabel=r"$P_{rad}$",SAR=True)

#%%----------------------------------------------------------------------------

""" DIRECTIVITY MEASUREMENTS """

deg_val = [0,20,40,60,80,90,100,120,140,160,180]
folder = "DATA/MESURES PF/DIREC/D{}.csv"
freq_val = [30,63, 200,1000,2000,3200]

PlotDirect(folder, deg_val, freq_val , dB=True)
PlotDirect(folder, deg_val, freq_val , dB=None)


