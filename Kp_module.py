# -*- coding: utf-8 -*-
"""
Created on Tue May 25 13:35:36 2021

@author: xie(508)
"""

import numpy as np
import Funcs
import time
from scipy import sparse
from scipy.sparse.linalg import eigs

#########################################
def KP_start():
    global a0;  global hbar;    global m0;  global q0; global Time0;
    q0 = 1.602176462e-19 
    hbar = 1.05457266e-34/q0 
    m0 = (9.10956e-31)*(1e-18)/q0 
    Material = ['InAs','GaSb','InSb','AlSb']
    Time0 = time.process_time()    
    Params = Funcs.MaterParam("Material_Parameters.csv")
    Funcs._init()    
    for item in Material:
        Funcs.set_value('ac_'+item,Params[item]["ac"]) 
        Funcs.set_value('av_'+item,Params[item]["av"]) 
        Funcs.set_value('b_'+item,Params[item]["b"]) 
        Funcs.set_value('C11_'+item,Params[item]["C11"])
        Funcs.set_value('C12_'+item,Params[item]["C12"])
        Funcs.set_value('Eg0_'+item,Params[item]["Eg0"]) 
        Funcs.set_value('Eg77_'+item,Params[item]["Eg77"]) 
        Funcs.set_value('Lattice_'+item,0.1*Params[item]["Lattice"]) 
        Funcs.set_value('mc_'+item,m0*Params[item]["mc"]) 
        Funcs.set_value('mcc_'+item,m0*Params[item]["mcc"]) 
        Funcs.set_value('r1_'+item,Params[item]["r1"])
        Funcs.set_value('r2_'+item,Params[item]["r2"])
        Funcs.set_value('r3_'+item,Params[item]["r3"])
        Funcs.set_value('r1c_'+item,Params[item]["r1c"]) 
        Funcs.set_value('r2c_'+item,Params[item]["r2c"]) 
        Funcs.set_value('r3c_'+item,Params[item]["r3c"]) 
        Funcs.set_value('Ep_'+item,Params[item]["Ep"]) 
        Funcs.set_value('spin_'+item,Params[item]["spin"]) 
        Funcs.set_value('VBO_'+item,Params[item]["VBO"]) 
    for item in Material:
        Funcs.MatSetVal(item)
    return Time0

def Structures(T2SL):
    T2SLM = [ ]; T2SLL = [0]; T2SLAPer = []
    for item in T2SL:
        if item == 'Period':
            Period = T2SL[item]
            PeriodL = 0
            for item0 in T2SLL:
                PeriodL += item0
            MNumber = len(T2SLM)
            for j in range(Period-1):
                for i in range(MNumber):
                    T2SLM.append(T2SLM[i])
                    T2SLL.append(T2SLL[i+1])                  
            break
        T2SLM.append(T2SL[item][0])
        T2SLAPer.append(T2SL[item][0])
        a = Funcs.get_value("Lattice_"+T2SL[item][0]) 
        a0 = Funcs.get_value("Lattice_GaSb") 
        az = (2*a0-a) 
        T2SLL.append(0.5*T2SL[item][1]*az) 
    return T2SLM, T2SLL, T2SLAPer, PeriodL, T2SL['Period']

def FDM_initial(T2SLL,PeriodL,Period):
    dz = 0.01
    StepNum = [0]
    T2SLGL = [0]
    for i,item in enumerate(T2SLL,start=1):
        Step = round(T2SLL[i]/(dz))
        StepNum.append(Step+StepNum[i-1])
        T2SLGL.append(T2SLL[i]+T2SLGL[i-1])
        if i == len(T2SLL)-1:
            break
    TotalS = StepNum[-1]
    X = np.linspace(0,PeriodL*Period,TotalS) 
    return TotalS, StepNum, T2SLGL, X, dz    

def BulkBand(T2SLM,StepNum,TotalS,X):
    BulkEc = np.zeros(TotalS)
    BulkEv = np.zeros(TotalS)
    for i,item in enumerate(T2SLM):
        BulkEc[StepNum[i]:StepNum[i+1]] = Funcs.get_value("Ec_"+item) 
        BulkEv[StepNum[i]:StepNum[i+1]] = Funcs.get_value("Ev_"+item) 
    Funcs.set_valueloop("BulkEc",BulkEc)
    Funcs.set_valueloop("BulkEv",BulkEv)

def HamMatrix(T2SLAPer,dz,TotalS):
    global kx;  global ky;    global kt;    global kz;  
    kx = 0; ky = 0; kt = np.sqrt(kx*kx+ky*ky);  kz = 0  
    for item in T2SLAPer:
        MofBulk = Funcs.FDEBulk(kt,item,dz) 
        Funcs.set_valueloop('HamBulk_'+item,MofBulk)
    return kt,kz        

def BulkE_K(material,Startpoint,Endpoint,StepN):
    Funcs._initLoop() 
    HamMatrix(material,0.01,1)
    k_z = np.linspace(Startpoint,Endpoint,StepN)
    BulkE_k = np.zeros((len(k_z),4))
    BulkH_k = Funcs.get_value('BulkH_'+material[0])
    for i,item in enumerate(k_z):
        H_k = BulkH_k[2]*item**2+BulkH_k[1]*item+BulkH_k[0]
        E_k,V_k = np.linalg.eig(H_k)
        E_k = np.sort(E_k.real)
        BulkE_k[i,:] = E_k
    return k_z,BulkE_k

def SLE_V(T1,T2SLM,T2SLAPer,PeriodL,TotalS,StepNum,kt,kz,dz,E_e,BandNum = 8):
    Matrix = sparse.lil_matrix((TotalS*4,TotalS*4),dtype=complex)
    M10,M20,M30 = Funcs.FDEBoundary(T2SLM[-1],T2SLM[0],T2SLM[0],kt,dz) 
    HamS0 = np.concatenate((M20,M30),axis=1)
    Matrix[0:4,0:8] =  HamS0
    Matrix[0:4,4*(TotalS-1):4*TotalS] = np.exp(complex(0,-PeriodL*kz))*M10 
    Step = 0
    StepNumL = []
    Go = True
    for item in StepNum:
        StepNumL.append(item-1)
    for i in range(1,TotalS-1,1):   
        if i in StepNumL:
            M12,M22,M32 = Funcs.FDEBoundary(T2SLM[Step],T2SLM[Step],T2SLM[Step+1],kt,dz)
            Matrix[4*i:4*(i+1),4*(i-1):4*(i+2)] = np.concatenate((M12,M22,M32),axis=1)
            Go = False
        if Go:
            HamM2 = Funcs.get_valueloop('HamBulk_'+T2SLM[Step])
            Matrix[4*i:4*(i+1),4*(i-1):4*(i+2)] = HamM2 
        if i in StepNum:
            Step += 1
            M12,M22,M32 = Funcs.FDEBoundary(T2SLM[Step-1],T2SLM[Step],T2SLM[Step],kt,dz)
            Matrix[4*i:4*(i+1),4*(i-1):4*(i+2)] = np.concatenate((M12,M22,M32),axis=1)
            Go = True
    M13,M23,M33 = Funcs.FDEBoundary(T2SLM[-1],T2SLM[-1],T2SLM[0],kt,dz) 
    HamE3 = np.concatenate((M13,M23),axis=1)
    Matrix[4*(TotalS-1):4*TotalS,4*(TotalS-2):4*TotalS] =  HamE3
    Matrix[4*(TotalS-1):4*TotalS,0:4] = np.exp(complex(0,PeriodL*kz))*M33 
    EigE, EigV = eigs(Matrix,k=BandNum,sigma=E_e,tol=1e-5)         
    return EigE, EigV

def WaveF(T1,T2SLAPer,TotalS,StepNum,BulkEc,BulkEv,EigE,EigV,X,dz,E_e,BandNum = 4):
    EigE = EigE.real
    EigV = EigV.real
    ESort = np.argsort(EigE)
    EigE = EigE[ESort]
    EigV = EigV[:,ESort]    
    WF1 = np.zeros((TotalS,BandNum));  WF2 = np.zeros((TotalS,BandNum))
    WF3 = np.zeros((TotalS,BandNum));  WF4 = np.zeros((TotalS,BandNum))
    for i in range(TotalS):
        WF1[i,:] = EigV[4*i,:] 
        WF2[i,:] = EigV[4*i+1,:] 
        WF3[i,:] = EigV[4*i+2,:] 
        WF4[i,:] = EigV[4*i+3,:] 
    
    Wave_dic = {}
    Wave_dic["WF1"] = WF1;     Wave_dic["WF2"] = WF2
    Wave_dic["WF3"] = WF3;     Wave_dic["WF4"] = WF4
    
    WaveFunc = (np.multiply(WF1,WF1)+np.multiply(WF2,WF2)+np.multiply(WF3,WF3)+np.multiply(WF4,WF4))/dz
    WaveFunc = Funcs.Normalize(WaveFunc)
    WaveFunc = Funcs.SetMatBase(WaveFunc)
    return EigE,WaveFunc,Wave_dic

