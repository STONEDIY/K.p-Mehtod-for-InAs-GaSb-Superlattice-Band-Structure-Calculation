# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 08:40:18 2020

@author: xie(508)
"""
# Electronic band structures and optiacal properities of type-2 superlattice photodectors with interfacial effect. Optical Express. Vol,20,No,2(2012).PengFei Qiao.


import re
import sys
import csv
import numpy as np

#########################################
def _init():
    global _global_dict
    _global_dict = {}
    
def set_value(key,value):
    _global_dict[key] = value
    
def get_value(key, defValue = None):
    try:
        return _global_dict[key]
    except KeyError:
        return key+" Error"

def MaterParam(MatersCsv):
    Maters = {}
    Params = {}
    material = ""
    with open (MatersCsv,'r') as csvfile:
        Info = csv.reader(csvfile,delimiter = ",")
        for lines in Info:
            if len(lines) == 1:
                if len(material) != 0:
                    Maters[material[0]] = Params
                Params = {}
                material = re.findall(r'[A-z]+',lines[0])
            elif len(lines) == 0:
                continue
            else:
                for item in lines:
                    Temp_0 = re.split(r'\s+',item)
                    for i,Temp_1 in enumerate(Temp_0):
                        if len(Temp_1) == 0:
                            del Temp_0[i]
                    if len(Temp_0) == 2:
                        Params[Temp_0[0]] = float(Temp_0[1])
    return Maters

#########################################
global a0;  global hbar;    global m0;  global q0;
q0 = 1.602176462e-19 #electron charge(C); unit energy in eV(J)
hbar = 1.05457266e-34/q0 # reduced Planck const(eV*s)
m0 = (9.10956e-31)*(1e-18)/q0 #kg;ev*s^2/nm^2

def _initLoop():
    global _global_dictloop
    _global_dictloop = {}

def set_valueloop(key,value):
    _global_dictloop[key] = value
    
def get_valueloop(key,defValue = None):
    try:
        return _global_dictloop[key]
    except KeyError:
        return key+" Error"

def ShowPara(material):
    ac = get_value('ac_'+material);        av = get_value("av_"+material)
    b = get_value("b_"+material);          C11 = get_value("C11_"+material)
    C12 = get_value("C12_"+material);      Eg0 = get_value("Eg0_"+material)
    Eg77 = get_value("Eg77_"+material);    a = get_value("Lattice_"+material)
    mc = get_value("mc_"+material);        mcc = get_value("mcc_"+material)
    r1 = get_value("r1_"+material);        r2 = get_value("r2_"+material)
    r3 = get_value("r3_"+material);        r1c = get_value("r1c_"+material)
    r2c = get_value("r2c_"+material);      r3c = get_value("r3c_"+material)
    Ep = get_value("Ep_"+material);        spin = get_value("spin_"+material)
    VBO = get_value("VBO_"+material);
    
    Ev = get_value('Ev_'+material);         Ec = get_value("Ec_"+material)     
    Ae = get_value("Ae_"+material);         At = get_value("At_"+material)
    Pe = get_value("Pe_"+material);         Pt = get_value("Pt_"+material)
    Qe = get_value("Qe_"+material);         Qt = get_value("Qt_"+material)
    Rp = get_value("Rp_"+material);         Sp = get_value("Sp_"+material)         
    Vp = get_value("Vp_"+material);         U = get_value("U_"+material)
    spin = get_value("spin_"+material);     Pcv = get_value('Pcv_'+material)    
    
    HamBulk = get_value("HamBulk_"+material) 
    BulkH = get_value("BulkH_"+material) 
    Para_bulk = {}
    Para_bulk['ac'] = ac;        Para_bulk['av'] = av;      Para_bulk['b'] = b
    Para_bulk['C11'] = C11;      Para_bulk['C12'] = C12;    Para_bulk['Eg0'] = Eg0
    Para_bulk['Eg77'] = Eg77;    Para_bulk['a'] = a;        Para_bulk['mc'] = mc;
    Para_bulk['mcc'] = mcc;      Para_bulk['r1'] = r1;      Para_bulk['r2'] = r2;
    Para_bulk['r3'] = r3;        Para_bulk['r1c'] = r1c;    Para_bulk['r2c'] = r2c;
    Para_bulk['r3c'] = r3c;      Para_bulk['Ep'] = Ep;      Para_bulk['spin'] = spin;
    Para_bulk['VBO'] = VBO;
    return Para_bulk
def MatSetVal(material):
# =============================================================================
     ac = get_value('ac_'+material);        av = get_value("av_"+material)
     b = get_value("b_"+material);          C11 = get_value("C11_"+material)
     C12 = get_value("C12_"+material);      Eg0 = get_value("Eg0_"+material)
     Eg77 = get_value("Eg77_"+material);    a = get_value("Lattice_"+material)
     mc = get_value("mc_"+material);        mcc = get_value("mcc_"+material)
     r1 = get_value("r1_"+material);        r2 = get_value("r2_"+material)
     r3 = get_value("r3_"+material);        r1c = get_value("r1c_"+material)
     r2c = get_value("r2c_"+material);      r3c = get_value("r3c_"+material)
     Ep = get_value("Ep_"+material);        spin = get_value("spin_"+material)
     VBO = get_value("VBO_"+material);      a0 = get_value("Lattice_GaSb")
# =============================================================================
     Ev = VBO
     Ec = Ev + Eg0
     exx = (a0-a)/a; eyy = exx; ezz = (-2*C12/C11)*exx
     Ae = ac*(exx+eyy+ezz); At = hbar*hbar/(2*mcc)
     Pe = -av*(exx+eyy+ezz); Pt = hbar*hbar*r1c/(2*m0)
     Qe = -b/2*(exx+eyy-2*ezz); Qt = hbar*hbar*r2c/(2*m0)
     Pcv = np.sqrt(Ep*hbar*hbar/(2*m0))
     Rp = -(hbar*hbar/(2*m0))*np.sqrt(3)*((r2c+r3c)/2)
     Sp = (hbar*hbar/(2*m0))*2*np.sqrt(3)*r3c
     Vp = Pcv/np.sqrt(6)
     U = Pcv/np.sqrt(3)
# =============================================================================
     set_value('Ev_'+material,Ev);          set_value('Ec_'+material,Ec)
     set_value('exx_'+material,exx);        set_value('eyy_'+material,eyy)
     set_value('ezz_'+material,ezz);        set_value('Ae_'+material,Ae)
     set_value('At_'+material,At);          set_value('Pe_'+material,Pe)
     set_value('Pt_'+material,Pt);          set_value('Qe_'+material,Qe)
     set_value('Qt_'+material,Qt);          set_value('Rp_'+material,Rp)
     set_value('Sp_'+material,Sp);          set_value('Vp_'+material,Vp)
     set_value('U_'+material,U);            set_value('Pcv_'+material,Pcv)
# =============================================================================   
    
def Ham(kt,material):   
# =============================================================================
     Ev = get_value('Ev_'+material);         Ec = get_value("Ec_"+material)     
     Ae = get_value("Ae_"+material);         At = get_value("At_"+material)
     Pe = get_value("Pe_"+material);         Pt = get_value("Pt_"+material)
     Qe = get_value("Qe_"+material);         Qt = get_value("Qt_"+material)
     Rp = get_value("Rp_"+material);         Sp = get_value("Sp_"+material);         
     Vp = get_value("Vp_"+material);         U = get_value("U_"+material);
     spin = get_value("spin_"+material)  
# =============================================================================
     MatrixU0 = np.zeros((4,4),dtype=complex)
     MatrixU0[0] = [Ec+Ae+At*kt*kt,             -np.sqrt(3)*Vp*kt,              -Vp*kt,                         -np.sqrt(2)*Vp*kt]
     MatrixU0[1] = [np.conj(MatrixU0[0,1]),     Ev-Pe-Pt*kt*kt-Qe-Qt*kt*kt,     Rp*kt*kt,                       np.sqrt(2)*Rp*kt*kt]
     MatrixU0[2] = [np.conj(MatrixU0[0,2]),     np.conj(MatrixU0[1,2]),         Ev-Pe-Pt*kt*kt+Qe+Qt*kt*kt,     -np.sqrt(2)*(Qe+Qt*kt*kt)]
     MatrixU0[3] = [np.conj(MatrixU0[0,3]),     np.conj(MatrixU0[1,3]),         np.conj(MatrixU0[2,3]),         Ev-Pe-Pt*kt*kt-spin]
     MatrixU1 = np.zeros((4,4),dtype=complex)
     MatrixU1[0] = [0,                          0,                              complex(0,np.sqrt(2)*U),        complex(0,-U)]
     MatrixU1[1] = [np.conj(MatrixU1[0,1]),     0,                              complex(0,Sp*kt),               complex(0,-Sp*kt/np.sqrt(2))]
     MatrixU1[2] = [np.conj(MatrixU1[0,2]),     np.conj(MatrixU1[1,2]),         0,                              complex(0,-Sp*kt*np.sqrt(3/2))]
     MatrixU1[3] = [np.conj(MatrixU1[0,3]),     np.conj(MatrixU1[1,3]),         np.conj(MatrixU1[2,3]),         0]
     MatrixU2 = np.zeros((4,4),dtype=complex)
     MatrixU2[0] = [At,                         0,                              0,                              0]
     MatrixU2[1] = [np.conj(MatrixU2[0,1]),     -Pt+2*Qt,                       0,                              0]
     MatrixU2[2] = [np.conj(MatrixU2[0,2]),     np.conj(MatrixU2[1,2]),         -Pt-2*Qt,                       -np.sqrt(2)*(-2)*Qt]
     MatrixU2[3] = [np.conj(MatrixU2[0,3]),     np.conj(MatrixU2[1,3]),         np.conj(MatrixU2[2,3]),         -Pt]
# =============================================
     return MatrixU0,MatrixU1,MatrixU2 #eV*nm^2  eV*nm eV
 
def FDEBulk(kt,material,dz):
    MM = Ham(kt,material)
    set_value('BulkH_'+material,MM) 
    MM0 = MM[0];    MM1 = MM[1];    MM2 = MM[2]
    M1Bu = -MM2/(dz*dz)+complex(0,1)*MM1/(2*dz) #eV
    M2Bu = (2*MM2)/(dz*dz)+MM0 #eV
    M3Bu = -MM2/(dz*dz)-complex(0,1)*MM1/(2*dz) #eV
    HamS = np.concatenate((M1Bu,M2Bu,M3Bu),axis=1)
    return HamS

def FDEBoundary(material1,material2,material3,kt,dz):
    MM1 = get_value('BulkH_'+material1)
    MM10 = MM1[0];    MM11 = MM1[1];    MM12 = MM1[2]
    MM2 = get_value('BulkH_'+material2)
    MM20 = MM2[0];    MM21 = MM2[1];    MM22 = MM2[2]
    MM3 = get_value('BulkH_'+material3)
    MM30 = MM3[0];    MM31 = MM3[1];    MM32 = MM3[2]
    M1Bo = -(MM12+MM22)/(2*dz*dz)+complex(0,1)*(MM11+MM21)/(4*dz)
    M2Bo = (MM32+2*MM22+MM12)/(2*dz*dz)+MM20
    M3Bo = -(MM32+MM22)/(2*dz*dz)-complex(0,1)*(MM31+MM21)/(4*dz)
    return M1Bo,M2Bo,M3Bo
    
def EigSort(E,W,ExecValue):
    E0 = abs(E-ExecValue)
    EigS = np.argsort(E0)
    return E[EigS], W[:, EigS]

def EigenEV(EigE,EigV):
    R,C = EigV.shape
    V = np.zeros((R,C))
    for i in range(R):
        V[i,:] = EigV[i,:]+EigE
    return V

def Normalize(Matrix):
    row,col = np.shape(Matrix)
    for i in range(col):
        temp = 0
        for j in range(row):
            temp += Matrix[j,i]**2
        for j in range(row):
            Matrix[j,i] = Matrix[j,i]/(temp**0.5)
    for i in range(col):
        temp = 0
        for j in range(row):
            temp += Matrix[j,i]**2
    return Matrix

def SetMatBase(Matrix):
    row,col = np.shape(Matrix)
    for i in range(col):
        temp = Matrix[0,i]
        Matrix[:,i] = Matrix[:,i]-temp
    return Matrix

def WaveAmplify(Matrix,c=[6,5],am=[10,15]):
    Matrix = 1*Matrix
    Matrix[:,c[0]] = am[0]*Matrix[:,c[0]]
    Matrix[:,c[1]] = am[1]*Matrix[:,c[1]]
    return Matrix

def InAsGaSb_WhichEcEv(StepNum,X,BulkEc,BulkEv,E,Waves):    
    Wrow,Wcol = np.shape(Waves)    
    x11 = round((StepNum[0]+StepNum[1])/2)
    x12 = round((x11+StepNum[2])/2)
    x21 = round((StepNum[1]+StepNum[2])/2)
    x22 = round((x21+StepNum[2])/2)    
    if BulkEc[x11]>BulkEc[x21]:
        for i in range(Wcol):
            try:
                Temp11 = Waves[x11,i]-Waves[x12,i]
                Temp12 = Waves[x12,i]-Waves[x21,i]
                Temp21 = Waves[x21,i+1]-Waves[x22,i+1]
                Temp22 = Waves[x22,i+1]-Waves[x11,i+1]
            except IndexError:
                print('Error: You need to increase BandNums ')
                sys.exit(1)
            if Temp11>0 and Temp12>0:
                if Temp21>0 and Temp22>0:
                    if E[i]>min(BulkEv) and E[i]<max(BulkEv):
                        if E[i+1]>min(BulkEc) and E[i+1]<max(BulkEc):
                            NumEv = i
                            break
    else:
        for i in range(Wcol):
            try:
                Temp11 = Waves[x21,i]-Waves[x22,i]
                Temp12 = Waves[x22,i]-Waves[x11,i]
                Temp21 = Waves[x11,i+1]-Waves[x12,i+1]
                Temp22 = Waves[x12,i+1]-Waves[x21,i+1]
            except IndexError:
                print('Error: You need to increase BandNums ')
                sys.exit(1)
            if Temp11>0 and Temp12>0:
                if Temp21>0 and Temp22>0:
                   if E[i]>min(BulkEv) and E[i]<max(BulkEv):
                        if E[i+1]>min(BulkEc) and E[i+1]<max(BulkEc):
                            NumEv = i
                            break
    return NumEv,NumEv+1    