# -*- coding: utf-8 -*-
"""
Created on Thu May 20 08:24:16 2021

@author: xie(508)
"""

import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext 
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import time
import Kp_module

##########################################################################################################
#Kp method 
#M1 E1 Name material A;M2 E2 Name material B;M3 E3 Layers material A;M4 E4 Layers material B;M5 E5 Periods
#L1 Label Periods;     L2 Bulkband Plot;     T1 State showing       ;M6 E6 Band Numbers;     M7 E7 Expected band gap 
#B1 RUN single click:KP_one_pre(T2SL); double click:KP_one(T2SL,T1)      
def KP_MODULE_START():
    T1.delete('1.0','end')
    F1_x = (1/15)*Window_w;     F1_y = (1/20)*Window_h; F1_w = (41/30)*Window_w;    F1_h = (11/20)*Window_h  #Frame_1 
    E1_x = F1_x;                E1_y = F1_y-30;         E1_w = 110;                 E1_h = 40  #Name material A 
    E2_x = E1_x+150;            E2_y = E1_y;            E2_w = E1_w;                E2_h = E1_h  #Name material B
    E3_x = E1_x;                E3_y = E1_y+60;         E3_w = E1_w;                E3_h = E1_h  #Layers material A 
    E4_x = E1_x+150;            E4_y = E3_y;            E4_w = E1_w;                E4_h = E1_h  #Layers material B
    E5_x = E2_x+180;            E5_y = E1_y+60;         E5_w = E1_w;                E5_h = E1_h  #Periods 
    E6_x = F1_w-6*E1_x;         E6_y = F1_h-4*E1_y;     E6_w = 50;                  E6_h = E1_h  #Band Numbers
    L1_x = E5_x;                L1_y = E1_y;            L1_w = E1_w;                L1_h = E1_h  #Label Periods    
    L2_x = E1_x+20;             L2_y = E1_y+139;        L2_w = (2/3)*Window_w;      L2_h = (1/3)*Window_h  # Bulkband Plot
    T1_x = L1_x+138;            T1_y = E1_y;            T1_w = (1/2)*Window_w;      T1_h = (1/2)*Window_h-40  #State showing
    B1_x = F1_w-4*E1_x;         B1_y = F1_h-4*E1_y-5;   B1_w = 130;                 B1_h = E1_h+10  #RUN button
    B2_x = T1_x+T1_w-43;        B2_y = T1_y+T1_h-23;    B2_w = 0.5*E6_w;            B2_h = 0.5*E1_h  #B2 Clear the state window 
   
    Frame_1.place(x=F1_x, y=F1_y, width=F1_w, height=F1_h)
    M1.set("InAs") # Name material A
    E1.place(x=E1_x, y=E1_y, width=E1_w, height=E1_h)
    M2.set('GaSb') # Name material B
    E2.place(x=E2_x, y=E2_y, width=E2_w, height=E2_h) 
    M3.set('14') # Layers material A 
    E3.place(x=E3_x, y=E3_y, width=E3_w, height=E3_h)
    M4.set('7') # Layers material B 
    E4.place(x=E4_x, y=E4_y, width=E4_w, height=E4_h)
    M5.set('5') # Periods 
    E5.place(x=E5_x, y=E5_y, width=E5_w, height=E5_h)
    M6.set('8') # Band Numbers M6*4
    E6.place(x=E6_x, y=E6_y, width=E6_w, height=E6_h)        
    M8.set('Period')# label Periods 
    L1.place(x=L1_x, y=L1_y, width=L1_w, height=L1_h)
    M9.set('Plot')# label Bandplot
    L2.place(x=L2_x, y=L2_y, width=L2_w, height=L2_h)
    T1.place(x=T1_x, y=T1_y, width=T1_w, height=T1_h)# 1 State showing
    B1.place(x=B1_x, y=B1_y, width=B1_w, height=B1_h) # run
    
    B2.place(x=B2_x, y=B2_y, width=B2_w, height=B2_h) #B2 Clear the state window
    B2.config(image=pic1,relief="flat",bg = 'lightblue')
    
    Temp_out = "T2SL Bandstructure Calculation..."
    T1insert(Temp_out)
    
    B1.bind("<Button-1>",B1_RUN_C1)
    B1.bind("<Double-Button-1>",B1_RUN_C2)    
    B2.bind("<Button-1>",State_Clear)
    
    V_temp_1.set(L2_x);   V_temp_2.set(L2_y);
    V_temp_3.set(L2_w);   V_temp_4.set(L2_h);

def KP_one_pre(T2SL,T2SL_Period = 1):
    Kp_module.KP_start()
    T2SL['Period'] = T2SL_Period
    T2SLM,T2SLL,T2SLAPer,PeriodL,Period = Kp_module.Structures(T2SL)
    TotalS, StepNum, T2SLGL, X, dz = Kp_module.FDM_initial(T2SLL,PeriodL,Period)
    Kp_module.Funcs._initLoop()
    kt,kz = Kp_module.HamMatrix(T2SLAPer,dz,TotalS)
    Kp_module.BulkBand(T2SLM,StepNum,TotalS,X)
    BulkEc = Kp_module.Funcs.get_valueloop("BulkEc")
    BulkEv = Kp_module.Funcs.get_valueloop("BulkEv")
    return T2SLM,T2SLL,T2SLAPer,PeriodL,TotalS,StepNum,T2SLGL,X,dz,kt,kz,BulkEc,BulkEv

def Pre_run(T1,T2SL,T2SLAPer,T2SLM,PeriodL,TotalS,StepNum,X,BulkEc,BulkEv,kt,kz,dz,BandNums,T2SL_Period = 1):
    T2SL['Period'] = T2SL_Period
    E_expected = []
    for item in T2SLAPer:
        E_expected.append(Kp_module.Funcs.get_value("Ec_"+item)) 
    E_e = np.min(E_expected)
    # constructe the 8X8 matrix and obtain eigen values
    EigE, EigV = Kp_module.SLE_V(T1,T2SLM,T2SLAPer,PeriodL,TotalS,StepNum,kt,kz,dz,E_e,BandNum = BandNums)
    # Bandgap energy, Band energies, Wavefunction 
    EigE_sorted,Wave_sorted,Wave_dic = Kp_module.WaveF(T1,T2SLAPer,TotalS,StepNum,BulkEc,BulkEv,EigE,EigV,X,dz,E_e,BandNum = BandNums)
    Kp_module.Funcs.set_valueloop('X',X)
    Kp_module.Funcs.set_valueloop('EigE_sorted',EigE_sorted)
    Kp_module.Funcs.set_valueloop('Wave_sorted',Wave_sorted)
    Kp_module.Funcs.set_valueloop('Wave_dic',Wave_dic)
    # Find hh1 and c1)
    h1,c1 = Kp_module.Funcs.InAsGaSb_WhichEcEv(StepNum,X,BulkEc,BulkEv,EigE_sorted,Wave_sorted)
    Eg = EigE_sorted[c1]-EigE_sorted[h1]    
    return Eg,EigE_sorted,c1  

def KP_one(T2SL,T1):
    Temp_out = "Program starts:"+ str(0) +'s'
    T1insert(Temp_out)
    Time0 = Kp_module.KP_start()
    period = T2SL['Period']
    BandNums = int(E6.get())*period
    if BandNums < 10:
        BandNums = 10
    if period == 1:
        T2SLM,T2SLL,T2SLAPer,PeriodL,TotalS,StepNum,T2SLGL,X,dz,kt,kz,BulkEc,BulkEv = KP_one_pre(T2SL)
        Time1 = time.process_time()
        Temp_out = "Structure constructed: "+str(format(Time1-Time0,'.3f'))+'s'
        T1insert(Temp_out)        
        Time2 = time.process_time() 
        Temp_out = "Hamiltonian matrices constructed:"+ str(format(Time2-Time0,'.3f'))+'s'
        T1insert(Temp_out)
        Temp_out = "Matrix eigen problem solving ..."
        T1insert(Temp_out)
        Eg_pre,EigE_pre,c1band_pre = Pre_run(T1,T2SL,T2SLAPer,T2SLM,PeriodL,TotalS,StepNum,X,BulkEc,BulkEv,kt,kz,dz,BandNums=20,T2SL_Period = 1)
        Time3 = time.process_time()
        Temp_out = "Matrix eigen problem solved"+ str(format(Time3-Time0,'.3f'))+'s'
        T1insert(Temp_out)
        Temp_out = "BandGap:"+ str(format(Eg_pre,'.3f'))+'eV'
        T1insert(Temp_out,colour="snow")
        Temp_out = "Wavelength: "+ str(format(1.24/Eg_pre,'.3f'))+'μm'
        T1insert(Temp_out,colour="snow")
    else:
        T2SLM,T2SLL,T2SLAPer,PeriodL,TotalS,StepNum,T2SLGL,X,dz,kt,kz,BulkEc,BulkEv = KP_one_pre(T2SL)
        Eg_pre,EigE_pre,c1band_pre = Pre_run(T1,T2SL,T2SLAPer,T2SLM,PeriodL,TotalS,StepNum,X,BulkEc,BulkEv,kt,kz,dz,BandNums=20,T2SL_Period = 1)
        Time1 = time.process_time()
        Temp_out = "Structure constructed: "+str(format(Time1-Time0,'.3f'))+'s'
        T1insert(Temp_out)
        T2SLM,T2SLL,T2SLAPer,PeriodL,TotalS,StepNum,T2SLGL,X,dz,kt,kz,BulkEc,BulkEv = KP_one_pre(T2SL,T2SL_Period = period) 
        Time2 = time.process_time() 
        Temp_out = "Hamiltonian matrices constructed:"+ str(format(Time2-Time0,'.3f'))+'s'
        T1insert(Temp_out)
        Temp_out = "Matrix eigen problem solving ..."
        T1insert(Temp_out)
        EigE, EigV = Kp_module.SLE_V(T1,T2SLM,T2SLAPer,PeriodL,TotalS,StepNum,kt,kz,dz,round(EigE_pre[c1band_pre],3),BandNum = BandNums)
        Kp_module.Funcs.set_valueloop('EigE',EigE)
        Kp_module.Funcs.set_valueloop('EigV',EigV)
        Time3 = time.process_time()
        Temp_out = "Matrix eigen problem solved: "+ str(format(Time3-Time0,'.3f'))+'s'
        T1insert(Temp_out)
        EigE_sorted,Wave_sorted,Wave_dic = Kp_module.WaveF(T1,T2SLAPer,TotalS,StepNum,BulkEc,BulkEv,EigE,EigV,X,dz,round(EigE_pre[c1band_pre],3),BandNum = BandNums)
        Temp_out = "BandGap:"+ str(format(Eg_pre,'.3f'))+'eV'
        T1insert(Temp_out,colour="snow")
        Temp_out = "Wavelength: "+ str(format(1.24/Eg_pre,'.3f'))+'μm'
        T1insert(Temp_out,colour="snow")
        Kp_module.Funcs.set_valueloop('X',X)
        Kp_module.Funcs.set_valueloop('EigE_sorted',EigE_sorted)
        Kp_module.Funcs.set_valueloop('Wave_sorted',Wave_sorted)
        Kp_module.Funcs.set_valueloop('Wave_dic',Wave_dic)
        Kp_module.Funcs.set_valueloop('X',X)
        #Kp_module.Funcs.set_valueloop('StepNum',StepNum)
    c1_pre = c1band_pre
    X = Kp_module.Funcs.get_valueloop('X')
    EigE_Sorted =Kp_module.Funcs.get_valueloop('EigE_sorted')
    Wave_sorted = Kp_module.Funcs.get_valueloop('Wave_sorted')
    Wave_dic = Kp_module.Funcs.get_valueloop('Wave_dic')
    BulkEc = Kp_module.Funcs.get_valueloop("BulkEc")
    BulkEv = Kp_module.Funcs.get_valueloop("BulkEv")
    c1 = round(EigE_pre[c1_pre],3)
    for i,item in enumerate(EigE_Sorted):
        if c1 == round(item,3):
            c1 = i
            break
    mini_bands = {}
    mini_bands['c2'] = [c1+period,c1+2*period-1];       mini_bands['c1'] = [c1,c1+period-1]
    mini_bands['hh1'] = [c1-period,c1-1];               mini_bands['lh1'] = [c1-2*period,c1-period-1]
    mini_bands['hh2'] = [c1-3*period,c1-2*period-1];    mini_bands['so'] = [c1-5*period,c1-4*period-1]
    Kp_module.Funcs.set_valueloop('mini_bands',mini_bands)
    mini_bands_value_low = 1.24/(EigE_Sorted[mini_bands['c1'][0]]-EigE_Sorted[mini_bands['hh1'][1]])
    mini_bands_value_upe = 1.24/(EigE_Sorted[mini_bands['c1'][1]]-EigE_Sorted[mini_bands['hh1'][0]])
    Temp_out = "mini-bands: "
    T1insert(Temp_out,colour="snow")
    Temp_out = str(mini_bands_value_low)+"μm ~ "+str(mini_bands_value_upe)+"μm"
    T1insert(Temp_out,colour="snow")
    c1Amp = 20; h1Amp = -10 # Amplify the and c1 and h1 wavefunctions
    WaveFuncEV = np.zeros(np.shape(Wave_sorted))
    WaveFuncEV = Kp_module.Funcs.WaveAmplify(Wave_sorted,c=[c1,c1-1],am=[c1Amp,h1Amp])      
    WaveFuncEV = Kp_module.Funcs.EigenEV(EigE_Sorted,WaveFuncEV)  
    
    Figure2 = Figure()
    Band_plot2 = Figure2.add_subplot(111)
    Band_plot2.plot(X,BulkEc,color='r',linewidth=2)
    Band_plot2.plot(X,BulkEv,color='b',linewidth=2)
    Band_plot2.plot(X,WaveFuncEV[:,c1-1],linewidth=2,color='c',linestyle='-.') #hh
    Band_plot2.plot(X,WaveFuncEV[:,c1],linewidth=2,color='c',linestyle='-.') #c
    Band_plot2.plot([min(X),max(X)],[EigE_Sorted[mini_bands['c1'][0]],EigE_Sorted[mini_bands['c1'][0]]],color='r',linestyle=":")
    Band_plot2.plot([min(X),max(X)],[EigE_Sorted[mini_bands['c1'][1]],EigE_Sorted[mini_bands['c1'][1]]],color='r',linestyle=":")
    Band_plot2.plot([min(X),max(X)],[EigE_Sorted[mini_bands['hh1'][0]],EigE_Sorted[mini_bands['hh1'][0]]],color='b',linestyle=":")
    Band_plot2.plot([min(X),max(X)],[EigE_Sorted[mini_bands['hh1'][1]],EigE_Sorted[mini_bands['hh1'][1]]],color='b',linestyle=":")
    Band_plot2.set_xlabel('Z (nm)',fontsize=14)
    Band_plot2.set_ylabel('Energy (eV)',fontsize=14)
    Top2 = tk.Toplevel()
    Top2.title('Energy Band')
    Top2.geometry("645x550") 
    Top2.resizable(0,0)
    canvas2 = FigureCanvasTkAgg(Figure2,master=Top2)
    canvas2.draw()
    toolbar = NavigationToolbar2Tk(canvas2,Top2)
    toolbar.pack(anchor="sw")      
    toolbar.update()
    canvas2.get_tk_widget().pack(pady = 20)
    
    Figure3 = Figure()
    Band_plot3 = Figure3.add_subplot(111)
    WF1 = Wave_dic["WF1"];      WF2 = Wave_dic["WF2"]
    WF3 = Wave_dic["WF3"];      WF4 = Wave_dic["WF4"]
    A1 = abs(WF1)/np.sqrt(sum(WF1*WF1*dz))
    A2 = abs(WF2)/np.sqrt(sum(WF2*WF2*dz))
    A3 = abs(WF3)/np.sqrt(sum(WF3*WF3*dz))
    #A4 = abs(WF4)/np.sqrt(sum(WF4*WF4*0.01))    
    c2 = A1[:,mini_bands['c2'][1]]
    c1 = A1[:,mini_bands['c1'][0]]
    hh1 = A2[:,mini_bands['hh1'][1]] 
    lh1 = A3[:,mini_bands['lh1'][1]]
    hh2 = A2[:,mini_bands['hh2'][0]]
    #so = A4[:,mini_bands['so'][1]]    
    c1max = 3*max(c1); c1min = min(c1)
    hh1max = 3*max(hh1); hh1min = min(hh1)
    lh1max = 3*max(lh1); lh1min = min(lh1)
    hh2max = 3*max(hh2); hh2min = min(hh2)
    #somax = max(so)    
    Band_plot3.plot(X,c2+c1max,label='c2')
    Band_plot3.plot(X,c1,label='c1')
    Band_plot3.plot(X,hh1-hh1max+c1min,label='hh1')
    Band_plot3.plot(X,lh1-hh1max+c1min-lh1max+hh1min,label='lh1')
    Band_plot3.plot(X,hh2-hh1max+c1min-lh1max+hh1min-hh2max+lh1min,label='hh2')    
    Band_plot3.set_yticks([])
    Band_plot3.legend(fontsize=12,loc='best')
    Band_plot3.set_xlabel('Z (nm)',fontsize=14)
    Band_plot3.set_ylabel('Namornalized WaveFunctions',fontsize=14)
    #plt.plot(X,so-hh1max+c1min-lh1max+hh1min-hh2max+lh1min-somax+hh2min,label='so') 
    Top3 = tk.Toplevel()
    Top3.title('Wavefunction')
    Top3.geometry("645x550") 
    Top3.resizable(0,0)
    canvas3 = FigureCanvasTkAgg(Figure3,master=Top3)
    canvas3.draw()
    toolbar1 = NavigationToolbar2Tk(canvas3,Top3)
    toolbar1.pack(anchor="sw")      
    toolbar1.update()
    canvas3.get_tk_widget().pack(pady = 20)

def showerror():
    Temp_out = "Error:"
    T1.tag_add('tag','end')
    T1.tag_config('tag',foreground='red')
    T1.insert('end',Temp_out+'\n'+'\n','tag')
    T1.see("end")
  
def B1_RUN_C1(event):
    #Which material 
    Material = ["InAs","GaSb","AlSb","GaAs"]
    M1_material = E1.get()
    M2_material = E2.get()
    if M1_material in Material and M2_material in Material:
        Temp1 = True
    else:
        showerror()
        Temp_out = "Please input materials properly \n:[InAs GaSb]"
        T1insert(Temp_out)
        Temp1 = False
    #How long
    try:
        M1_layers = int(E3.get())
        M2_layers = int(E4.get())
        Temp2 = True
    except:
        showerror()
        Temp_out = "Please input layer information properly \n:only integral numbers are acceptable"
        T1insert(Temp_out)
        Temp2 = False
    #Period
    try:
        Periods = int(E5.get())
        Temp3 = True
    except:
        showerror()
        Temp_out = "Please input period information properly \n:only integral numbers are acceptable"
        T1insert(Temp_out)
        Temp3 = False             
    #Simulation~~~~
    if Temp1 and Temp2 and Temp3:
        Temp_out = "Structure:  [" + M1_material+":"+str(M1_layers)+"MLs "+M2_material+":"+str(M2_layers)+"MLs "+"Period:"+str(Periods)+"]"
        T1insert(Temp_out)
        T2SL = {'A':[M1_material,M1_layers], 'B':[M2_material,M2_layers],'Period':Periods}
        T2SLM,T2SLL,T2SLAPer,PeriodL,TotalS,StepNum,T2SLGL,X,dz,kt,kz,BulkEc,BulkEv = KP_one_pre(T2SL,T2SL_Period=Periods)
        E_CV_plot(X,BulkEc,BulkEv)
    else:
        showerror()
        Temp_out = "Please correct the input information ..."
        T1insert(Temp_out)

def B1_RUN_C2(event):
    #Which material 
    Material = ["InAs","GaSb","AlSb","GaAs"]
    M1_material = E1.get()
    M2_material = E2.get()
    if M1_material in Material and M2_material in Material:
        Temp1 = True
    else:
        Temp1 = False
    #How long
    try:
        M1_layers = int(E3.get())
        M2_layers = int(E4.get())
        Temp2 = True
    except:
        Temp2 = False
    #Period
    try:
        Periods = int(E5.get())
        Temp3 = True
    except:
        Temp3 = False
    try:
        BandNumber = int(E6.get())
        Temp4 = True
    except:
        showerror()
        Temp_out = "Please input BandNumber Cofficient properly \n:only integral numbers are acceptable"
        T1insert(Temp_out)
        Temp4 = False
    try:
        if BandNumber > Periods:
            Temp5 = True
        else:
            showerror()
            Temp5 = False
            Temp_out = "Please conform: BandNumber Cofficient > = Periods"
            T1insert(Temp_out)
    except:
        Temp5 = False           
    #Simulation~~~~
    if Temp1 and Temp2 and Temp3 and Temp4 and Temp5:
        T2SL = {'A':[M1_material,M1_layers], 'B':[M2_material,M2_layers],'Period':Periods}
        Temp_out = "Simulation is running ......"
        Gif = tk.PhotoImage(file=r'pic/run.gif')
        B1.config(image=Gif)
        KP_one(T2SL,T1)
        B1.config(image="")    
        
def State_Clear(event):
    T1.delete('1.0','end')

def E_CV_plot(X,BulkEc,BulkEv):
    Fig = Figure()
    E_CV = Fig.add_subplot(111)
    E_CV.plot(X,BulkEc,color='r',linewidth=2)
    E_CV.plot(X,BulkEv,color='b',linewidth=2)
    canvas = FigureCanvasTkAgg(Fig,master=Frame_1)
    canvas.draw()
    canvas.get_tk_widget().place(x=V_temp_1.get(), y=V_temp_2.get(), width=V_temp_3.get(), height=V_temp_4.get())
###########################################################################################################
def Bulk_MODULE_START():
    T1.delete('1.0','end')
    F1_x = (1/15)*Window_w;     F1_y = (1/20)*Window_h; F1_w = (41/30)*Window_w;    F1_h = (11/20)*Window_h  #Frame_1 
    E1_x = F1_x-20;             E1_y = F1_y+50;         E1_w = 110;                 E1_h = 40  #Name material
    E2_x = E1_x;                E2_y = E1_y+60;         E2_w = E1_w;                E2_h = E1_h  #start point
    E3_x = E1_x+120;            E3_y = E1_y+60;         E3_w = E1_w;                E3_h = E1_h  #end point
    E4_x = E2_x+240;            E4_y = E1_y+60;         E4_w = E1_w;                E4_h = E1_h  #stepNumber 
    L2_x = E1_x+370;            L2_y = E1_y-16;         L2_w = (2/3)*Window_w;      L2_h = (2/7)*Window_h+50  # Bulkband Plot
    T1_x = E2_x;                T1_y = E1_y+140;        T1_w = 0.6*Window_w;        T1_h = (1/5)*Window_h-15  #State showing
    B1_x = F1_w-8*E1_x;         B1_y = F1_h-E1_y+1;    B1_w = 130;                 B1_h = E1_h+10  #RUN button
    B2_x = T1_x+T1_w-43;        B2_y = T1_y+T1_h-23;    B2_w = 0.5*50;              B2_h = 0.5*E1_h  #B2 Clear the state window 
    Frame_1.place(x=F1_x, y=F1_y, width=F1_w, height=F1_h)
    M1.set("InAs") # Name material
    E1.place(x=E1_x, y=E1_y, width=E1_w, height=E1_h)
    M2.set('-2') # start point
    E2.place(x=E2_x, y=E2_y, width=E2_w, height=E2_h) 
    M3.set('2') # end point
    E3.place(x=E3_x, y=E3_y, width=E3_w, height=E3_h)
    M4.set('500') # stepNumber
    E4.place(x=E4_x, y=E4_y, width=E4_w, height=E4_h)
    M9.set('Plot')# label Bandplot 
    L2.place(x=L2_x, y=L2_y, width=L2_w, height=L2_h)
    T1.place(x=T1_x, y=T1_y, width=T1_w, height=T1_h)# 1 State showing
    B1.place(x=B1_x, y=B1_y, width=B1_w, height=B1_h) # run
    B2.place(x=B2_x, y=B2_y, width=B2_w, height=B2_h) #B2 Clear the state window
    B2.config(image=pic1,relief="flat",bg = 'lightblue')
    Temp_out = "Bulk Material Bandstructure Calculation..."
    T1insert(Temp_out)
    B1.bind("<Button-1>",B1_RUN_C1_Bulk)
    B1.bind("<Double-Button-1>",B1_RUN_C2_Bulk)    
    B2.bind("<Button-1>",State_Clear)
    V_temp_1.set(L2_x);     V_temp_2.set(L2_y);
    V_temp_3.set(L2_w);     V_temp_4.set(L2_h);

def Bulk_pre(event):
    Kp_module.KP_start()
    material = E1.get()
    Pars = Kp_module.Funcs.ShowPara(material)
    for i,item in enumerate(Pars):
        if i % 2 == 0:
            Temp_out = item+': '+str(Pars[item])+'     '
            T1.insert('end',Temp_out)
            T1.see("end")
        else:
            Temp_out = item+': '+str(Pars[item])
            T1.insert('end',Temp_out+'\n'+'\n')
            T1.see("end")    

def B1_RUN_C1_Bulk(event):
    #Which material
    ifgo = False
    Material = ["InAs","GaSb","AlSb","GaAs"]
    M1_material = E1.get()
    if M1_material in Material:
        Temp1 = True
    else:
        showerror()
        Temp_out = "Please input materials properly \n:[InAs GaSb AlSb GaAs]"
        T1.insert('end',Temp_out+'\n'+'\n')
        T1.see("end")
        Temp1 = False
    #How long
    try:
        Startpoint = float(E2.get())
        Endpoint = float(E3.get())
        StepN = int(E4.get())
        if Endpoint>Startpoint and StepN>3:
            Temp2 = True
        else:
            showerror()
            Temp_out = "Please input the scan parameters properly"
            T1.insert('end',Temp_out+'\n'+'\n')
            T1.see("end")
            Temp2 = False         
    except:
        showerror()
        Temp_out = "Please input the scan parameters properly"
        T1.insert('end',Temp_out+'\n'+'\n')
        T1.see("end")
        Temp2 = False      
    #Simulation~~~~
    if Temp1 and Temp2:
        Temp_out = "Material: "+M1_material+" \nScan: "+str(Startpoint)+"--"+str(Endpoint)+"  (nm-1) \nStep Number: "+str(StepN)
        T1.insert('end',Temp_out+'\n'+'\n')
        T1.see("end")
        Bulk_pre(event)
        ifgo = True
    else:
        showerror()
        Temp_out = "Please correct the input information ..."
        T1.insert('end',Temp_out+'\n'+'\n')
        T1.see("end")
    V_temp_5.set(ifgo);


def B1_RUN_C2_Bulk(event):
    ifgo=V_temp_5.get()
    if ifgo:
        M1_material = [E1.get()]
        Startpoint = float(E2.get())
        Endpoint = float(E3.get())
        StepN = int(E4.get())
        k_z,BulkE_k = Kp_module.BulkE_K(M1_material,Startpoint,Endpoint,StepN)
        E_CV_plot_Bulk(k_z,BulkE_k)
        
def E_CV_plot_Bulk(k_z,BulkE_k):
    Fig = Figure()
    E_CV = Fig.add_subplot(111)
    E_CV.plot(k_z,BulkE_k,linewidth=2)
    canvas = FigureCanvasTkAgg(Fig,master=Frame_1)
    canvas.draw()
    canvas.get_tk_widget().place(x=V_temp_1.get(), y=V_temp_2.get(), width=V_temp_3.get(), height=V_temp_4.get())
    toolbar = NavigationToolbar2Tk(canvas,Frame_1)
    toolbar.place(x=0, y=450,width=800,height=40)      
    toolbar.update()     
  
###########################################################################################################    
    
Window = tk.Tk()
Window_w = 600; Window_h = 900
Screen_w = 0.6*(Window.winfo_screenwidth()-Window_w)
Screen_h = 0.3*(Window.winfo_screenheight()-Window_h)

Window.title('NanoPhotonics')
Window.geometry("%dx%d+%d+%d" % (Window_h,Window_w,Screen_w,Screen_h))
Window.resizable(0,0)

O1_x = 40;                  O1_y = 550;             O1_w = 150;             O1_h = 30 #Optionmenu

Frame_1 = tk.LabelFrame(Window,height=200,width=400)
M1 = tk.StringVar() #M1 E1 Entry1
E1 = tk.Entry(Frame_1,textvariable=M1,font="Helvetica 20 bold",justify="center",width=20)
M2 = tk.StringVar() #M2 E2 Entry2
E2 = tk.Entry(Frame_1,textvariable=M2,font="Helvetica 20 bold",justify="center",width=20)
M3 = tk.StringVar() #M3 E3 Entry3
E3 = tk.Entry(Frame_1,textvariable=M3,font="Helvetica 20 bold",justify="center",width=20)
M4 = tk.StringVar(value='10')  #M4 E4 Entry4
E4 = tk.Entry(Frame_1,textvariable=M4,font="Helvetica 20 bold",justify="center",width=20)
M5 = tk.StringVar()  #M5 E5 Entry5
E5 = tk.Entry(Frame_1,textvariable=M5,font="Helvetica 20 bold",justify="center",width=20)
M6 = tk.StringVar() #M6 E6 Entry6
E6 = tk.Entry(Frame_1,textvariable=M6,font="Helvetica 20 bold",justify="center",width=20)
M7 = tk.StringVar()  #M7 E7 Entry7
E7 = tk.Entry(Frame_1,textvariable=M7,font="Helvetica 20 bold",justify="center",width=20)
M8 = tk.StringVar() #M8 L1 label 1
L1 = tk.Label(Frame_1,textvariable=M8,bg='lightyellow',font="Helvetica 20 bold",justify="center",width=20)#L1 Label Periods 
M9 = tk.StringVar() #M9 L2 label 2
L2 = tk.Label(Frame_1,textvariable=M9,justify="center",bg = 'lightblue')#L2 Bulkband Plot
T1 = scrolledtext.ScrolledText(Frame_1,bg = 'lightblue') #T1 ScrolledText

B1 = tk.Button(Frame_1,text='RUN', bg='Lavender',font="Helvetica 20 bold",justify="center",width=20,relief="raised") #B1 RUN
pic1 = tk.PhotoImage(file = r'pic/DELETE.gif')#Clear the state window 
B2 = tk.Button(Frame_1) #B2 Clear the state window

V_temp_1 = tk.StringVar();    V_temp_3 = tk.StringVar()
V_temp_2 = tk.StringVar();    V_temp_4 = tk.StringVar()
V_temp_5 = tk.StringVar()            

def T1insert(Temp_out,colour=""):
    if len(colour) == 0:
        T1.insert('end',Temp_out+'\n'+'\n')
        T1.see("end")
        T1.update()
    else:                        
        T1.tag_add('tag1','end')
        T1.tag_config('tag1',foreground=colour)
        T1.insert('end',Temp_out+'\n'+'\n',"tag1")
        T1.see("end")
        T1.update()
    
def WidgetDestory(Frame):
    for widget in Frame.winfo_children():
        widget.place_forget() 
    Frame.place_forget()

def NP(event):   
    if Optionmenu.get() == 'Kp':
        WidgetDestory(Frame_1)
        KP_MODULE_START()
    if Optionmenu.get() == 'Bulk':
        WidgetDestory(Frame_1)
        Bulk_MODULE_START()
    if Optionmenu.get() == 'NanoPhotonics':
        WidgetDestory(Frame_1)
        print("1")

#Optionmenu       
V1 = tk.StringVar() 
Optionmenu = ttk.Combobox(Window,textvariable=V1)
Optionmenu['value'] = ('NanoPhotonics','Kp','Bulk')
Optionmenu.current(0)
Optionmenu.bind("<<ComboboxSelected>>",NP)
Optionmenu.place(x=O1_x, y=O1_y, width=O1_w, height=O1_h)

def copyJob():
    T1.event_generate("<<Copy>>")
def pasteJob():
    T1.event_generate("<<Paste>>")

def showPopupMenu(event):
    Popupment.post(event.x_root,event.y_root)

#popupmenu
Popupment = tk.Menu(Frame_1,tearoff=False)
Popupment.add_command(label="Copy",command=copyJob)
Popupment.add_command(label="Paste",command=pasteJob)

T1.bind("<Button-3>",showPopupMenu)

Window.mainloop()









































