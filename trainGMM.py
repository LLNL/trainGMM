#Copyright (c) 2019, Lawrence Livermore National Security, LLC. All rights reserved.

# This work was produced at the Lawrence Livermore National Laboratory
# (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between
# the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  Copyright is
# reserved to Lawrence Livermore National Security, LLC for purposes
# of controlled dissemination, commercialization through formal
# licensing, or other disposition under terms of Contract 44; DOE
# policies, regulations and orders; and U.S. statutes.  The rights
# of the Federal Government are reserved under Contract 44.


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import numpy as np
from sklearn import mixture
from scipy import stats
import time
import sys
import os
import platform
import gc

TRAINING_PATH = 'sis3316_Cf_V02504730.txt'
OUTPUT_PATH = 'gmmVars.txt'
number_of_pulses = 10000
tailOffset =17

#DO NOT CHANGE IF OUTPUT IS FOR SIS3316 FIRMWARE
baseline_start = 43
baseline_end  = 50
gate_1 = 51
gate_2 = 51+128


def newFile(*args):
    print("NOT CURRENTLY IMPLEMENTED")

def closeFile(*args):
    print("NOT CURRENTLY IMPLEMENTED")

def process_binary_file(input_path,num_pulses,channel):
    val = int(0)
    myWord=b''
    iEvt=0
    channelID=int(0)
    iCounts=[0]*16
    samples = []
    pulse_number =  0
    global progBarVal
    global root
    progBarVal.set(0)
    root.update()
    with open (input_path,"rb") as f:
        myWord=f.read(4)
        while myWord:
            val=int.from_bytes(myWord,byteorder=sys.byteorder)
            if(val==0xdeadbeef):
                f.read(4*3) #NBuffer,NEvents,Length                                                             
                val = int.from_bytes(f.read(4),byteorder=sys.byteorder) #TS,ChID,FMt                            
            elif(val==0x0E0F0E0F): break
            formatBits = (0xF & val)
            channelID = 4*(0x3 & (val>>6))+(0x3 & (val>>4)) #  4*iADC+iChan                                    
            iCounts[channelID] += 1
            f.read(4) #TS, PH                                                                                  
            if(0x1 & formatBits): f.read(4*7) #green section                                                   
            if(0x2 & formatBits): f.read(4*2) #yellow                                                          
            if(0x4 & formatBits): f.read(4*3) #orange                                                          
            val = int.from_bytes(f.read(4),byteorder=sys.byteorder) #word before raw samples                   
            NsampleWords=val & 0x3FFFFFF
            if(channelID==channel):
                samples.append([])
                samples[pulse_number].append(0) #to match unused fields in text
                samples[pulse_number].append(0) #to match unused fields in text  
                for i in range(2*NsampleWords):
                    if(i<gate_2):samples[pulse_number].append(int.from_bytes(f.read(2),byteorder=sys.byteorder))
                    else:f.read(2)
                pulse_number += 1
            else:
                f.read(4*NsampleWords)
            myWord=f.read(4)
            if(pulse_number>num_pulses):break
            if(pulse_number%500==0):
                progBarVal.set(pulse_number/num_pulses*30.0)
                root.update()
            iEvt += 1
    X = np.array(samples,dtype='float')
    print("Counts in each channel")
    print(iCounts)
    print('FILE WAS READ')
    progBarVal.set(30)
    root.update()
    return X

def process_text_file(input_path,num_pulses):
  ##### Parsing the input data
  X = []
  j = 0
  pulse_number =  0
  
  global progBarVal
  global root
  progBarVal.set(0)
  root.update()
  for line in open(input_path,'r'):
        pulse_number = pulse_number + 1
        if(pulse_number%500==0):
            progBarVal.set(pulse_number/num_pulses*30.0)
            root.update()
        if((pulse_number <=num_pulses)):
                X.append([])
                for ii in range(2): X[j].append(int(0)) #Cleanup
                obs = line.rstrip()
                counter = 0
                temp = ''
                
                for i in range(len(obs)):
                        temp = temp + obs[i]
                        if obs[i]==' ':
                                if(counter<3): temp="0"
                                if(counter <= gate_2):X[j].append(int(temp))
                                temp = ''
                                counter = counter +1
                X[j].append(int(temp))
                j = j + 1
        else:
                break

  X = np.array(X,dtype='float')
  print('FILE WAS READ')
  progBarVal.set(30)
  root.update()
  return X

def preprocess_input(dataArray,gate1,gate2,baseline_start,baseline_end):
  # Baseline subtract and Anscombe 
  global progBarVal
  global root
  X = dataArray
  qRatio = []
  qTotal = []
  maxADC = []
  for j in range(len(X[:,0])):
        if(j%500==0):
            progBarVal.set(30+j*60.0/len(X[:,0]))
            root.update()
        maxADC.append(max(X[j,2:]))
        X[j,:] = X[j,:]-np.median(X[j,baseline_start:baseline_end])
        qRatio.append(sum(X[j,gate1+int(spinvalTailOffset.get()):gate2]))
        qTotal.append(sum(X[j,gate1:gate2]))
        for k in range(len(X[j,])):
                if(X[j,k]<0):
                        X[j,k] = 0
                X[j,k] = 2*np.sqrt(X[j,k]+.375)

  # Z-score along each row (pulse)
  for j in range(len(qRatio)):
      if(qTotal[j]!=0): qRatio[j]=qRatio[j]/qTotal[j]
  Y=[]
  minE = int(spinvalMinE.get())
  global qT
  global qR
  qT=[]
  qR=[]

  for j in range(len(qTotal)):
      if((qTotal[j]>minE) and (maxADC[j]<16383)): #Don't train with pulses that are below E cut of saturated
          Y.append(X[j,gate1:gate2])
          qT.append(qTotal[j])
          qR.append(qRatio[j])
  X_zscore = (stats.zscore(Y,axis=1,ddof=1))
  progBarVal.set(90)
  root.update()
  return(X_zscore)

def output_probability_scores(output_path,preprocessed_data,num_pulses): #slight misnomer as first used to write to file
  g = mixture.GaussianMixture(n_components = 2, covariance_type = 'tied',random_state=31410)
  print('TRAINING MODEL')
  g.fit(preprocessed_data)
  predictions_probs = g.predict_proba(preprocessed_data)
  root.update()
  NeutronIndex=1
  GammaIndex=0
  if(g.means_[0,0]<g.means_[1,0]):
      NeutronIndex=0
      GammaIndex=1

  global ax
  ax[0][0].set_xlabel('sample')
  ax[0][0].set_ylabel('scaled amplitude')
  ax[0][0].plot(range(len(g.means_[NeutronIndex])), g.means_[NeutronIndex])
  ax[0][0].plot(range(len(g.means_[GammaIndex])), g.means_[GammaIndex])

  global qT
  global qR
  ax[0][1].scatter(qT,qR,s=0.5,c=predictions_probs[:,NeutronIndex],alpha=0.5)
  ax[0][1].set_ylim(-0.2,0.6)
  ax[0][1].set_xlabel('Integral')
  ax[0][1].set_ylabel('Qtail/Qtotal')

  detConst=1
  for i in range(len(g.means_[0])): detConst = detConst*g.covariances_[i][i]
  detConst = np.log(1.0/np.sqrt(np.power(2.0*np.pi, len(g.means_[0])) * detConst))

  global gmmNmeans
  global gmmGmeans
  global gmmInvCovar
  global gmmDconst
  gmmDconst=detConst
  gmmNmeans=[]
  gmmGmeans=[]
  gmmInvCovar=[]
  for i in range(len(g.means_[NeutronIndex])):
      gmmNmeans.append(g.means_[NeutronIndex,i])
  for i in range(len(g.means_[GammaIndex])):
      gmmGmeans.append(g.means_[GammaIndex,i])
  for i in range(len(g.means_[1])):
      gmmInvCovar.append(1.0/g.covariances_[i,i])

  return (0)

def calc_draw_logVals(preprocessed_data,num_pulses):
  logN=[]
  logG=[]
  global gmmNmeans
  global gmmGmeans
  global gmmInvCovar
  global gmmDconst
  global progBarVal
  global root
  global saveButton
  progBarVal.set(100)
  saveButton.state(['!disabled'])
  root.update()

  for j in range(len(preprocessed_data[:,0])):
      scalarN=0.0
      scalarG=0.0
      for i in range(len(gmmNmeans)):
          scalarN = scalarN + (preprocessed_data[j,i]-gmmNmeans[i])*(preprocessed_data[j,i]-gmmNmeans[i])*gmmInvCovar[i]
          scalarG = scalarG + (preprocessed_data[j,i]-gmmGmeans[i])*(preprocessed_data[j,i]-gmmGmeans[i])*gmmInvCovar[i]

      gmmDconst=0 #We only want deal with the exponent of the likelihood for now
      logN.append(gmmDconst/2.0 - scalarN)
      logG.append(gmmDconst/2.0 - scalarG)

  global logDiff
  logDiff=[]
  logSum=[]
  for i in range(len(logN)):
      logDiff.append((logG[i]-logN[i])/(logN[i]+logG[i]))

  ax[1][1].cla()
  ax[1][1].set_xlabel('Integral')
  ax[1][1].set_ylabel('(LogG-LogN)/(LogN+LogG)')
  ax[1][1].callbacks.connect('xlim_changed', zoomChanged)
  ax[1][1].callbacks.connect('ylim_changed', zoomChanged)
  ax[1][1].scatter(qT,logDiff,s=0.5,c='blue',alpha=0.3)

  return(0)

def selectFile(*args):
    global TRAINING_PATH
    global root
    homeDir = os.path.expanduser("~")
    TRAINING_PATH =  filedialog.askopenfilename(
        initialdir = homeDir,title = "Select file",filetypes = (("text data files","*.txt"),("binary files","*.dat"),("all files","*.*")))
    if(TRAINING_PATH != ''):
        trainButton.state(['!disabled'])
    global chanSB
    chanSB.configure(state="normal")
    if (TRAINING_PATH .split('.')[-1] == "txt"): chanSB.configure(state="disabled")
    print (TRAINING_PATH+' selected')
    root.update()

def zoomChanged(*args):
    global logDiff
    global logSum
    global qT
    global qR
    qTsubset=[]
    qRsubset=[]
    qTother=[]
    qRother=[]
    xLo = ax[1][1].viewLim.bounds[0]
    xHi = xLo + ax[1][1].viewLim.bounds[2]
    yLo = ax[1][1].viewLim.bounds[1]
    yHi = yLo + ax[1][1].viewLim.bounds[3]
    for i in range(len(qT)):
        if((qT[i] > xLo) and (qT[i]<xHi) and (logDiff[i] > yLo) and (logDiff[i] < yHi)):
            qTsubset.append(qT[i])
            qRsubset.append(qR[i])
        else:
            qTother.append(qT[i])
            qRother.append(qR[i])
    ax[1][0].cla()
    ax[1][0].scatter(qTother,qRother,s=0.5, c='black',alpha=1)
    ax[1][0].scatter(qTsubset,qRsubset,s=0.5, c='red',alpha=0.25)
    ax[1][0].set_xlabel('Intergral')
    ax[1][0].set_ylabel('Qtail/Qtotal')
    ax[1][0].set_ylim(-0.2,0.6)

def trainAction(*args):
    global TRAINING_PATH
    global fig
    global ax
    global testButton
    global spinvalChannel
    fig, ax = plt.subplots(2,2)
    fig.set_tight_layout(True)
    fig.set_size_inches(10, 7)
    print ('File '+TRAINING_PATH+' will be used for training')
    number_of_pulses = int(spinvalMaxN.get())
    if (TRAINING_PATH .split('.')[-1] == "txt"): X = process_text_file(TRAINING_PATH,number_of_pulses)
    else: X = process_binary_file(TRAINING_PATH,number_of_pulses,int(spinvalChannel.get()))
    gmm_processed_input = preprocess_input(X,gate_1,gate_2,baseline_start,baseline_end)
    del X
    output_probability_scores(OUTPUT_PATH,gmm_processed_input,number_of_pulses)
    calc_draw_logVals(gmm_processed_input,number_of_pulses)
    del gmm_processed_input
    print ('File '+TRAINING_PATH+' processed')
    testButton.state(['!disabled'])
    plt.show()
    root.update()
    gc.collect()

def exitAction(*args):
    sys.exit()

def testAction(*args):
    global TRAINING_PATH
    global spinvalChannel
    print ('File '+TRAINING_PATH+' will be used for testing')
    number_of_pulses = int(spinvalMaxN.get())
    if (TRAINING_PATH .split('.')[-1] == "txt"): X = process_text_file(TRAINING_PATH,number_of_pulses)
    else: X = process_binary_file(TRAINING_PATH,number_of_pulses,int(spinvalChannel.get()))
    gmm_processed_input = preprocess_input(X,gate_1,gate_2,baseline_start,baseline_end)
    del X
    calc_draw_logVals(gmm_processed_input,number_of_pulses)
    del gmm_processed_input
    print ('File '+TRAINING_PATH+' processed')
    plt.show()
    root.update()
    gc.collect()

def saveButtonAction(*args):
  global gmmNmeans
  global gmmGmeans
  global gmmInvCovar
  global gmmDconst
 
  homeDir = os.path.expanduser("~")
  OUTPUT_PATH =  filedialog.asksaveasfilename(
      initialdir = homeDir,title = "Select file",filetypes = (("text data files","*.txt"),("all files","*.*")))
  if OUTPUT_PATH == '': 
      print("NO OUTPUT FILE SELECTED")
      return ()
  output = open(OUTPUT_PATH,'w')
  
  for i in range(len(gmmNmeans)):
      output.write(str(gmmNmeans[i])+'\n')
  for i in range(len(gmmGmeans)):
      output.write(str(gmmGmeans[i])+'\n')
  for i in range(len(gmmInvCovar)):
      output.write(str(gmmInvCovar[i])+'\n') 
  output.write(str(gmmDconst)+'\n')
  output.close()


###################################################
#This is all the main GUI and loop setup.

if getattr(sys, 'frozen', False):
    os.chdir(sys._MEIPASS)

root = Tk()
root.title("GMM trainer")
root.option_add('*tearOff', FALSE)

menubar = Menu(root)
appmenu = Menu(menubar, name='apple')
menubar.add_cascade(menu=appmenu)
appmenu.add_command(label='About My Application')
appmenu.add_separator()
menu_file = Menu(menubar)
menu_edit = Menu(menubar)
menubar.add_cascade(menu=menu_file, label='File')
menubar.add_cascade(menu=menu_edit, label='Edit')
menu_file.add_command(label='New', command=newFile)
menu_file.add_command(label='Open...', command=selectFile)
menu_file.add_command(label='Exit', command=exitAction)
root['menu'] = menubar


mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)
ttk.Label(mainframe, text="Max Events").grid(column=3, row=3, sticky=W)
ttk.Label(mainframe, text="Min Integral").grid(column=3, row=4, sticky=W)
ttk.Label(mainframe, text="Tail Start").grid(column=3, row=5, sticky=W)
ttk.Button(mainframe, text="Choose File", command=selectFile).grid(column=1, row=1, sticky=W)
trainButton = ttk.Button(mainframe, text="Train", command=trainAction)
trainButton.grid(column=1, row=2, sticky=W)
testButton = ttk.Button(mainframe, text="Test", command=testAction)
testButton.grid(column=3, row=2, sticky=E)
saveButton = ttk.Button(mainframe, text="Save Templates", command=saveButtonAction)
saveButton.grid(column=3, row=1, sticky=E)
spinvalMaxN = StringVar()
Spinbox(mainframe, from_=500.0, to=1000000.0, width=15, textvariable=spinvalMaxN).grid(column=1, row=3, sticky=W, padx=50)
spinvalMaxN.set(number_of_pulses)
spinvalMinE = StringVar()
sb = Spinbox(mainframe, from_=0.0, to=1000000.0, width=15, textvariable=spinvalMinE)
sb.grid(column=1, row=4, sticky=W, padx=50)
spinvalTailOffset = StringVar()
Spinbox(mainframe, from_=0, to=128, width=15,  textvariable=spinvalTailOffset).grid(column=1, row=5, sticky=W, padx=50)
spinvalTailOffset.set(tailOffset)
spinvalChannel = StringVar()
chanSB = Spinbox(mainframe, from_=0, to=15, width=15,  textvariable=spinvalChannel)
chanSB.grid(column=1, row=6, sticky=W, padx=50)
spinvalChannel.set(0)
ttk.Label(mainframe, text="Channel (bin only)").grid(column=3, row=6, sticky=W)

progBarVal=StringVar()
progBarVal.set(25)
progBar = ttk.Progressbar(mainframe, variable=progBarVal, value=0, mode="determinate", length=225, maximum=100)
progBar.grid(column=1, row=7, sticky=W+E, padx=50, columnspan=3)
progBarVal.set(0)
for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
mainframe.columnconfigure(0, weight=1)
mainframe.columnconfigure(1, weight=1,minsize=120)
mainframe.columnconfigure(2, weight=1)
mainframe.rowconfigure(1, weight=1,minsize=35)
mainframe.rowconfigure(2, weight=1,minsize=35)
mainframe.rowconfigure(3, weight=1,minsize=35)

trainButton.state(['disabled'])
testButton.state(['disabled'])
saveButton.state(['disabled'])
chanSB.configure(state="disabled")

fig=0
ax=0

qT=[]
qR=[]

logDiff=[]
logSum=[]

gmmNmeans=[]
gmmGmeans=[]
gmmInvCovar=[]
gmmDconst=0.0

root.mainloop()

#Copyright (c) 2019, Lawrence Livermore National Security, LLC. All rights reserved.                                                         

# This work was produced at the Lawrence Livermore National Laboratory
# (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between  
# the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  Copyright is   
# reserved to Lawrence Livermore National Security, LLC for purposes
# of controlled dissemination, commercialization through formal 
# licensing, or other disposition under terms of Contract 44; DOE
# policies, regulations and orders; and U.S. statutes.  The rights 
# of the Federal Government are reserved under Contract 44.
