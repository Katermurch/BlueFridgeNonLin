import numpy as np
import sys
file_dir= r"C:\Users\crow104\Documents\Python Scripts\sequence_generator"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\data_acquisition_scripts"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Data\2020\201021_thermal_readout\Ramsey vs thermal noise"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\sequence_generator\py_sequences"
if file_dir not in sys.path: sys.path.append(file_dir)

from generator import *
import os
pi = np.pi
import pyvisa
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
#import daq_programs as dp
from tkinter import Tk
import tkinter as tki
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
import wx_programs
import scipy.io as io
import time
from IPython import get_ipython
from scipy import optimize as opt
from matplotlib.patches import Rectangle
import warnings
from scipy.optimize import fsolve


#import data_analysis_pythonRead as dap\

ge_amp_setting = 1.2#1.2#1.9 for 1in driving
ssm_ge_setting = -0.15 #-.162#-0.17 #0.0875#-0.17 #0.0875 #0.170#0.1675 # 0.150 #0.1779 + 0.004#0.1751
pi_ge_time_setting = 25 #33#24#76#17 #176 #102   #50#18 #64
pi2_ge_time_setting = 16#12 #88 #51

ssm_ef_setting = -0.0675#-0.20555
pi_ef_time_setting = 48
pi2_ef_time_setting = 24

ro_pulse_dur =5000 #5000
readout_amp_1 = 1#1
readout_amp_2 =1#1.0#0.3077#0.769#0.134*0.6

mixer_offset = 0#7.7
mixer_offset_ef = 20#25.7

_original_stdout = sys.stdout


global seqInfo

def loading(num_steps = 51, sweeptime =5000):
    #file_length = 20000
    file_length = 16000# 10000 * (int(sweeptime/10000)+2)
#    totlength = rabi_time + 4000
#    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    ringupdown_seq = Sequence(file_length, num_steps)      
#    write_dir = r"C:\Data\2019\encircling\rabi_ef"
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
##END geom


def test_single_generator(ssm_ge = 0.005, ph = 0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 8000
    num_steps = 5
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
    
    ssm_ef = 0.1
    readout_amp = 1 
    
    weaksignal = Pulse(start=file_length-6000, duration=6000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0+ph) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=weaksignal)
    weaksignal.phase = 90+ph
    ringupdown_seq.add_sweep(channel=4,  sweep_name='none', initial_pulse=weaksignal)
    
    #pump = Pulse(start=file_length-6000, duration=1000, amplitude=ge_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
    #ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pump)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-6000-1000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-6000-300:file_length-6000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1.25,1,0.2,0.2])
##END geom
    
def freqsweep_rabi(freq = np.arange(0.1, 0.12, 0.001), singleave=5, ifsave=0):
    
    num_steps = 101
    prob = np.zeros([np.size(freq), num_steps])
    
    for k in np.arange(np.size(freq)):
                
        
        rabi(ssm_ge = freq[k])
        sweepparam, currentprob = readout_python(num_patterns = 101, num_records = 200, 
               sweepparam = np.linspace(0,100,101), ave =singleave, 
               ifprint = 0, ifsave = 0, ifplot=0)

        plt.figure()
        
        plt.plot(sweepparam, currentprob, 'k.-')
        plt.xlabel('Time')
        plt.ylabel('Pe')
        plt.show()

        prob[k, :] = currentprob
        
    if ifsave:
        data = {'prob':prob, 'freq':freq, 
                'sweepparam':sweepparam}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\04_redFridge_test',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        
        
    plt.figure()
    x = freq
    y = sweepparam
    xd = np.append(x-(x[1]-x[0])/2,x[-1]+(x[1]-x[0])/2)
    yd = np.append(y-(y[1]-y[0])/2,y[-1]+(y[1]-y[0])/2)
    X, Y = np.meshgrid(xd, yd)

    print(np.shape(prob))
    plt.pcolormesh(X, Y, prob.T)
    
    plt.xlabel('SSM frequency (GHz)')
    plt.ylabel('Time (ns)')
    
    
def freqsweep_rabi_directChange(freq = np.arange(4.172, 4.1725, 0.001), 
                                singleave=5, ifsave=0, ifloadsequence=1):
    
    num_steps = 101
    prob = np.zeros([np.size(freq), num_steps])
    
    if ifloadsequence:
        rabi()
    
    for k in np.arange(np.size(freq)):
                
        set_qubit_BNC_freq(setfreq=freq[k])
        
        sweepparam, currentprob = readout_python(num_patterns = 101, num_records = 200, 
               sweepparam = np.linspace(0,100,101), ave =singleave, 
               ifprint = 0, ifsave = 0, ifplot=0)

        plt.figure()
        
        plt.plot(sweepparam, currentprob, 'k.-')
        plt.xlabel('Time')
        plt.ylabel('Pe')
        plt.show()

        prob[k, :] = currentprob
        
    if ifsave:
        data = {'prob':prob, 'freq':freq, 
                'sweepparam':sweepparam}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\04_redFridge_test',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        
        
    plt.figure()
    x = freq
    y = sweepparam
    xd = np.append(x-(x[1]-x[0])/2,x[-1]+(x[1]-x[0])/2)
    yd = np.append(y-(y[1]-y[0])/2,y[-1]+(y[1]-y[0])/2)
    X, Y = np.meshgrid(xd, yd)

    print(np.shape(prob))
    plt.pcolormesh(X, Y, prob.T)
    
    plt.xlabel('SSM frequency (GHz)')
    plt.ylabel('Time (ns)')
    
    
def ramsey_ge(num_steps= 101, ifrotphase=1, sweeptime = 15000, ifload = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
   totlength = sweeptime + 4000
   file_length = 10000 * (int(totlength/10000)+1)
   
   if file_length > 80000:
       raise ValueError('File length too long. Make it less than 80000')
       file_length = 80000
       

   ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class 
   
   ge_amp = ge_amp_setting
   ef_amp = 1
   pi_ge_time = pi_ge_time_setting
   pi2_ge_time = pi2_ge_time_setting
   pi_ef_time = pi_ef_time_setting
   pi2_ef_time = pi2_ef_time_setting
   ssm_ef = ssm_ef_setting
   readout_amp = 1 
   oscNum = 6
   
   pi2_ge = Pulse(start=file_length-3050-pi2_ge_time, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
   ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi2_ge)
   pi2_ge.phase = 90
   ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi2_ge)
   
   
   if ifrotphase == 0:
       pi2_ge = Pulse(start=file_length-3050, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
       ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi2_ge)
       pi2_ge.phase = 90
       ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi2_ge)
   elif ifrotphase == 1:
       pi2_ge = Pulse(start=file_length-3050, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
       ringupdown_seq.add_sweep(channel=1, sweep_name='phase', start=0, stop=360*oscNum, initial_pulse=pi2_ge)
       pi2_ge.phase = 90
       ringupdown_seq.add_sweep(channel=2, sweep_name='phase', start=0, stop=360*oscNum, initial_pulse=pi2_ge)
   
   
   
   main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
   ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
   
   ## markers
   alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
   ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
   
   ##create the gate for ch1 an ch2

   ## view output
   if True:
       channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
       channel2_ch = ringupdown_seq.channel_list[1][0]
       channel3_ch = ringupdown_seq.channel_list[2][0]
       channel4_ch = ringupdown_seq.channel_list[3][0]
       marker1 = ringupdown_seq.channel_list[0][2]
       
       channel = channel1_ch + channel3_ch + marker1
       plt.figure()
       plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
       plt.show()
       
       plt.figure()
       plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
       plt.show()
       
   if ifload:
       write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
       ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
       ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])


def wx_ge(num_steps= 101, sweeptime = 30000, ifload = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime 
    file_length = 10000 * (int(totlength/10000)+1)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
                         
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting

    
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1 # 1
    readout_dur = ro_pulse_dur#ro_pulse_dur#8000 #1000
   
    buffer = 50#6000
    
    pi_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    pi_ge.phase = 0
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    #trigger pulse to open switch gate
    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur-500,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000-500, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.358, 1.5, 1.5], offset=[-0.076, 0.025 , 0.09, -0.086])


def repeat_ramsey_ge(repnum = 20, ifsave = 0, ifrunsequence=1):
    
    sweeptime = 15000
    if ifrunsequence:
        ramsey_ge(ifrotphase=1, sweeptime = sweeptime, ifload=1)
    
    totprob = np.zeros([repnum, 101])
    
    for k in np.arange(repnum):
        sweepparam, prob = readout_python(num_patterns = 101, num_records = 200, 
                   sweepparam = np.linspace(0,sweeptime*1e-3,101), ave = 1, 
                   ifprint = 0, ifsave = 0, ifplot=0)
        
        plt.figure()
        plt.plot(sweepparam, prob, 'k.-')
        plt.xlabel('Time ($\\mu$s)')
        plt.ylabel('$P_e$')
        
        plt.title('# {} in {}'.format(k+1, repnum))
        plt.show()
        totprob[k,:] = prob
    
    
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob, 'repnum':repnum}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
        
    plt.figure()
    x = sweepparam
    y = np.arange(repnum)
    xd = np.append(x-(x[1]-x[0])/2,x[-1]+(x[1]-x[0])/2)
    yd = np.append(y-(y[1]-y[0])/2,y[-1]+(y[1]-y[0])/2)
    X, Y = np.meshgrid(xd, yd)

    print(np.shape(totprob))
    plt.pcolormesh(X, Y, totprob)
    
    plt.xlabel('Time ($\\mu$s)')
    plt.ylabel('data #')
    plt.show()
        
def T1_ge(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,q=0, ifload = 1,ROIF=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)

        
    readout_dur = ro_pulse_dur
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    

    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    #Readout
    if q == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2 
    
    #Readout
    if q == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_1, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=50, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE
        
def T1_ge_modulation(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,q=0, ifload = 1,ROIF=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1) #80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        mod_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='width', start=0, stop=-sweeptime,initial_pulse=mod_pulse)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        

    
   
    #Readout
    if q == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2 
    
    #Readout
    if q == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    if ifload:
        
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE

def T1_ge_noise(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,q=0, ifload = 1,ROIF=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1) #80000

    readout_dur = ro_pulse_dur
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        #Adds pi pulse to reach excited state. 
        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        #Adds DL Noise pulse
        noise_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='width', start=0, stop=-sweeptime,initial_pulse=noise_pulse)
        
        #Adds modulation 
        mod_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='width', start=0, stop=-sweeptime,initial_pulse=mod_pulse)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        

    
   
    #Readout
    if q == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2 
    
    #Readout
    if q == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE

def T1_ge_DLnoise_FluxNoise(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,q=0, ifload = 1,ROIF=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1) #80000
#    totlength = t1_time + 4000
#    file_length = 10000 * (int(totlength/10000)+2)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        noise_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='width', start=0, stop=-sweeptime,initial_pulse=noise_pulse)
        
        Flux_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='width', start=0, stop=-sweeptime,initial_pulse=mod_pulse)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        

    
   
    #Readout
    if q == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2 
    
    #Readout
    if q == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE

def T1_ge_noise_only(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,q=0, ifload = 1,ROIF=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1) #80000
#    totlength = t1_time + 4000
#    file_length = 10000 * (int(totlength/10000)+2)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        noise_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='width', start=0, stop=-sweeptime,initial_pulse=noise_pulse)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        

    
   
    #Readout
    if q == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2 
    
    #Readout
    if q == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE

def T1_ge_2q_RO(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,ROIF1=0, ROIF2=0, q=0,ifload = 1 ):
    totlength = sweeptime + 4000
    
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    #file_length = 16000
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur-100, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-100, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    
    
    
    #Readout
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2,ssm_freq=ROIF2,phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE_2q_RO

def T1_ge_Raman_2q_RO(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,ROIF1=0, ROIF2=0, q=0,ifload = 1 ):
    totlength = sweeptime + 4000
    
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    ifsideband = 1
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    #ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
#    pi_ge_time = pi_ge_time_setting

##    ssm_ge = ssm_ge_setting\
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    detun = -0.03 
    ssm_coax = ROIF2 + detun
    amp_coax = 1#*3.317#0#1
    ssm_Q2 = ssm_ge + detun + 0.001 # +0.001 for amp_Q2=0.5 # excitation -detun; decay + detun
    amp_Q2 = 0.5#*3.98#0#0.5 #0.5 #0
    the_seq = Sequence(file_length, num_steps)
    coaxtime = sweeptime#-5000 #500*4
    pulse_len_add = 2000
    
    if q == 0: #qubit 2
               
         # off-resonance drive for RO (of Q1)
        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
#        coax_drive.phase = 90
#        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q0
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1
    
    
        sigma = 40# 200+0*300 #100
        width = int(sigma*5/2)

        sweeptime = coaxtime
        if ifsideband == 1:
            for k in np.arange(num_steps):
                if k >= 0:  # Note (WC): The first pulse with nonzero duration should larger than 2*width

                # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                    pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                    pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                    pulseStart = file_length - readout_dur - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                    pulseStart2 = file_length - readout_dur - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                    tdata = np.arange(pulseLength)
                    tdata2 = np.arange(pulseLength2)
                    pdata1 = np.zeros(pulseLength)
                    pdata2 = np.zeros(pulseLength2)

                    # gaussian shape                    
                    pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
                    pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                    
                    # sine shape
#                    pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata[0:width]) 
#                    pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata[0:width])
                    
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                    secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                    secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                    pdata1[width:secondStop] = amp_Q2
                    pdata2[width:secondStop2] = amp_coax
                    
                    pdata1[secondStop:pulseLength] = (amp_Q2 *
                          np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
                    pdata2[secondStop2:pulseLength2] = (amp_coax *
                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                          
                    
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                    the_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                    
                    the_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
                
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()
        #pi_ge = Pulse(start=file_length-readout_dur-100, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        #ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=0*ge_amp+0.251, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    #        pi_ge.phase = 90+phase_offset
    #        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-100, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    
    
    
    #Readout
    #amplitude = *0 + 0.035
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1*0+0.035,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    #amplitude = *0 + 0.035

    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp_2*0+0.035,ssm_freq=ROIF2,phase=-file_length*ROIF2*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_GE_2q_RO



def T1_M_ge_2q_RO(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,ROIF1=0, ROIF2=0, q=0,ifload = 1 ):
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
#   The goal of this sequence is to stick a readout in the sequence, after an initial pi pulse. 
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    pre_readout_dur = 1000
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur-pre_readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-pre_readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    #Pre-readout: Q1
    main_pulse = Pulse(start = file_length- readout_dur-pre_readout_dur,duration= pre_readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF1, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime,initial_pulse=main_pulse) 
    
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    

    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#END T1_M_GE_2q_RO



def T1_ef(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,ssm_ef=-0.30,pi_ge_time=20,pi_ef_time=20,q=0,ifload = 1,ROIF=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = ge_amp_setting
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
#    ssm_ef = ssm_ef_setting
    readout_amp = 1 
#    oscNum = 6
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q == 0: #0 for qubit 2
        pi_ge = Pulse(start=file_length-readout_dur-0*pi_ge_time-1*pi_ef_time-50, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        pi_ef = Pulse(start=file_length-readout_dur-0*pi_ge_time-0*pi_ef_time-50, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
        pi_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
        
    if q == 1: #0 for qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-0*pi_ge_time-1*pi_ef_time-50, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        # pi_ge.phase = 90+phase_offset
        # ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        pi_ef = Pulse(start=file_length-readout_dur-0*pi_ge_time-0*pi_ef_time-50, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
        # pi_ef.phase = 90+phase_offset
        # ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
    
#    pi_ef = Pulse(start=file_length-readout_dur-pi_ge_time-50, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ef)
#    pi_ef.phase = 90+phase_offset
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ef)
    
#    pi_ge = Pulse(start=file_length-readout_dur-50, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
#    pi_ge.phase = 90+phase_offset
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    #Readout
    if q == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2 
    
    #Readout
    if q == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_1, ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #original readout_amp=1, duration = 1000     -1000
        ringupdown_seq.add_sweep(channel=2,sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=50, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    # #Readout
    # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    # ringupdown_seq.add_sweep(channel=1, marker = 2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    # ## markers
    # alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    # ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
    
#        plt.colorbar()
#        plt.show()
#        plt.plot(channel1_ch[0,31500:32000],'b--o')
#        plt.show()

#        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
        
#end T1_ef
        
def T1_ge_modified(num_steps= 101, sweeptime = 15000, ifload = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(totlength/10000)+1)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    ssm_x=-.1356
    ssm_y=.0542
    phase_offset = mixer_offset
    
    pi_ge = Pulse(start=file_length-readout_dur-50-5, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    pi_ge.phase = 90+phase_offset
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    pi_ge = Pulse(start=file_length-readout_dur-50, duration=0, amplitude=ge_amp, ssm_freq=ssm_x, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    pi_ge.phase = 90+phase_offset
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweeptime, initial_pulse=pi_ge)
      
    
    pi_ge = Pulse(start=file_length-readout_dur-50, duration=0, amplitude=ge_amp, ssm_freq=ssm_y, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    pi_ge.phase = 90+phase_offset
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweeptime, initial_pulse=pi_ge)    
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])


#def ramsey_ef(num_steps= 101, ifrotphase=1, sweeptime = 15000, ifload = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
#    totlength = sweeptime + 4000
#    file_length = 10000 * (int(totlength/10000)+1)
#    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
#        
#                         
##    num_steps = 101
#    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
#    
##    sweep_time = 100
#    ## channels   
#    
#    ge_amp = ge_amp_setting
#    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting
#    ssm_ef = ssm_ef_setting
#    readout_amp = 1 
#    oscNum = 6
#    
#    pi_ge = Pulse(start=file_length-3050-2*pi2_ef_time-pi_ge_time, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#    pi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#    
#    
#    pi2_ef = Pulse(start=file_length-3050-pi2_ef_time-pi_ge_time, duration=-pi2_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi2_ef)
#    pi2_ef.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi2_ef)
#    
#    
#    if ifrotphase == 0:
#        pi2_ef = Pulse(start=file_length-3050-pi_ge_time, duration=-pi2_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi2_ef)
#        pi2_ef.phase = 90
#        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi2_ef)
#        
#    elif ifrotphase == 1:
#        pi2_ef = Pulse(start=file_length-3050-pi_ge_time, duration=-pi2_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='phase', start=0, stop=360*oscNum, initial_pulse=pi2_ef)
#        pi2_ef.phase = 90
#        ringupdown_seq.add_sweep(channel=2, sweep_name='phase', start=0, stop=360*oscNum, initial_pulse=pi2_ef)
#    
#    pi_ge = Pulse(start=file_length-3050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
#    pi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
#    
#    
#    
#    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
#    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#
#    ## view output
#    if True:
#        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = ringupdown_seq.channel_list[1][0]
#        channel3_ch = ringupdown_seq.channel_list[2][0]
#        channel4_ch = ringupdown_seq.channel_list[3][0]
#        marker1 = ringupdown_seq.channel_list[0][2]
#        
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
##        plt.colorbar()
#        plt.show()
#        
#    if ifload:
#        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
#        
#    
#    global seqInfo
#    sweeplabel = 'Time ($\\mu$s)'
#    seqInfo = {'num_steps':num_steps, 'sweepparam':sweeptime,
#                    'sweeplabel':sweeplabel}



def repeat_ramsey_ef(repnum = 20, ifsave = 0, ifrunsequence=1):
    
    sweeptime = 10000
    num_steps = 101
    if ifrunsequence:
        ramsey_ef(num_steps= num_steps, ifrotphase=0, sweeptime = sweeptime, ifload = 1)
    
    totprob = np.zeros([repnum, 101])
    
    for k in np.arange(repnum):
        sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                   sweepparam = np.linspace(0,sweeptime*1e-3,101), ave = 1, 
                   ifprint = 0, ifsave = 0, ifplot=0)
        
        plt.figure()
        plt.plot(sweepparam, prob, 'k.-')
        plt.xlabel('Time ($\\mu$s)')
        plt.ylabel('$P_e$')
        
        plt.title('# {} in {}'.format(k+1, repnum))
        plt.show()
        totprob[k,:] = prob
    
    
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob, 'repnum':repnum}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
        
    plt.figure()
    x = sweepparam
    y = np.arange(repnum)
    xd = np.append(x-(x[1]-x[0])/2,x[-1]+(x[1]-x[0])/2)
    yd = np.append(y-(y[1]-y[0])/2,y[-1]+(y[1]-y[0])/2)
    X, Y = np.meshgrid(xd, yd)

    print(np.shape(totprob))
    plt.pcolormesh(X, Y, totprob)
    
    plt.xlabel('Time ($\\mu$s)')
    plt.ylabel('data #')
    plt.show()
def ch12_leakage(ssb = -0.205, offset=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 8000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
    ## channels  
#    pi_ge=34
#    pi_ef=28
#    pi_hf=50
    pulse_time = 7000
#    ssm_ge = 0.3885
#    ssm_ef = 0.1105
#    ssm_hf = 0.22265
    
    p2_pi_ge = Pulse(start=7000, duration=-pulse_time, amplitude=1, ssm_freq=ssb, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-pulse_time,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+offset
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-pulse_time,initial_pulse=p2_pi_ge)


    #main readout

#    main_pulse = Pulse(start = 7000,duration = 1000, amplitude= 1 )
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
#    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])    
    
def pi_nopi(coef=0,coefpief=0,off=0,pi_ge=24,pi_ef=20,ssm_ge = -0.04575, q=0,r=0, ROIF=0.01): #this is pulsed readout to ring up and ring down cavity dfor e state
    #(coef=0,coefpief=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,pi_ef=20,ssm_ge = -0.04575,ssm_ef=-0.04575-0.15)

    file_length = 16000
    num_steps = 3
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200
    ## channels   
#    pi_ge = pi_ge_time_setting    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_dur = ro_pulse_dur
    phase_offset = mixer_offset
    #phase_offset_ef = mixer_offset_ef
 
    if q==1:
        pi_ge = Pulse(start=file_length-readout_dur-coefpief*pi_ef, duration=-(pi_ge)*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_ge)
        
        
        
        #pi_ef pulse
        pi_ef = gen.Pulse(start=file_length-readout_dur, duration=-pi_ef*coefpief, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_ef)
    # if q==0:
    #     rabi_ge = Pulse(start=file_length-readout_dur, duration=-(pi_ge)*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    #     ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    if q==0:
        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-2*pi_ge,initial_pulse=rabi_ge)
# Different start time    
#    pi_2x = Pulse(start=file_length-readout_dur-coefpief*pi_ef, duration=-(pi_ge)/2*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_2x)

#    rabi_ef = Pulse(start=file_length-readout_dur, duration=-pi_ef*coefpief, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#    rabi_ef.phase = 90+phase_offset_ef
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    
#    rabi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    
     #HET 
     
    #Q1 Readout
    if r == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF, phase=-file_length*ROIF*360) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse) 
        
    #Q2 Readout
    if r == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if False:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
def pi_nopi_noise(coef=0,off=0,pi_ge=24,ssm_ge = -0.04575, q=0,r=0, ROIF=0.01): #this is pulsed readout to ring up and ring down cavity dfor e state
    #(coef=0,coefpief=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,pi_ef=20,ssm_ge = -0.04575,ssm_ef=-0.04575-0.15)

    file_length = 16000
    num_steps = 3
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200
    ## channels   
#    pi_ge = pi_ge_time_setting    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_dur = ro_pulse_dur
    phase_offset = mixer_offset
    #phase_offset_ef = mixer_offset_ef
 
    if q==1:
        rabi_ge = Pulse(start=file_length-readout_dur, duration=-(pi_ge)*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        noise_pulse = Pulse(start=file_length-readout_dur, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=noise_pulse)
    if q==0:
        rabi_ge = Pulse(start=file_length-readout_dur, duration=-(pi_ge)*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)

# Different start time    
#    pi_2x = Pulse(start=file_length-readout_dur-coefpief*pi_ef, duration=-(pi_ge)/2*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_2x)

#    rabi_ef = Pulse(start=file_length-readout_dur, duration=-pi_ef*coefpief, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#    rabi_ef.phase = 90+phase_offset_ef
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    
#    rabi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    
     #HET 
     
    #Q1 Readout
    if r == 1:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF, phase=-file_length*ROIF*360) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse) 
        
    #Q2 Readout
    if r == 0:
        main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if False:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])


def nopi_plus_minus(coef=0,coefES=0,coef_rev=1,off=0,amp_ef=1,pi_ge=24,pi_ef=20,ssm_ge = -0.04575,ssm_ef=-0.04575-0.15,es_phase=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    num_steps = 3
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
#    readout_amp = ro_amp#0.5# 1
    readout_dur = ro_pulse_dur#ro_dur #8000#13000 #1000
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
 
    #pi_ge pulse
    rabi_ge = Pulse(start=file_length-readout_dur-coefES*pi_ef, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90+phase_offset
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #pi_ef/2 pulse to prepare |+-> eigenstate
    rabi_ef = Pulse(start=file_length-(pi_ef/2)*coefES-readout_dur, duration=-(pi_ef/2)*coefES, amplitude=amp_ef, ssm_freq=ssm_ef, phase=0+es_phase) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    rabi_ef.phase = 90+phase_offset_ef+es_phase
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    
    #pi_ef/2 pulse to rotate eigenstate to |e> for |+> (coef_rev=-1) and |f> for |-> (coef_rev=+1)
    rabi_ef = Pulse(start=file_length-readout_dur, duration=-(pi_ef/2)*coefES, amplitude=amp_ef*coef_rev, ssm_freq=ssm_ef, phase=0+es_phase) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    rabi_ef.phase = 90+phase_offset_ef+es_phase
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if False:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

##END geom
def pipi_pi_nopi(coef=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,ssm_ge = -0.04575): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    num_steps = 3
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = ro_amp#0.5# 1
    readout_dur = ro_pulse_dur#ro_dur #8000#13000 #1000
    
 
    
#    rabi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    rabi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1, marker = 2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom    
  
def n_pi_2(coef=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,pi_ef=20,ssm_ge = -0.04575, q=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    num_steps = 2
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = ro_amp#0.5# 1
    readout_dur = ro_pulse_dur#ro_dur #8000#13000 #1000
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q==0: #qubit 2
        pi_2x = Pulse(start=file_length-readout_dur-5, duration=-(pi_ge/2)*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_2x)
    
    if q==1: #qubit 1
        pi_2x = Pulse(start=file_length-readout_dur-5, duration=-(pi_ge/2)*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_2x)

    

    
    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1 ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if False:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
    
def x_measurement(num_steps=51,sweep_time=200,ssm_ge=-0.15,q=0, pi_ge_time=50): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    #    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    #    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
    #    ssm_ge = ssm_ge_setting
       
    readout_amp = 1
    readout_dur = ro_pulse_dur
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-pi_ge_time - sweep_time -5, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        pi_2x = Pulse(start=file_length-readout_dur-pi_ge_time/2 -5, duration= pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=pi_2x)
        
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-pi_ge_time/2 -10, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        pi_2x = Pulse(start=file_length-readout_dur-pi_ge_time/2 -5, duration= pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) 
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_2x)
        
#        rabi_ge = Pulse(start=file_length-readout_dur-pi_ge_time/2 - sweep_time -5, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
#        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        
#        pi_2x = Pulse(start=file_length-readout_dur- -5, duration= -pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) 
#        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_2x)
#        
        
    #    #trigger pulse to open switch gate
    #    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
    #    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    
    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

    
def y_measurement(num_steps=51,sweep_time=200,ssm_ge=-0.15,q=0, pi_ge_time=50): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    #    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    #    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
    #    ssm_ge = ssm_ge_setting
       
    readout_amp = 1
    readout_dur = ro_pulse_dur
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-pi_ge_time - sweep_time -200, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        pi_2x = Pulse(start=file_length-readout_dur-pi_ge_time/2 -100, duration= pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=180) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=pi_2x)
        
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-pi_ge_time - sweep_time -200, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        pi_2x = Pulse(start=file_length-readout_dur-pi_ge_time/2 -100, duration= pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=180) 
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=pi_2x)
        
    #    #trigger pulse to open switch gate
    #    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
    #    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    
    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
    
def rabi_ge(num_steps=51,sweep_time=200,ssm_ge=-0.15,ROIF1=0,ROIF2 = 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        # rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=1, ssm_freq=ssm_ge, phase=0) 
        # ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        # rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=1.9, ssm_freq=ssm_ge, phase=0)
        # ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)


#    
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    #HOMO
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    #HET   
        #Q1 Readout
    #if q == 1:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
   
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=50, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom

def rabi_ge_mod(num_steps=51,sweep_time=200,ssm_ge=-0.15,ROIF1=0,ROIF2 = 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        mod_pulse = Pulse(start=file_length-readout_dur-100, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=mod_pulse)   


#    
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    #HOMO
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    #HET   
        #Q1 Readout
    #if q == 1:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
   
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom



 
def vacuum_rabi(num_steps=51,sweep_time=200,ssm_ge=-0.15,ROIF1=0,ROIF2 = 0,q=0,pi_ge_time=50): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time
    coupler_amp = 1 #1
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        
#        pi_ge = Pulse(start=file_length-readout_dur-20, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge)
#        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=-readout_dur, stop=-readout_dur-sweep_time, initial_pulse=pi_ge)
#        
#        coupler = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=coupler_amp) 
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=-readout_dur, stop=-readout_dur-sweep_time,initial_pulse=coupler)
        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge)
        
        coupler = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=coupler_amp) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=coupler)

   
    if q == 1: #qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-200, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge)

        coupler = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=coupler_amp) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=coupler)
#    
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    #HOMO
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    #HET   
        #Q1 Readout
    #if q == 1:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
   
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
def teaching_test(num_steps=51,sweep_time=200,pi_ge_time=49,ssm_ge=-0.15,ROIF2 = 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 16000
    totlength = sweep_time + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        pi_ge = Pulse(start=file_length-readout_dur-pi_ge_time, duration=-pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time,initial_pulse=pi_ge)
        
        pi_ge=Pulse(start=file_length-readout_dur-pi_ge_time/2, duration=-pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=4*360,initial_pulse=pi_ge)
        
        
    if q == 1: #qubit 1
#        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge)
        pi_ge = Pulse(start=file_length-readout_dur-pi_ge_time, duration=-pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time,initial_pulse=pi_ge)
        
        pi_ge=Pulse(start=file_length-readout_dur-pi_ge_time/2, duration=-pi_ge_time/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=4*360,initial_pulse=pi_ge)
        
#    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF2, phase=0 ) 
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    
##    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
##    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
#    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#
#    ## view output
#    if True:
#        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = ringupdown_seq.channel_list[1][0]
#        channel3_ch = ringupdown_seq.channel_list[2][0]
#        channel4_ch = ringupdown_seq.channel_list[3][0]
#        marker1 = ringupdown_seq.channel_list[0][2]
#        
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,6000:8000], aspect='auto')
#        plt.show()
#        
#    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
def rabi_ge_stest(num_steps=51,sweep_time=200,ssm_ge=-0.15,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        rabi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    #if q == 1: #qubit 1
        #rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        #ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        #rabi_ge.phase = 90+phase_offset
        #ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
    
    
    
    
    
    
    
    
def rabi_ge_orginal(num_steps=51,sweep_time=200,ssm_ge=-0.15,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        rabi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        rabi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
def Always_On(num_steps=1,ssm_ge=-0.15,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    ge_amp = ge_amp_setting
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    phase_offset = mixer_offset
    if q == 0: #qubit 2
        rabi_ge = Pulse(start = 100,duration= file_length-100, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ge)
        rabi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start = 100,duration= file_length-100, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=rabi_ge)
        rabi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=rabi_ge)
   
    
    
    #    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = 100,duration= file_length-100, amplitude= 1 ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])   
def rabi_ge_2qubit(num_steps=51,sweep_time=200,ssm_ge1=-0.15,ssm_ge2=-0.15, ROIF1 =0.1, ROIF2 = 0.1): #for herterodyne
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
#    qubit 1 rabi
    rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge1, phase=0)
    ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    
    #qubit 2 rabi
    rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge2, phase=0)
    ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)

#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker = 2, sweep_name='none',initial_pulse=main_pulse)
    #Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=-file_length*ROIF1*360, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=-file_length*ROIF2*360, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    # main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
 
    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom    
    
def rabi_ef(num_steps=51,sweep_time=200,pi_ge = 20,ssm_ge = -0.150,ssm_ef=-0.300,ef_amp=1,ROIF1= 0,ROIF2= 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#20000# #300 #1000 #3000
    ## channels   
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
#    ef_amp = ge_amp_setting
#    ssm_ef = ssm_ef_setting
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    buffer = 50
    if q == 0: #qubit 2
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        rabi_ef.phase = 90+phase_offset_ef
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
    if q == 1: #qubit 1
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        rabi_ef.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        # second pi_ge-pulse    
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
# #        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
  
#    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker = 1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
    # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    

    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
        
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
def ef_test(num_steps=51,sweep_time=200,pi_ge = 20,ssm_ge = -0.150,ssm_ef=-0.300,ef_amp=1,q=0):
    file_length = 30000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#20000# #300 #1000 #3000
    ## channels   
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
#    ef_amp = ge_amp_setting
#    ssm_ef = ssm_ef_setting
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    buffer = 50
    if q == 0: #qubit 2
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        rabi_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
#        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
    if q == 1: #qubit 1
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        rabi_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
  
#    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker = 1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    

def rabi_J(num_steps=51,sweep_time=200,ef_amp=1,pi_ge = 20,pi_ef=20,ssm_ge = -0.150,ssm_ef=-0.300,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweep_time + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#20000# #300 #1000 #3000
    ## channels   
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
#    ef_amp = ge_amp_setting
#    ssm_ef = ssm_ef_setting
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    buffer = 50
    if q == 0: #qubit 2
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-0*pi_ge-1*pi_ef-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #first pi_ef pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-0*pi_ge-buffer, duration=-pi_ef, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-0*pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        rabi_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
#        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
    if q == 1: #qubit 1
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        rabi_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge_pulse)
        pi_ge_pulse.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
  
#    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1,marker = 2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom

## TO DO rewrite so offsets cumulate sequence
#def rabi_J(num_steps=51,sweep_time=200,ef_amp=1,pi_ge = 20,pi_ef=20,ssm_ge = -0.150,ssm_ef=-0.300,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
#    totlength = sweep_time + 4000
#    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
##    num_steps = 101
#    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
#    
##    sweep_time = 200#20000# #300 #1000 #3000
#    ## channels   
#    
#    ge_amp = ge_amp_setting
##    ssm_ge = ssm_ge_setting
##    pi_ge = pi_ge_time_setting
##    ef_amp = ge_amp_setting
##    ssm_ef = ssm_ef_setting
#    phase_offset = mixer_offset
#    readout_amp = 1#0.5# 1
#    readout_dur = ro_pulse_dur#8000 #13000 #1000
#    buffer = 50
#    if q == 0: #qubit 2
#        #first pi_ge pulse
#        pi_ge_pulse = Pulse(start=file_length-readout_dur-0*pi_ge-1*pi_ef-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        
#        #first pi_ef pulse
#        pi_ge_pulse = Pulse(start=file_length-readout_dur-0*pi_ge-buffer, duration=-pi_ef, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        
#        #drive rabi e-f
#        rabi_ef = Pulse(start=file_length-readout_dur-0*pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        rabi_ef.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        
#        #second pi_ge-pulse    
##        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
##        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge_pulse)
##        pi_ge_pulse.phase = 90+phase_offset
##        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
#        
#    if q == 1: #qubit 1
#        #first pi_ge pulse
#        pi_ge_pulse = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        
#        #drive rabi e-f
#        rabi_ef = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        rabi_ef.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        
#        #second pi_ge-pulse    
#        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
#  
##    #Readout
##    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
##    ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
#    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#
#    ## view output
#    if True:
#        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = ringupdown_seq.channel_list[1][0]
#        channel3_ch = ringupdown_seq.channel_list[2][0]
#        channel4_ch = ringupdown_seq.channel_list[3][0]
#        marker1 = ringupdown_seq.channel_list[0][2]
#        
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
#    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
#    
###END geom

def rabi_ef_ES(ssm_ge = .3925,ssm_ef= .1, num_steps = 51, amp_ge = 0.5, amp_ef=0.5, amp_rabi = 0.5, rabi_time = 100,off_set=0,pi_ge=34,pi_ef=0,ph=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000    
    readout_dur = ro_pulse_dur
    orth_ph = mixer_offset
    orth_ph_ef = mixer_offset_ef
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    #first pi-pulse
    pi_pulse = Pulse(start=file_length-readout_dur-pi_ef/2-3*5, duration=-pi_ge, amplitude=amp_ge, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance

    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time, initial_pulse=pi_pulse)
    pi_pulse.phase = 90+orth_ph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time, initial_pulse=pi_pulse)

    #pulse to prepare eigenstate
    pi_pulse = Pulse(start=file_length-readout_dur-2*5, duration=-pi_ef/2, amplitude=amp_ef, ssm_freq=ssm_ef, phase=0+ph) #pulse is also a class p is an instance
    
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time, initial_pulse=pi_pulse)
    pi_pulse.phase = 90+ph+orth_ph_ef
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time, initial_pulse=pi_pulse)

#    #Rabi Pulse
    rabief = Pulse(start=file_length-readout_dur-5, duration=0, amplitude=amp_rabi, ssm_freq=ssm_ef, phase=90) #pulse is also a class p is an instance
    
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabief)
    rabief.phase = 180+orth_ph_ef
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabief)
    
    #Rabi Pulse
#    rabief = Pulse(start=file_length-readout_dur-5, duration=-rabi_time, amplitude=amp_rabi, ssm_freq=ssm_ef, phase=90) #pulse is also a class p is an instance
#    
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabief)
#    rabief.phase = 180+orth_ph_ef
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabief)
    
#    #main readout
#    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate

    ## view output
#    if True:
#        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = the_seq.channel_list[1][0]
#        channel3_ch = the_seq.channel_list[2][0]
#        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,16000:17000], aspect='auto', extent=[16000,17000,200,0])
#        plt.plot(channel1_ch[3,16800:17000],'b--o')
#        plt.plot(channel2_ch[50,16800:17000],'r--o')
#        plt.show()
#        
    ## write output
    write_dir = r"C:\Data\2023\Jarzynski\sequences"
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
#    wx_programs.wx_set_and_amplitude_and_offset()
##END geom

#def 2_qubit_engine(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,pi_ge_time=20,coupled_mA ifload = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
#    totlength = sweeptime + 4000
#    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
#    
##    if file_length > 80000:
##        raise ValueError('File length too long. Make it less than 80000')
##        file_length = 80000
#        
#    readout_dur = ro_pulse_dur#8000
##    num_steps = 101
#    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
#    
##    sweep_time = 100
#    ## channels   
#    
#    ge_amp = ge_amp_setting
#    ef_amp = 1
##    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
###    ssm_ge = ssm_ge_setting\
#    ssm_ef = ssm_ef_setting
#    readout_amp = 1 
#    oscNum = 6
#    phase_offset = mixer_offset
#    
#    ##pi pulse prepare qubit 1 in |1>
#    pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#    pi_ge.phase = 90+phase_offset
#    ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#    
#    ##bring qubit 2 to qubit 1 frequency to turn on coupling
#    keithley2401.set_current(coupled_mA, step_size_mA=0.001) ## use awg to send a small pulse to ffl
#     
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
#    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#
#    ## view output
#    if True:
#        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = ringupdown_seq.channel_list[1][0]
#        channel3_ch = ringupdown_seq.channel_list[2][0]
#        channel4_ch = ringupdown_seq.channel_list[3][0]
#        marker1 = ringupdown_seq.channel_list[0][2]
#        
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
##        plt.colorbar()
#        plt.show()
#        
#    if ifload:
#        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END engine

def dynamics_variables_EP(ssm_ge = .3925,ssm_ef= .05, num_steps = 51,amp_ge = 1, amp_ES=1,amp_tomo=1, amp_factor = 1, amp_max = 0.5, 
                          detun = 0.0, period = 1000, rabi_time = 1000, off_set=0, pi_ge=24, pi_ES = 30, pi_tomo = 30, phase_tomo = 90, ph_ini = 0,phase_es=0,delta_ef=.00): #this is pulsed readout to ring up and ring down cavity dfor e state
#    totlength = rabi_time + 4000
#    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    file_length = 18000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    readout_dur = ro_pulse_dur
    
    sph=0
#    off_set=0
    ortho_ph = mixer_offset
    orth_ph_ef = mixer_offset_ef

    #first pi-pulse
    pi_pulse = Pulse(start=file_length-readout_dur-pi_tomo/2-pi_ES/2-4*5, duration=-pi_ge, amplitude=amp_ge, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop= -rabi_time,initial_pulse=pi_pulse)
    pi_pulse.phase = 90 + ortho_ph
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop= -rabi_time,initial_pulse=pi_pulse)    
    
    #pulse for prepare initial state
    pi_pulse = Pulse(start=file_length-readout_dur-pi_tomo/2-3*5, duration=-pi_ES/2, amplitude=amp_ES, ssm_freq=ssm_ef, phase=0+phase_es) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop= -rabi_time, initial_pulse=pi_pulse)
    pi_pulse.phase = 90+phase_es + orth_ph_ef
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop= -rabi_time, initial_pulse=pi_pulse)
    
    #variable pulse
    pt_pi_ge_r = Pulse(start=file_length-readout_dur-pi_tomo/2-2*5, duration=0, amplitude=amp_max, ssm_freq=ssm_ef+delta_ef, phase=90,phase_ini=ph_ini, t_loop=period, ff=1,detun_NH=detun,jmin=amp_factor)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1,sweep_name='width', start=0, stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180+sph + orth_ph_ef
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop= -rabi_time,initial_pulse=pt_pi_ge_r)    
    
#    if delta_ef == 0:
        
        #unwrapping and tomography pulse for measurement axis
    p2_pi_ge = Pulse(start=file_length-readout_dur-5, duration=-pi_tomo/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=period,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo+sph + orth_ph_ef
    ringupdown_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
        
#    else:
#        
#        #constant detuning
#        p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=90+phase_tomo,phase_ini=ph_ini, t_loop=period,ff=None,detunlinear=delta_ef,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=1, sweep_name='phase_linear_detun', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
#        p2_pi_ge.phase = 0+phase_tomo+sph
#        the_seq.add_sweep(channel=2,  sweep_name='phase_linear_detun'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)


    #main readout
#    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,16000:17000], aspect='auto', extent=[16000,17000,200,0])
    #plt.plot(channel1_ch[3,:],'b--o')
#        plt.plot(channel2_ch[50,16800:17000],'r--o')
#        plt.show()
#        
    ## write output
    write_dir = r"C:\Data\2023\Jarzynski\sequences"
    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom    
    return channel1_ch

def const_Jdelta_drive_EP(ssm_ge = .3925,ssm_ef= .05, num_steps = 51,amp_ge = 1, amp_ES=1,amp_tomo=1, amp_max = 0.5, rabi_time = 1000, 
                          off_set=0, pi_ge=24, pi_ES = 30, pi_tomo = 30, phase_tomo = 90,phase_es=0,delta_ef=.00): #this is pulsed readout to ring up and ring down cavity dfor e state
#    totlength = rabi_time + 4000
#    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    file_length = 18000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    readout_dur = ro_pulse_dur
    
#    off_set=0
    ortho_ph = mixer_offset
    orth_ph_ef = mixer_offset_ef

    #first pi-pulse
    pi_pulse = Pulse(start=file_length-readout_dur-pi_tomo/2-pi_ES/2-4*5, duration=-pi_ge, amplitude=amp_ge, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop= -rabi_time,initial_pulse=pi_pulse)
    pi_pulse.phase = 90 + ortho_ph
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop= -rabi_time,initial_pulse=pi_pulse)    
    
    #pulse for prepare initial state
    pi_pulse = Pulse(start=file_length-readout_dur-pi_tomo/2-3*5, duration=-pi_ES/2, amplitude=amp_ES, ssm_freq=ssm_ef, phase=0+phase_es) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop= -rabi_time, initial_pulse=pi_pulse)
    pi_pulse.phase = 90+phase_es + orth_ph_ef
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop= -rabi_time, initial_pulse=pi_pulse)
    
    #variable pulse
    pt_pi_ge_r = Pulse(start=file_length-readout_dur-pi_tomo/2-2*5, duration=0, amplitude=amp_max, ssm_freq=ssm_ef+delta_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1,sweep_name='width', start=0, stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180 + orth_ph_ef
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop= -rabi_time,initial_pulse=pt_pi_ge_r)    
    
#    if delta_ef == 0:
        
        #unwrapping and tomography pulse for measurement axis
    p2_pi_ge = Pulse(start=file_length-readout_dur-5, duration=-pi_tomo/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo + orth_ph_ef
    ringupdown_seq.add_sweep(channel=2,  sweep_name='start'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
        
#    else:
#        
#        #constant detuning
#        p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=90+phase_tomo,phase_ini=ph_ini, t_loop=period,ff=None,detunlinear=delta_ef,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=1, sweep_name='phase_linear_detun', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
#        p2_pi_ge.phase = 0+phase_tomo+sph
#        the_seq.add_sweep(channel=2,  sweep_name='phase_linear_detun'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)


    #main readout
#    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate

    ## view output
#    if True:
#        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = the_seq.channel_list[1][0]
#        channel3_ch = the_seq.channel_list[2][0]
#        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,16000:17000], aspect='auto', extent=[16000,17000,200,0])
#        plt.plot(channel1_ch[3,16800:17000],'b--o')
#        plt.plot(channel2_ch[50,16800:17000],'r--o')
#        plt.show()
#        
    ## write output
    write_dir = r"C:\Data\2023\Jarzynski\sequences"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom    
    return the_seq



##BEGIN SPECTROSCOPY Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF2,which_qubit)

def spectroscopy_ge(num_steps,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns

    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

    readout_dur = ro_pulse_dur#2000 #13000 #1000

    
    if q == 1: #qubit 1
        #drive through RO
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)

        # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        # ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
  
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
        # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        # ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)

#READOUT:::HET    

#    Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq

def spectroscopy_pump(num_steps,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns

    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

    readout_dur = ro_pulse_dur#2000 #13000 #1000
    ge_amp = ge_amp_setting
    
    if q == 1: #qubit 1
        
        #qubit freq
        # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
        #pump freq
        pump_pulse = Pulse(start=file_length-readout_dur-10, duration=-700, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=pump_pulse)
  
    if q == 0: #qubit 2
    
        pi_ge = Pulse(start=file_length-readout_dur-720, duration=-59.9431, amplitude=ge_amp, ssm_freq=-.111, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=0, initial_pulse=pi_ge)
        # pi_ge_2 = Pulse(start=file_length-readout_dur-520, duration=-73, amplitude=ge_amp, ssm_freq=-0.152, phase=0) #pulse is also a class p is an instance
        # ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=0, initial_pulse=pi_ge_2)
        #qubit freq
        # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
        #pump freq
        pump_pulse = Pulse(start=file_length-readout_dur-10, duration=-700, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=pump_pulse)

#READOUT:::HET    

#    Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq
##END geom


def tomo_pump(num_steps,off=0,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,pi_ge_time_q1=148,pi_ge_time_q2=73,pi_q1=0,pi_q2=0,ssm_geq1=-.110,ssm_geq2=-.152,pump=0,pump_time=500): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns

    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

    readout_dur = ro_pulse_dur#2000 #13000 #1000
    ge_amp = ge_amp_setting
    
    # if q == 1: #qubit 1
        
    #     #qubit freq
    #     # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
    #     # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
    #     #pump freq
    #     pump_pulse = Pulse(start=file_length-readout_dur-10, duration=-500, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
    #     ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=pump_pulse)
  
    # if q == 0: #qubit 2
    
    pi_ge = Pulse(start=file_length-readout_dur-pump_time-20, duration=-pi_ge_time_q1*pi_q1, amplitude=ge_amp, ssm_freq=ssm_geq1, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=0, initial_pulse=pi_ge)
    pi_ge_2 = Pulse(start=file_length-readout_dur-pump_time-20, duration=-pi_ge_time_q2*pi_q2, amplitude=ge_amp, ssm_freq=ssm_geq2, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=0, initial_pulse=pi_ge_2)
        #qubit freq
        # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
        #pump freq
    pump_pulse = Pulse(start=file_length-readout_dur-10, duration=-pump_time*pump, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=pump_pulse)

#READOUT:::HET    

#    Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        # channel = channel1_ch + channel3_ch + marker1
        # plt.figure()
        # plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        # plt.show()
        
        # plt.figure()
        # plt.imshow(channel[:,6000:8000], aspect='auto')
        # plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    # ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq
##END geom
def pump_versus_time(num_steps,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0, pump_freq = -.075, pulse_time = 500): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns

    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

    readout_dur = ro_pulse_dur#2000 #13000 #1000
    ge_amp = ge_amp_setting
    
    # if q == 1: #qubit 1
        
    #     #qubit freq
    #     # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
    #     # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
    #     #pump freq
    #     pump_pulse = Pulse(start=file_length-readout_dur-10, duration=-500, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
    #     ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=pump_pulse)
  
    if q == 0: #qubit 2
    
        pi_ge = Pulse(start=file_length-readout_dur-20, duration=-148, amplitude=ge_amp, ssm_freq=-.110, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-pulse_time, initial_pulse=pi_ge)
        # pi_ge_2 = Pulse(start=file_length-readout_dur-20, duration=-73, amplitude=ge_amp, ssm_freq=-0.152, phase=0) #pulse is also a class p is an instance
        # ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-500, initial_pulse=pi_ge_2)
        #qubit freq
        # rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
        #pump freq
        pump_pulse = Pulse(start=file_length-readout_dur-10, duration=0, amplitude=spec_amp, ssm_freq=pump_freq, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-pulse_time,initial_pulse=pump_pulse)

#READOUT:::HET    

#    Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq











def spectroscopy_ge_TekAWG520(num_steps,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    
#    ssm_ge = ssm_ge_setting
   
#    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#2000 #13000 #1000
    
#    phase_offset = mixer_offset
    
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
##    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
        
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    #HET - removed single qubit readout kwm
   # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
   # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
##    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    # ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq
def spectroscopy_ge_modulation(num_steps=101,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    
#    ssm_ge = ssm_ge_setting
   
#    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#2000 #13000 #1000
    
#    phase_offset = mixer_offset
    
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
       
        mod_pulse = Pulse(start=file_length-readout_dur-10, amplitude=1, duration=-10000, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='none', initial_pulse=mod_pulse)
        
       
        
##    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
        
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    #HET - removed single qubit readout kwm
   # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
   # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
##    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq
##END geom


def spectroscopy_ge_DLnoise(num_steps=101,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    
#    ssm_ge = ssm_ge_setting
   
#    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#2000 #13000 #1000
    
#    phase_offset = mixer_offset
    
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
       
        noise_pulse = Pulse(start=file_length-readout_dur-10, amplitude=1, duration=-10000, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='none',initial_pulse=noise_pulse)
        mod_pulse = Pulse(start=file_length-readout_dur-10, amplitude=1, duration=-10000, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='none', initial_pulse=mod_pulse)
        
       
        
##    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-10000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
        
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    #HET - removed single qubit readout kwm
   # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
   # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
##    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq


def spectroscopy_ge_coupler_switch(num_steps=101,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF1 =0,ROIF2 =0,q=0,time_coupler=2000): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    
#    ssm_ge = ssm_ge_setting
    sweep_time_ssm = 5000
    coupler_amp =-1
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#2000 #13000 #1000
    
    phase_offset = mixer_offset
    ##    
    if q == 0: #qubit 2
        coupler = Pulse(start=file_length-readout_dur-sweep_time_ssm-time_coupler, duration=+readout_dur+sweep_time_ssm+time_coupler, amplitude=coupler_amp) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=coupler)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-sweep_time_ssm, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    
    if q == 1: #qubit 1
        coupler = Pulse(start=file_length-readout_dur-sweep_time_ssm-time_coupler, duration=+readout_dur+sweep_time_ssm+time_coupler, amplitude=coupler_amp) 
        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=coupler)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-sweep_time_ssm, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)

     

        
        
        
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    #HET - removed single qubit readout kwm
   # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
   # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq
##END geom
    
    
def spectroscopy_ef(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    sweep_time = 200#200#6000 #300 #1000 #3000
    # totlength = sweep_time + 4000
    # file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    file_length = 16000

#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    ## channels       
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    

    #first pi_ge pulse
    if q == 1: #qubit 1
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
    #    pi_ge_pulse.phase = 90+phase_offset
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    #    rabi_ge.phase = 90+phase_offset_ef
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#    
    
    if q == 0: #qubit 2
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
    #    pi_ge_pulse.phase = 90+phase_offset
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    #    rabi_ge.phase = 90+phase_offset_ef
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    
    
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
##    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=50, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq

def spectroscopy_ROcavity(num_steps=3,ssm_ge=-.15,pi_ge_time=20,spec_amp=.05,estate = 1, start_freq=0.1,stop_freq=0.1,ROIF1=0.1,ROIF2 = 0.1,q = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 10000
    the_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

   
#    readout_amp = spec_amp #1#0.5# 1  
    readout_dur = ro_pulse_dur#2000 #13000 #1000
    phase_offset = mixer_offset
    
    ge_amp = ge_amp_setting
#    ge_amp = 0.1

    if estate == 1:
        if q == 1: #qubit 1
            rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            the_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=rabi_ge)

    #    
        if q == 0: #qubit 2
            rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            the_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=rabi_ge)
        
    #Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1*spec_amp,ssm_freq=ROIF1, phase=-file_length*ROIF1*360)
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
 
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2*spec_amp,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[2][1]
        
#        plt.plot(marker1[0])
#        plt.show()
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,6000:8000], aspect='auto')
#        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom

def parametric_coupling(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0):
    
    sweep_time = 1000#200#6000 #300 #1000 #3000
    totlength = sweep_time + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    ## channels       
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q==0:
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
                
        parametric_12 = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=parametric_12)
    if q ==1:
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge_pulse)
                
        parametric_12 = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=parametric_12)
        
    
    
     #HET
    #main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    #ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    #Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
  
   
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
       ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
        ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom
    
  
def parametric_coupling_time_domain(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_para=0,spec_amp=.5,ROIF1=0,ROIF2=0,q=0,sweep_time=0):
    
    #sweep_time = 1000#200#6000 #300 #1000 #3000
    totlength = sweep_time + 4000
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    #num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    ## channels       
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q==0:
        pi_ge_pulse = Pulse(start=file_length-readout_dur-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=4, sweep_name='start',start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=0, amplitude=spec_amp, ssm_freq=ssm_para, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    if q==1:
        pi_ge_pulse = Pulse(start=file_length-readout_dur-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='start',start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=0, amplitude=spec_amp, ssm_freq=ssm_para, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)   
    #HET
    #main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    #ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
  
   
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1*readout_amp_2,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
        ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom
 
    
    
def parametric_coupling_time_domain_Raman_Q2(num_steps=101,ssm_ge_Q1 = -0.2,ssm_ge_Q2 = -0.2,pi_ge_Q1 =20,pi_ge_Q2 =20,ssm_para=0,spec_amp=.5,ROIF1=0,ROIF2=0,q=0,sweep_time=0):
    
    totlength = sweep_time + 4000
    sweeptime = sweep_time
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    readout_amp = 1 
    oscNum = 6
    phase_offset = mixer_offset
    detun = -0.03 
    ssm_coax = ROIF2 + detun
    amp_coax = 1#*3.317#0#1
    ssm_Q2 = ssm_ge_Q2 + detun + 0.001 # +0.001 for amp_Q2=0.5 # excitation -detun; decay + detun
    amp_Q2 = 0.5#*3.98#0#0.5 #0.5 #0
    the_seq = Sequence(file_length, num_steps)
    coaxtime = sweeptime#-5000 #500*4
    pulse_len_add = 0
    readout_dur = ro_pulse_dur
    ge_amp = ge_amp_setting
     #qubit 2
    ifsideband = 1           
         # off-resonance drive for RO (of Q2)
    coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
#        coax_drive.phase = 90
#        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q2
    q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1
    
    
    sigma = 40# 200+0*300 #100
    width = int(sigma*5/2)

    sweeptime = coaxtime
    if ifsideband == 1:
        for k in np.arange(num_steps):
            if k >= 0:  # Note (WC): The first pulse with nonzero duration should larger than 2*width

            # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                pulseStart = file_length - readout_dur - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                pulseStart2 = file_length - readout_dur - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                tdata = np.arange(pulseLength)
                tdata2 = np.arange(pulseLength2)
                pdata1 = np.zeros(pulseLength)
                pdata2 = np.zeros(pulseLength2)

                # gaussian shape                    
                pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
                pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                
                # sine shape
#                    pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata[0:width]) 
#                    pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata[0:width])
                
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                pdata1[width:secondStop] = amp_Q2
                pdata2[width:secondStop2] = amp_coax
                
                pdata1[secondStop:pulseLength] = (amp_Q2 *
                      np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
                pdata2[secondStop2:pulseLength2] = (amp_coax *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                      
                
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                the_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                
                the_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
                
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()
        #pi_ge = Pulse(start=file_length-readout_dur-100, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        #ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
     #qubit 1 pi pulse
    pi_ge = Pulse(start=file_length-readout_dur-400, duration=-pi_ge_Q1, amplitude=ge_amp, ssm_freq=ssm_ge_Q1, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    parametric_12 = Pulse(start=file_length-readout_dur-10, duration=0, amplitude=spec_amp, ssm_freq=ssm_para, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=parametric_12)
    
    
    #Readout
    #amplitude = *0 + 0.035
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1*0+0.035*1.525,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    #amplitude = *0 + 0.035

    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp_2*0+0.035,ssm_freq=ROIF2,phase=-file_length*ROIF2*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        

    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
def parametric_coupling_Q1Ramsey_Raman_Q2(num_steps=101,ssm_ge_Q1 = -0.2,ssm_ge_Q2 = -0.2,pi_ge_Q1 =20,pi_ge_Q2 =20,ssm_para=0,spec_amp=.5,ROIF1=0,ROIF2=0,q=0,sweep_time=0):
    
    totlength = sweep_time + 4000
    sweeptime = sweep_time
    file_length = 10000 * (int(np.ceil(totlength/10000))+1)
    readout_amp = 1 
    oscNum = 0
    phase_offset = mixer_offset
    detun = -0.03 
    ssm_coax = ROIF2 + detun
    amp_coax = 1#*3.317#0#1
    ssm_Q2 = ssm_ge_Q2 + detun + 0.001 # +0.001 for amp_Q2=0.5 # excitation -detun; decay + detun
    amp_Q2 = 0.5#*3.98#0#0.5 #0.5 #0
    the_seq = Sequence(file_length, num_steps)
    coaxtime = sweeptime#-5000 #500*4
    pulse_len_add = 0
    readout_dur = ro_pulse_dur
    ge_amp = ge_amp_setting
     #qubit 2
    ifsideband = 1           
         # off-resonance drive for RO (of Q2)
    coax_drive = Pulse(start=file_length-readout_dur-pi_ge_Q1/2-10, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
#        coax_drive.phase = 90
#        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q2
    q2_off_resonance_drive = Pulse(start=file_length-readout_dur-pi_ge_Q1/2-10, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1
    
    
    sigma = 40# 200+0*300 #100
    width = int(sigma*5/2)

    sweeptime = coaxtime
    if ifsideband == 1:
        for k in np.arange(num_steps):
            if k >= 0:  # Note (WC): The first pulse with nonzero duration should larger than 2*width

            # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                pulseStart = file_length - readout_dur - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                pulseStart2 = file_length - readout_dur - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                tdata = np.arange(pulseLength)
                tdata2 = np.arange(pulseLength2)
                pdata1 = np.zeros(pulseLength)
                pdata2 = np.zeros(pulseLength2)

                # gaussian shape                    
                pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
                pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                
                # sine shape
#                    pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata[0:width]) 
#                    pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata[0:width])
                
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                pdata1[width:secondStop] = amp_Q2
                pdata2[width:secondStop2] = amp_coax
                
                pdata1[secondStop:pulseLength] = (amp_Q2 *
                      np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
                pdata2[secondStop2:pulseLength2] = (amp_coax *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                      
                
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                the_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                
                the_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
                
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()
        #pi_ge = Pulse(start=file_length-readout_dur-100, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        #ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
 
    #qubit 1 pi/2 pulse
    pi_ge = Pulse(start=file_length-readout_dur-pi_ge_Q1/2-100, duration=-pi_ge_Q1/2, amplitude=ge_amp, ssm_freq=ssm_ge_Q1, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)

    #qubit 1 second pi/2 pulse
    if oscNum !=0:
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_Q1/2, amplitude=ge_amp, ssm_freq=ssm_ge_Q1, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=oscNum*360, initial_pulse=pi_ge)
    elif oscNum == 0:
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_Q1/2, amplitude=ge_amp, ssm_freq=ssm_ge_Q1, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge)
    
    parametric_12 = Pulse(start=file_length-readout_dur-20-pi_ge_Q1/2, duration=0, amplitude=spec_amp, ssm_freq=ssm_para, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=parametric_12)
    
    
    #Readout
    #amplitude = *0 + 0.035
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp_1*0+0.035*1.525,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    #amplitude = *0 + 0.035

    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp_2*0+0.035,ssm_freq=ROIF2,phase=-file_length*ROIF2*360 ) 
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        

    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
def e_char(ge_coef=1,ssm_ef=-0.065): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200#20000# #300 #1000 #3000
    ## channels   
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
    pi_ge = pi_ge_time_setting
    ef_amp = ge_amp_setting
#    ssm_ef = ssm_ef_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    buffer = 50
    
    #first pi_ge pulse
    pi_ge_pulse = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=-pi_ge*ge_coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
    pi_ge_pulse.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
    
    #drive rabi e-f
    rabi_ef = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    rabi_ef.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    
    #second pi_ge-pulse    
    pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge_pulse)
    pi_ge_pulse.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
  
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom

def repeat_ramsey_ef(repnum = 20, ifsave = 0, ifrunsequence=1):
    
    sweeptime = 10000
    num_steps = 101
    if ifrunsequence:
        ramsey_ef(num_steps= num_steps, ifrotphase=0, sweeptime = sweeptime, ifload = 1)
    
    totprob = np.zeros([repnum, 101])
    
    for k in np.arange(repnum):
        sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                   sweepparam = np.linspace(0,sweeptime*1e-3,101), ave = 1, 
                   ifprint = 0, ifsave = 0, ifplot=0)
        
        plt.figure()
        plt.plot(sweepparam, prob, 'k.-')
        plt.xlabel('Time ($\\mu$s)')
        plt.ylabel('$P_e$')
        
        plt.title('# {} in {}'.format(k+1, repnum))
        plt.show()
        totprob[k,:] = prob
    
    
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob, 'repnum':repnum}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
        
    plt.figure()
    x = sweepparam
    y = np.arange(repnum)
    xd = np.append(x-(x[1]-x[0])/2,x[-1]+(x[1]-x[0])/2)
    yd = np.append(y-(y[1]-y[0])/2,y[-1]+(y[1]-y[0])/2)
    X, Y = np.meshgrid(xd, yd)

    print(np.shape(totprob))
    plt.pcolormesh(X, Y, totprob)
    
    plt.xlabel('Time ($\\mu$s)')
    plt.ylabel('data #')
    plt.show()
    
def rabi_ge_pi_2(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 300
    ## channels   
    pi_ge = pi_ge_time_setting
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp =0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    
    pi_ge = Pulse(start=file_length-readout_dur-2050, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=pi_ge)
    pi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-sweep_time,initial_pulse=pi_ge)
 
    
    rabi_ge = Pulse(start=file_length-readout_dur-2050, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    

    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur-1000,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur+1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom

def ramsey(num_steps = 101,t1_time = 100000,pi_echo_coef=0,osc=0,ssm_ge=-0.15,pi_ge=20,ROIF=0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 10000*(int(t1_time/10000)+1) #64000
    # totlength = t1_time + 4000
    # file_length = 10000 * (int(totlength/10000)+2)
    file_length = 16000
    #file_length = 100000#100000# 
#    num_steps = 51
    ge_amp = ge_amp_setting
    
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
    phase_offset = mixer_offset

#    pi_echo=0
    readout_amp = 1#0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    buffer =0#550#2050
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  
#    pi_ge=17
#    t1_time = 1000
#    ssm_ge = 0.3885
    if q == 0: #qubit 2
        readout_amp= readout_amp_2
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2-pi_echo_coef*pi_ge, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=2,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=t2_ge)
            
            
            
    if q == 1: #qubit 1
        readout_amp= readout_amp_1
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2-pi_echo_coef*pi_ge, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=4,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=t2_ge)
# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=50, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
def ramsey_quantum_effiency(num_steps = 101,ssm_ge=-0.15,pi_ge=20,ROIF=0,q=0, RO_ram_amp = 0.0): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 10000*(int(t1_time/10000)+1) #64000
    # totlength = t1_time + 4000
    #file_length = 10000 * (int(totlength/10000)+2)
    file_length = 16000#16000
    #file_length = 100000#100000# 
#    num_steps = 51
    ge_amp = ge_amp_setting
    
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
    phase_offset = mixer_offset

#    pi_echo=0
    readout_dur = 5000#2000#8000 #1000
    buffer =1000#500#550#2050
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  
#    pi_ge=17
#    t1_time = 1000
#    ssm_ge = 0.3885
    if q == 0: #qubit 2
        readout_amp= readout_amp_2
        t2_ge = Pulse(start=file_length- 2*readout_dur-buffer-pi_ge/2-100, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    
        #READOUT 
        main_pulse = Pulse(start = file_length- readout_dur-50- pi_ge/2 - buffer,duration= -readout_dur, amplitude= RO_ram_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
        ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)


        t2_ge = Pulse(start=file_length-readout_dur -50, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=360,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#         if osc !=0:
#             t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#             ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
# #            t2_ge.phase = 90+phase_offset
# #            ringupdown_seq.add_sweep(channel=2,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
#         elif osc == 0:
#             t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#             ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)
# #            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=t2_ge)
            
            
            
#     if q == 1: #qubit 1
#         readout_amp= readout_amp_1
#         t2_ge = Pulse(start=file_length- 2*readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#         ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=t2_ge)
# #        t2_ge.phase = 90+phase_offset
# #        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    
#         #READOUT 
#         main_pulse = Pulse(start = file_length- 2*readout_dur- pi_ge/2 - buffer,duration= readout_dur, amplitude= RO_ram_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
#         ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)



#         t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#         ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=360,initial_pulse=t2_ge)
# #        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#         if osc !=0:
#             t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#             ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
# #            t2_ge.phase = 90+phase_offset
# #            ringupdown_seq.add_sweep(channel=4,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
#         elif osc == 0:
#             t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#             ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=t2_ge)
# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=readout_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=50, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])    

def ramsey_modulation(num_steps = 101,t1_time = 100000,osc=0,ssm_ge=-0.15,pi_ge=20,ROIF=0,off=0,q=0,ifload = 0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = t1_time + 4000
    file_length = 10000 * (int(totlength/10000)+2) #80000

    ge_amp = ge_amp_setting

    phase_offset = mixer_offset

    readout_amp = 1#0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    buffer =10
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  

    if q == 0: #qubit 2
        readout_amp= readout_amp_2
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)


        

            
            
            
    if q == 1: #qubit 1
        readout_amp= readout_amp_1
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)


        mod_pulse = Pulse(start=file_length-readout_dur-buffer/2-pi_ge/2, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='width', start=0, stop=-t1_time,initial_pulse=mod_pulse)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ge)


# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

def ramsey_modulation_step(num_steps = 101,t1_time = 100000,osc=0,ssm_ge=-0.15,pi_ge=20,ROIF=0,off=0,q=0,ifload = 0):
    totlength = t1_time + 4000
    file_length = 10000 * (int(totlength/10000)+2) #80000

    ge_amp = ge_amp_setting

    phase_offset = mixer_offset

    readout_amp = 1#0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    buffer =20
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  

    if q == 0: #qubit 2
        readout_amp= readout_amp_2
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)

            
            
    if q == 1: #qubit 1
        readout_amp= readout_amp_1
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2-t1_time, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=t2_ge)


        mod_pulse = Pulse(start=file_length-readout_dur-buffer/2-pi_ge/2-t1_time, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='none',initial_pulse=mod_pulse)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ge)



# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

def ramsey_noise(num_steps = 101,t1_time = 100000,osc=0,ssm_ge=-0.15,pi_ge=20,ROIF=0,off=0,q=0,ifload = 0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = t1_time + 4000
    file_length = 10000 * (int(totlength/10000)+2) #80000

    ge_amp = ge_amp_setting

    phase_offset = mixer_offset

    readout_amp = 1#0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    buffer =10
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  

    if q == 0: #qubit 2
        readout_amp= readout_amp_2
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)


        

            
            
            
    if q == 1: #qubit 1
        readout_amp= readout_amp_1
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)


        noise_pulse = Pulse(start=file_length-readout_dur-buffer/2-pi_ge/2, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='width', start=0, stop=-t1_time,initial_pulse=noise_pulse)
        mod_pulse = Pulse(start=file_length-readout_dur-buffer/2-pi_ge/2, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=3, marker=2, sweep_name='width', start=0, stop=-t1_time,initial_pulse=mod_pulse)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ge)


# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

def ramsey_noise_only(num_steps = 101,t1_time = 100000,osc=0,ssm_ge=-0.15,pi_ge=20,ROIF=0,off=0,q=0,ifload = 0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = t1_time + 4000
    file_length = 10000 * (int(totlength/10000)+2) #80000

    ge_amp = ge_amp_setting

    phase_offset = mixer_offset

    readout_amp = 1#0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    buffer =10
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  

    if q == 0: #qubit 2
        readout_amp= readout_amp_2
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)


        

            
            
            
    if q == 1: #qubit 1
        readout_amp= readout_amp_1
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)


        noise_pulse = Pulse(start=file_length-readout_dur-buffer/2-pi_ge/2, amplitude=1, duration=0, phase=0) 
        ringupdown_seq.add_sweep(channel=1, marker=1, sweep_name='width', start=0, stop=-t1_time,initial_pulse=noise_pulse)
        
        
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)

        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ge)


# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF, phase=-file_length*ROIF*360 ) #amplitude= 1.3*readout_amp_1,ssm_freq=ROIF1, phase=-file_length*ROIF1*360
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    if ifload:
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
def ramsey_ef(num_steps = 101,t1_time = 100000,pi_echo_coef=0,osc=0,ssm_ge=-0.15,ssm_ef=-0.3,pi_ge=20,pi_ef = 30,ef_amp=1,q=0):
#    file_length = 10000*(int(t1_time/10000)+1) #64000
    totlength = t1_time + 4000
    file_length = 10000 * (int(totlength/10000)+2)
    #file_length = 100000#100000# 
#    num_steps = 51
    ge_amp = ge_amp_setting
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef

    readout_dur = ro_pulse_dur#8000 #1000
    buffer =0#550#2050
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class

    if q == 0: #qubit 2
        #first pi_ge pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ef/2-pi_echo_coef*pi_ef-pi_ef/2, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #first pi_ef/2 pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ef/2-pi_echo_coef*pi_ef, duration=-pi_ef/2, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #echo pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ef/2, duration=-pi_ef*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        
        #second pi_ef/2 pulse
        if osc !=0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ef/2, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=1, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=2,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
        elif osc == 0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ef/2, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=t2_ef)
            
        #second pi_ge pulse
#        t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
#        t2_ef.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
    if q == 1: #qubit 1
        #first pi_ge pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge-pi_ef/2-pi_echo_coef*pi_ef-pi_ef/2, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #first pi_ef/2 pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge-pi_ef/2-pi_echo_coef*pi_ef, duration=-pi_ef/2, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #echo pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge-pi_ef/2, duration=-pi_ef*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        
        #second pi_ef/2 pulse
        if osc !=0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge, duration=-pi_ef/2, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=4,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
        elif osc == 0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge, duration=-pi_ef/2, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=t2_ef)
            
        #second pi_ge pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4,  sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

#def rabi_ef(ssm_ge = -0.150,ssm_ef = -0.300): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 8000
#    num_steps = 101
#    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
#    
#    sweep_time = 100
#    ## channels   
#    
#    ge_amp = ge_amp_setting
#    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting
##    ssm_ef = ssm_ef_setting
#    readout_amp = 1 
#    
#    
#    pi_ge = Pulse(start=file_length-3050-pi_ge_time, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='start',start=0,stop=-sweep_time, initial_pulse=pi_ge)
#    pi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='start',start=0,stop=-sweep_time ,initial_pulse=pi_ge)
#    
#    rabi_ef = Pulse(start=file_length-3050-pi_ge_time, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#    rabi_ef.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#    
#    
#    pi_ge = Pulse(start=file_length-3050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
#    pi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
#    
#    
#    
#    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
#    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#
#    ## view output
#    if True:
#        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = ringupdown_seq.channel_list[1][0]
#        channel3_ch = ringupdown_seq.channel_list[2][0]
#        channel4_ch = ringupdown_seq.channel_list[3][0]
#        marker1 = ringupdown_seq.channel_list[0][2]
#        
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
#    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
###END geom

def check_f0_e1_transition_GaussianPulse(num_steps= 101, rotangleEF = np.pi,
                                         iffinalpiGE=1, sweeptime = 2000,
                                         ifload = 1, sigma = 40, ifsideband=1,
                                         sideAmp=2): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(totlength/10000)+1)
    
    if file_length > 80000:
        raise ValueError('File length too long. Make it less than 80000')
        file_length = 80000
        
                         
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
 
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    
    rotDurationEF = int(rotangleEF/np.pi * pi_ef_time)
    width = int(sigma*5/2)
    
    
    pi_ge = Pulse(start=file_length-3050-rotDurationEF-pi_ge_time-2*width - 1000, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    pi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    
    
    rot_ef = Pulse(start=file_length-3050-pi_ge_time-2*width - 1000, duration=-rotDurationEF, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=rot_ef)
    rot_ef.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=rot_ef)
    
    
    
    
    amp = sideAmp
    
    if ifsideband == 1:
        for k in np.arange(num_steps):
            pulseLength = int(2*width + sweeptime * k/(num_steps-1))
            pulseStart = file_length-3050-pi_ge_time-pulseLength - 1000
            tdata = np.arange(pulseLength)
            pdata = np.zeros(pulseLength)
            
            pdata[0:width] = amp * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
            
            secondStop = int(width+sweeptime * k/(num_steps-1))
            pdata[width:secondStop] = amp
            pdata[secondStop:pulseLength] = (amp *
                  np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
                  
            
            ringupdown_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] = pdata
            ringupdown_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] = pdata
        
        
#    print(np.shape(ringupdown_seq.channel_list))
    
    
    if iffinalpiGE == 1:
        pi_ge = Pulse(start=file_length-3050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
        pi_ge.phase = 90
        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    
    
    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
        

def meas_f0toe1_population(ave=5, rotangleEF = 1*np.pi, ifsideband=1, 
                           ifgaussian=1, ifsave=1, sideAmp=2):
    
    num_steps = 101
    sweeptime = 2000
    
    if ifgaussian:
        check_f0_e1_transition_GaussianPulse(num_steps= num_steps, rotangleEF = rotangleEF,
                                             iffinalpiGE=1, sweeptime = sweeptime,
                                             ifload = 1, sigma = 40, ifsideband=ifsideband,
                                             sideAmp=sideAmp)
        
    sweepparam, prob0 = readout_python(num_patterns = num_steps, num_records = 200, 
                   sweepparam = np.linspace(0,sweeptime,num_steps), ave = ave, 
                   ifprint = 0, ifsave = 0, ifplot=1)
    
    
    if ifgaussian:
        check_f0_e1_transition_GaussianPulse(num_steps= num_steps, rotangleEF = rotangleEF,
                                             iffinalpiGE=0, sweeptime = sweeptime,
                                             ifload = 1, sigma = 40, ifsideband=ifsideband,
                                             sideAmp=sideAmp)
        
    sweepparam, prob1 = readout_python(num_patterns = num_steps, num_records = 200, 
                   sweepparam = np.linspace(0,sweeptime,num_steps), ave = ave, 
                   ifprint = 0, ifsave = 0, ifplot=1)
    
    prob_g = 1 - prob1
    prob_e = 1 - prob0
    prob_f = 1 - prob_g - prob_e
    
    plt.figure()
    plt.plot(sweepparam, prob_g, 'r')
    plt.plot(sweepparam, prob_e, 'b')
    plt.plot(sweepparam, prob_f, 'g')
    
    plt.xlabel('Time (us)')
    plt.ylabel('Population')
    plt.legend(['$P_g$','$P_e$','$P_f$'], loc='lower left',bbox_to_anchor=(1,0))
    
    if ifsave:
        data = {'sweepparam':sweepparam, 'prob_g':prob_g,
                'prob_e':prob_e, 'prob_f':prob_f}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
    


def check_ge_dephasing_amp_coherentReadout(num_steps= 101, ifload = 1,
                                         ifplotsequence=1, amp = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 8000
    
    if file_length > 80000:
        raise ValueError('File length too long. Make it less than 80000')
        file_length = 80000
        
                         
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
 
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1
    oscNum = 6
    
    
    
    pi2_ge = Pulse(start=file_length-6050-pi2_ge_time, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi2_ge)
    
    main_pulse = Pulse(start = file_length-6050-pi2_ge_time,duration = 2000, amplitude= amp )
    ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    
    pi2_ge = Pulse(start=file_length-3050, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='phase',start=0,stop=360, initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='phase',start=0,stop=360, initial_pulse=pi2_ge)
    
    
    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
 

def check_ge_dephasing_amp_ThermalReadout(num_steps= 101, ifload = 1,
                                         ifplotsequence=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 8000
    
    if file_length > 80000:
        raise ValueError('File length too long. Make it less than 80000')
        file_length = 80000
        
                         
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
 
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1
    oscNum = 6
    
    thermal_amp = 2
    
    
    
    pi2_ge = Pulse(start=file_length-6050-pi2_ge_time, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi2_ge)
    
    main_pulse = Pulse(start = file_length-6050-pi2_ge_time,duration = 2000, amplitude= thermal_amp )
    ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    
    pi2_ge = Pulse(start=file_length-3050, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='phase',start=0,stop=360, initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='phase',start=0,stop=360, initial_pulse=pi2_ge)
    
    
    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
 


def check_f0_e1_transition_dephasing(num_steps= 101, rotangleEF = np.pi,
                                         iffinalpiEF=1, sweeptime = 2000,
                                         ifload = 1, sigma = 40, ifsideband=1,
                                         ifplotsequence=1, sideAmp=1.9): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(totlength/10000)+1)
    
    if file_length > 80000:
        raise ValueError('File length too long. Make it less than 80000')
        file_length = 80000
        
                         
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
 
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    
    rotDurationEF = int(rotangleEF/np.pi * pi_ef_time)
    width = int(sigma*5/2)
    
    
    pi2_ge = Pulse(start=file_length-3050-pi2_ge_time-3000-rotDurationEF, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi2_ge)
    
    
    
    rot_ef = Pulse(start=file_length-3050-pi2_ge_time-3000, duration=-rotDurationEF, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=rot_ef)
    rot_ef.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=rot_ef)
    
    
    
    
    amp = sideAmp
    
    if ifsideband == 1:
        for k in np.arange(num_steps):
            pulseLength = int(2*width + sweeptime)
            pulseStart = file_length-3050-pi2_ge_time-3000
            tdata = np.arange(pulseLength)
            pdata = np.zeros(pulseLength)
            
            pdata[0:width] = amp * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
            
            secondStop = int(width+sweeptime)
            pdata[width:secondStop] = amp
            pdata[secondStop:pulseLength] = (amp *
                  np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
                  
            
            ringupdown_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] = pdata
            ringupdown_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] = pdata
        
        
#    print(np.shape(ringupdown_seq.channel_list))
    
    
    if iffinalpiEF == 1:
        pi_ef = Pulse(start=file_length-3050-pi2_ge_time, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ef)
        pi_ef.phase = 90
        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ef)
    
    pi2_ge = Pulse(start=file_length-3050, duration=-pi2_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='phase',start=0,stop=360, initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='phase',start=0,stop=360, initial_pulse=pi2_ge)
    
    
    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
 
#    global seqInfo
#    sweeplabel = 'Phase'
#    seqInfo = {'num_steps':num_steps, 'sweepparam':np.linspace(0,2*np.pi,num_steps),
#                    'sweeplabel':sweeplabel}
#    


    
def check_f0_e1_transition_population(num_steps= 101, GErotangle = np.pi/2,
                                         iffinalpiGE=1, sweeptime = 2000, rotTime = 100,
                                         ifload = 1, sigma = 40, ifsideband=1,
                                         ifplotsequence=1, sideAmp=1.9): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (int(totlength/10000)+1)
    
    if file_length > 80000:
        raise ValueError('File length too long. Make it less than 80000')
        file_length = 80000
        
                         
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
 
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    oscNum = 6
    
    width = int(sigma*5/2)
    
    
    pi2_ge = Pulse(start=file_length-3050-pi_ge_time-3000, duration=-pi2_ge_time*int(GErotangle*2/np.pi), amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rotTime, initial_pulse=pi2_ge)
    pi2_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-rotTime, initial_pulse=pi2_ge)
    
    
    
    rot_ef = Pulse(start=file_length-3050-pi_ge_time-3000, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rotTime, initial_pulse=rot_ef)
    rot_ef.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-rotTime, initial_pulse=rot_ef)
    
    
    
    
    amp = sideAmp
    
    if ifsideband == 1:
        for k in np.arange(num_steps):
            pulseLength = int(2*width + sweeptime)
            pulseStart = file_length-3050-pi2_ge_time-3000
            tdata = np.arange(pulseLength)
            pdata = np.zeros(pulseLength)
            
            pdata[0:width] = amp * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
            
            secondStop = int(width+sweeptime)
            pdata[width:secondStop] = amp
            pdata[secondStop:pulseLength] = (amp *
                  np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
                  
            
            ringupdown_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] = pdata
            ringupdown_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] = pdata
        
        
#    print(np.shape(ringupdown_seq.channel_list))
    
    
    if iffinalpiGE == 1:
        pi_ge = Pulse(start=file_length-3050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
        pi_ge.phase = 90
        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    
    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
 
#    global seqInfo
#    sweeplabel = 'Phase'
#    seqInfo = {'num_steps':num_steps, 'sweepparam':np.linspace(0,2*np.pi,num_steps),
#                    'sweeplabel':sweeplabel}
#    

def angleSweep_f0toe1_dephasing(ave=5, singleave=5, transtime = 0, sigma=0, angle = np.linspace(0,np.pi,10),
                                ifsave=1, ifprint=0, sideAmp=1.9):
    
    num_steps = 101
    
    pi_ef_time = pi_ef_time_setting

    
    
    actualAngle = np.zeros(np.size(angle))
    totprobGF = np.zeros([ave, np.size(angle), num_steps])
    totprobGE = np.zeros([ave, np.size(angle), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for l in np.arange(ave):
        for k in np.arange(np.size(angle)):
            actualAngle[k] = np.pi/pi_ef_time*int(angle[k]/np.pi * pi_ef_time)
            
            check_f0_e1_transition_dephasing(num_steps= num_steps, rotangleEF = angle[k],
                                             iffinalpiEF=1, sweeptime = transtime,
                                             ifload = 1, sigma = sigma, ifsideband=1,
                                             ifplotsequence=1, sideAmp=sideAmp)
            
            sweepparam, prob0 = readout_python(num_patterns = num_steps, num_records = 200, 
                           sweepparam = np.linspace(0,2*np.pi,num_steps), ave = singleave, 
                           ifprint = 0, ifsave = 0, ifplot=0)
            
    
            check_f0_e1_transition_dephasing(num_steps= num_steps, rotangleEF = angle[k],
                                             iffinalpiEF=0, sweeptime = transtime,
                                             ifload = 1, sigma = sigma, ifsideband=1,
                                             ifplotsequence=0, sideAmp=sideAmp)
            
            sweepparam, prob1 = readout_python(num_patterns = num_steps, num_records = 200, 
                           sweepparam = np.linspace(0,2*np.pi,num_steps), ave = singleave, 
                           ifprint = 0, ifsave = 0, ifplot=0)
            
            totprobGF[l, k, :] = prob0
            totprobGE[l, k, :] = prob1
            
            plt.figure()
            plt.plot(sweepparam, prob0, 'r')
            plt.plot(sweepparam, prob1, 'b')
            plt.xlabel('Phase')
            plt.ylabel('$P_e$')
            plt.title('# {} in {}'.format(k+1+l*np.size(angle), np.size(angle)*ave))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprobGF':totprobGF,
                'totprobGE':totprobGE, 'actualAngle':actualAngle, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)  
        
        
def timeSweep_f0toe1_dephasing(ave=1, sigma=40, sweeptime = np.linspace(0,2000,51),
                                ifsave=1, ifprint=0):
    
    num_steps = 101

    
    
    totprobGF = np.zeros([np.size(sweeptime), num_steps])
    totprobGE = np.zeros([np.size(sweeptime), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    
    for k in np.arange(np.size(sweeptime)):
        
        check_f0_e1_transition_dephasing(num_steps= num_steps, rotangleEF = np.pi,
                                         iffinalpiEF=1, sweeptime = sweeptime[k],
                                         ifload = 1, sigma = sigma, ifsideband=1,
                                         ifplotsequence=1)
        
        sweepparam, prob0 = readout_python(num_patterns = num_steps, num_records = 200, 
                       sweepparam = np.linspace(0,2*np.pi,num_steps), ave = ave, 
                       ifprint = 0, ifsave = 0, ifplot=0)
        

        check_f0_e1_transition_dephasing(num_steps= num_steps, rotangleEF = np.pi,
                                         iffinalpiEF=0, sweeptime = sweeptime[k],
                                         ifload = 1, sigma = sigma, ifsideband=1,
                                         ifplotsequence=0)
        
        sweepparam, prob1 = readout_python(num_patterns = num_steps, num_records = 200, 
                       sweepparam = np.linspace(0,2*np.pi,num_steps), ave = ave, 
                       ifprint = 0, ifsave = 0, ifplot=0)
        
        totprobGF[k,:] = prob0
        totprobGE[k,:] = prob1
        
        plt.figure()
        plt.plot(sweepparam, prob0, 'r')
        plt.plot(sweepparam, prob1, 'b')
        plt.xlabel('Phase')
        plt.ylabel('$P_e$')
        plt.title('# {} in {}: {} ns'.format(k+1, np.size(sweeptime), sweeptime[k]))
        plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprobGF':totprobGF,
                'totprobGE':totprobGE, 'sweeptime':sweeptime}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 


def ampSweep_ge_dephasingAmp(ave=1, amp = np.arange(0.05,1.5,0.1),
                                ifsave=1, ifprint=0):
    
    num_steps = 101
       
    
    totprob = np.zeros([np.size(amp), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    
    for k in np.arange(np.size(amp)):
                
        check_ge_dephasing_amp_coherentReadout(num_steps= 101, ifload = 1,
                                         ifplotsequence=0, amp = amp[k])
        
        sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                       sweepparam = np.linspace(0,2*np.pi,num_steps), ave = ave, 
                       ifprint = 0, ifsave = 0, ifplot=0)
        

        
        totprob[k,:] = prob
        
        plt.figure()
        plt.plot(sweepparam, prob, 'k.-')
        plt.xlabel('Phase')
        plt.ylabel('$P_e$')
        plt.title('# {} in {}: {} amp'.format(k+1, np.size(amp), amp[k]))
        plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'amp':amp}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)  
        
def ampSweep_ge_dephasingAmp_withave(ave=5, amp = np.arange(0.05,1.41,0.05),
                            ifsave=1, ifprint=0):
    
    num_steps = 101
       
    
    totprob = np.zeros([ave, np.size(amp), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(ave):
        for k in np.arange(np.size(amp)):
                    
            check_ge_dephasing_amp_coherentReadout(num_steps= 101, ifload = 1,
                                             ifplotsequence=0, amp = amp[k])
            
            sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                           sweepparam = np.linspace(0,2*np.pi,num_steps), ave = 5, 
                           ifprint = 0, ifsave = 0, ifplot=0)
            
    
            
            totprob[m, k, :] = prob
            
            plt.figure()
            plt.plot(sweepparam, prob, 'k.-')
            plt.xlabel('Phase')
            plt.ylabel('$P_e$')
            plt.title('# {} in {}: {} amp'.format(k+1+m*np.size(amp), 
                      np.size(amp)*ave, amp[k]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'amp':amp, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        
def PowerSweep_ge_dephasingAmp_withave_thermalReadout(ave=5, att = np.arange(40,0-0.1,-2),
                            ifsave=1, ifprint=0):
    
    num_steps = 101
       
    
    totprob = np.zeros([ave, np.size(att), num_steps])
    
    
    check_ge_dephasing_amp_ThermalReadout(num_steps= num_steps, ifload = 1,
                                         ifplotsequence=1)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(ave):
        for k in np.arange(np.size(att)):
                    
            vc.set_attenuation(int(att[k]))
            
            sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                           sweepparam = np.linspace(0,2*np.pi,num_steps), ave = 5, 
                           ifprint = 0, ifsave = 0, ifplot=0)
            
    
            
            totprob[m, k, :] = prob
            
            plt.figure()
            plt.plot(sweepparam, prob, 'k.-')
            plt.xlabel('Phase')
            plt.ylabel('$P_e$')
            plt.title('# {} in {}: {} dB att'.format(k+1+m*np.size(att), 
                      np.size(att)*ave, att[k]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'att':att, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        
        
def PowerSweep_ge_ReadoutSNR_withave_thermalReadout(ave=5, att = np.arange(40,0-0.1,-2),
                            ifsave=1, ifprint=0):
    
    num_steps = 101
       
    
    totvoltage = np.zeros([ave, np.size(att), 49])
    tothist_g = np.zeros([ave, np.size(att), 49])
    tothist_e = np.zeros([ave, np.size(att), 49])
    
    
    check_readout_multiplestate(measure = 'thermal', num_steps=5, measTime = 2000)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(ave):
        for k in np.arange(np.size(att)):
                    
            vc.set_attenuation(int(att[k]))
            
            if att[k] > 15:
                _, voltage, hist_g, hist_e, _ = xr.compare_ge_readout_FT_simultaneously_nosequence(
                        num_steps=5, fJPA = 5.6205, numrec=2000, 
                        ifplot=0, ifplotSNR=0,ifsave=0)
            elif att[k] <= 15:
                _, voltage, hist_g, hist_e, _ = xr.compare_ge_readout_FT_simultaneously_nosequence_weight(
                        num_steps=5, fJPA = 5.6205, numrec=2000, 
                        ifplot=0, ifplotSNR=0, ifsave=0)
            
    
            
            totvoltage[m, k, :] = voltage
            tothist_g[m, k, :] = hist_g
            tothist_e[m, k, :] = hist_e
            
            plt.figure()
            plt.plot(voltage, hist_g, 'b.-')
            plt.plot(voltage, hist_e, 'r.-')
            plt.xlabel('voltage')
            plt.ylabel('Count')
            plt.title('# {} in {}: {} dB att'.format(k+1+m*np.size(att), 
                      np.size(att)*ave, att[k]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'totvoltage':totvoltage, 'tothist_g':tothist_g,
                'tothist_e':tothist_e,
                'att':att, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        
        
def PowerSweep_ge_ReadoutSNR_withave_thermalReadout_allWeightMethod(
        ave=5, att = np.arange(30,0-0.1,-1), cutAtt=21, ifsave=1, ifprint=0, numrec=4000):
    
    num_steps = 101
       
    
    totvoltage = np.zeros([ave, np.size(att), 49])
    tothist_g = np.zeros([ave, np.size(att), 49])
    tothist_e = np.zeros([ave, np.size(att), 49])
    
    
    check_readout_multiplestate(measure = 'thermal', num_steps=5, measTime = 2000)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(ave):
        for k in np.arange(np.size(att)):
                    
            vc.set_attenuation(int(att[k]))
            
            if att[k] > cutAtt:
                _, voltage, hist_g, hist_e, _, _, _ = xr.compare_ge_readout_FT_simultaneously_nosequence_weight_fixedWeight(
                        num_steps=5, fJPA = 5.6205,numrec=numrec, 
                        ifplot=0, ifplotSNR=0,ifsave=0,weightType=2,width=1)
            elif att[k] <= cutAtt:
                _, voltage, hist_g, hist_e, _ = xr.compare_ge_readout_FT_simultaneously_nosequence_weight(
                        num_steps=5, fJPA = 5.6205, numrec=numrec, 
                        ifplot=0, ifplotSNR=0, ifsave=0)
            
    
            
            totvoltage[m, k, :] = voltage
            tothist_g[m, k, :] = hist_g
            tothist_e[m, k, :] = hist_e
            
            plt.figure()
            plt.plot(voltage, hist_g, 'b.-')
            plt.plot(voltage, hist_e, 'r.-')
            plt.xlabel('voltage')
            plt.ylabel('Count')
            plt.title('# {} in {}: {} dB att'.format(k+1+m*np.size(att), 
                      np.size(att)*ave, att[k]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'totvoltage':totvoltage, 'tothist_g':tothist_g,
                'tothist_e':tothist_e,
                'att':att, 'ave':ave, 'cutAtt':cutAtt}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 


def single_readout_monoThermal(n = 1, totnum = 5000, ifplot=1, ifprint=0, 
                              phase = np.arange(0,2*np.pi,np.pi/4), ifsave=0):
    
    num_steps = np.size(phase)
    
    ch_a_g = np.zeros(0)
    ch_a_e = np.zeros(0)
    ch_b_g = np.zeros(0)
    ch_b_e = np.zeros(0)
    
    
    para, fitfunc = dap.plot_alpha_vs_amp(fileno=31, ifplot=0)
    step = 0.4
    alpha = np.arange(0,3.5-step/2,step) + step/2
#    alpha = np.array([3]) 
    amp = np.zeros(np.size(alpha))
    
    for k in np.arange(np.size(amp)):
        result = fsolve(dap.getampfromalpha, [0.5], args=(alpha[k], para))
        amp[k] = result[0]
    
    
    p = 1/(np.pi * n) * np.exp(-alpha**2/n) * ( 2 * 
                  np.pi * alpha * step)
    p = p/np.sum(p)
    
    totrepnum = np.rint(p * totnum * 1/10) * 10
    print(totrepnum)
    print(amp)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
    for k in np.arange(np.size(amp)):
        repnum = int(totrepnum[k])
     
        if repnum>0:
            check_readout_multiplestate_differentPhase(
                                    measTime = 2000, c_amp = amp[k], 
                                    phase=phase,
                                    ifplotsequence=0)
            
            [current_ch_a_g, current_ch_a_e, current_ch_b_g, 
             current_ch_b_e] = xr.compare_ge_readout_simultaneously_noSequence(
                    num_steps=num_steps, numrec=repnum,ifplotSNR=0,ifplotamp=0,
                    ifreturn=3, ifplot=0)
            
            ch_a_g = np.append(ch_a_g, current_ch_a_g)
            ch_a_e = np.append(ch_a_e, current_ch_a_e)
            ch_b_g = np.append(ch_b_g, current_ch_b_g)
            ch_b_e = np.append(ch_b_e, current_ch_b_e)
    
    ch_a_g_reshape = np.reshape(ch_a_g,[int(np.size(ch_a_g)/(num_steps)),num_steps])
    ch_b_g_reshape = np.reshape(ch_b_g,[int(np.size(ch_b_g)/(num_steps)),num_steps])
    ch_a_g_correct = np.zeros([int(np.size(ch_a_g)/(num_steps)),num_steps])
    ch_b_g_correct = np.zeros([int(np.size(ch_b_g)/(num_steps)),num_steps])
    
    ch_a_e_reshape = np.reshape(ch_a_e,[int(np.size(ch_a_e)/(num_steps)),num_steps])
    ch_b_e_reshape = np.reshape(ch_b_e,[int(np.size(ch_b_e)/(num_steps)),num_steps])
    ch_a_e_correct = np.zeros([int(np.size(ch_a_e)/(num_steps)),num_steps])
    ch_b_e_correct = np.zeros([int(np.size(ch_b_e)/(num_steps)),num_steps])
    
    for m in np.arange(np.size(phase)):
        rec_rot = xr.rotate_iq(phase[m]*180/np.pi, [ch_a_g_reshape[:,m], ch_b_g_reshape[:,m]]) 
        ch_a_g_correct[:,m] = rec_rot[0] 
        ch_b_g_correct[:,m] = rec_rot[1]
        
        rec_rot = xr.rotate_iq(phase[m]*180/np.pi, [ch_a_e_reshape[:,m], ch_b_e_reshape[:,m]]) 
        ch_a_e_correct[:,m] = rec_rot[0] 
        ch_b_e_correct[:,m] = rec_rot[1]
            
    ch_a_g_correct  = ch_a_g_correct.flatten()
    ch_b_g_correct  = ch_b_g_correct.flatten()
    
    ch_a_e_correct  = ch_a_e_correct.flatten()
    ch_b_e_correct  = ch_b_e_correct.flatten()
    
    
    
    if ifprint == 0:
        sys.stdout = _original_stdout
    
    
    lim = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e),
                                 np.abs(ch_b_g), np.abs(ch_b_e))))*1.2
    binnum = 50
    
    
    
    
    if ifplot:
        plt.figure()
        ax_scatter = plt.axes()
        ax_scatter.scatter(ch_a_g, ch_b_g, s=3, c='b', marker='.', alpha=0.2)
        ax_scatter.scatter(ch_a_e, ch_b_e, s=3, c='r', marker='.', alpha=0.2)
    
        ax_scatter.set_xlim([-lim, lim])
        ax_scatter.set_ylim([-lim, lim])
        
        
        ax_scatter.plot(np.mean(ch_a_g),np.mean(ch_b_g), 'bo')
        ax_scatter.plot(np.mean(ch_a_e),np.mean(ch_b_e), 'kx')

        plt.xlabel('ch1')
        plt.ylabel('ch2')
        
        plt.gca().set_aspect('equal')
        plt.show()
        
        plt.figure()
        ax_scatter = plt.axes()
        ax_scatter.scatter(ch_a_g_correct, ch_b_g_correct, s=3, c='b', marker='.', alpha=0.2)
        ax_scatter.scatter(ch_a_e, ch_b_e, s=3, c='r', marker='.', alpha=0.2)
    
        ax_scatter.set_xlim([-lim, lim])
        ax_scatter.set_ylim([-lim, lim])
        
        
        ax_scatter.plot(np.mean(ch_a_g),np.mean(ch_b_g), 'bo')
        ax_scatter.plot(np.mean(ch_a_e),np.mean(ch_b_e), 'kx')

        plt.xlabel('ch1')
        plt.ylabel('ch2')
        
        plt.gca().set_aspect('equal')
        plt.show()
        
        
    bins_y = np.linspace(-lim, lim, binnum)
    hist_g, bins = np.histogram(ch_b_g_correct, bins=bins_y)
    hist_e, bins = np.histogram(ch_b_e_correct, bins=bins_y)
    binwidth = bins[1]-bins[0]
    voltage = bins_y[0:-1]+binwidth/2
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplot:
        index = np.argmax(ydata_g)
        startpoint = [np.max(ydata_g), xdata[index], np.ptp(xdata)/6]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        index = np.argmax(ydata_e)
        startpoint = [np.max(ydata_e), xdata[index], np.ptp(xdata)/6]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        print(startpoint)
        print(para_e)
        plt.title('SNR is {}'.format(SNR))
        
        plt.xlabel('Ch2 reading')
        plt.ylabel('Count')
        plt.show()    
    
    bins_amp = np.linspace(-lim/3, lim, binnum)
    hist_g_amp, bins = np.histogram(np.sqrt(ch_b_g**2 + ch_a_g**2), bins=bins_amp)
    hist_e_amp, bins = np.histogram(np.sqrt(ch_b_e**2 + ch_a_e**2), bins=bins_amp)
    
    binwidth = bins[1]-bins[0]
    voltage_amp = bins_amp[0:-1]+binwidth/2
    xdata = voltage_amp    
    
    if ifplot:
        
        plt.figure()
#        plt.plot(xdata, hist_g_amp, 'b.-')
#        plt.plot(xdata, hist_e_amp, 'r.-')
        plt.xlabel('Amplitude')
        plt.ylabel('Count')
        
        ydata_g = hist_g_amp
        ydata_e = hist_e_amp
        
        index = np.argmax(ydata_g)
        startpoint = [np.max(ydata_g), xdata[index], np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        index = np.argmax(ydata_e)
        startpoint = [np.max(ydata_e), xdata[index], np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        plt.title('SNR is {}'.format(SNR))
        
        plt.show()
    
    if ifsave:
        data = {'ch_a_g_correct':ch_a_g_correct, 'ch_a_e_correct':ch_a_e_correct,
                'ch_b_g_correct':ch_b_g_correct,'ch_b_e_correct':ch_b_e_correct,
                'ch_a_g':ch_a_g, 'ch_a_e':ch_a_e,
                'ch_b_g':ch_b_g, 'ch_b_e':ch_b_e,
                'n':n, 'phase':phase}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
    
    return voltage_amp, hist_g_amp, hist_e_amp, voltage, hist_g, hist_e


def single_readout_monoThermal_FTmethod(n = 1, totnum = 5000, ifplot=1, ifprint=0, 
                              phase = np.arange(0,2*np.pi,np.pi/4), ifsave=0, 
                              fJPA = 5.6205):

    
    num_steps = np.size(phase)
    
    sig_gstate = np.zeros(0)
    sig_estate = np.zeros(0)
    
    
    para, fitfunc = dap.plot_alpha_vs_amp(fileno=31, ifplot=0)
    step = 0.4
    alpha = np.arange(0,3.5-step/2,step) + step/2
#    alpha = np.array([3]) 
    amp = np.zeros(np.size(alpha))
    
    for k in np.arange(np.size(amp)):
        result = fsolve(dap.getampfromalpha, [0.5], args=(alpha[k], para))
        amp[k] = result[0]
    
    
    p = 1/(np.pi * n) * np.exp(-alpha**2/n) * ( 2 * 
                  np.pi * alpha * step)
    p = p/np.sum(p)
    
    totrepnum = np.rint(p * totnum * 1/10) * 10
    print(totrepnum)
    print(amp)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
    for k in np.arange(np.size(amp)):
        repnum = int(totrepnum[k])
     
        if repnum>0:
            check_readout_multiplestate_differentPhase(
                                    measTime = 2000, c_amp = amp[k], 
                                    phase=phase,
                                    ifplotsequence=0)
            
            [_, _, _, _, _, current_sig_gstate, current_sig_estate
             ] = xr.compare_ge_readout_FT_simultaneously_nosequence_weight_fixedWeight(
             num_steps=num_steps, fJPA = fJPA, numrec=repnum, 
             ifplot=0, ifplotSNR=0,ifsave=0,weightType=1)
            
            sig_gstate = np.append(sig_gstate, current_sig_gstate)
            sig_estate = np.append(sig_estate, current_sig_estate)

    
    if ifprint == 0:
        sys.stdout = _original_stdout
    
    
    lim = np.max(np.abs(np.append(sig_estate, sig_gstate)))*1.2   
    binnum = 50
    
    
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    
    voltage = bins_y[0:-1]+binwidth/2
   
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplot:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
    
    if ifsave:
        data = {'hist_g':hist_g,'hist_e':hist_e,
                'voltage':voltage, 
                'n':n, 'phase':phase}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
    
    return voltage, hist_g, hist_e
        

def single_readout_monoThermal_FTmethod_getSpectrum(n = 1, totnum = 5000, ifplot=1, ifprint=0, 
                              phase = np.arange(0,2*np.pi,np.pi/4), ifsave=0, 
                              fJPA = 5.6205):

    
    num_steps = np.size(phase)
    
    sig_gstate = np.zeros(0)
    sig_estate = np.zeros(0)
    
    
    
    para, fitfunc = dap.plot_alpha_vs_amp(fileno=31, ifplot=0)
    step = 0.4
    alpha = np.arange(0,3.5-step/2,step) + step/2
#    alpha = np.array([3]) 
    amp = np.zeros(np.size(alpha))
    
    totspectrum_gstate = np.zeros([np.size(amp), 2999])
    totspectrum_estate = np.zeros([np.size(amp), 2999])
    totspectrum_gstate_bg = np.zeros([np.size(amp), 2999])
    totspectrum_estate_bg = np.zeros([np.size(amp), 2999])
    
    for k in np.arange(np.size(amp)):
        result = fsolve(dap.getampfromalpha, [0.5], args=(alpha[k], para))
        amp[k] = result[0]
    
    
    p = 1/(np.pi * n) * np.exp(-alpha**2/n) * ( 2 * 
                  np.pi * alpha * step)
    p = p/np.sum(p)
    
    totrepnum = np.rint(p * totnum * 1/10) * 10
    print(totrepnum)
    print(amp)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
    for k in np.arange(np.size(amp)):
        repnum = int(totrepnum[k])
     
        if repnum>0:
            check_readout_multiplestate_differentPhase(
                                    measTime = 2000, c_amp = amp[k], 
                                    phase=phase,
                                    ifplotsequence=0)
            
            [freq, mean_spec_gstate, mean_spec_estate, mean_spec_gstate_bg, 
             mean_spec_estate_bg, current_sig_gstate, current_sig_estate
             ] = xr.compare_ge_readout_FT_simultaneously_nosequence_weight_fixedWeight(
             num_steps=num_steps, fJPA = fJPA, numrec=repnum, 
             ifplot=0, ifplotSNR=0,ifsave=0,weightType=1,ifreturn=2)
            
            sig_gstate = np.append(sig_gstate, current_sig_gstate)
            sig_estate = np.append(sig_estate, current_sig_estate)
            totspectrum_gstate[k,:] = mean_spec_gstate
            totspectrum_estate[k,:] = mean_spec_estate
            totspectrum_gstate_bg[k,:] = mean_spec_gstate_bg
            totspectrum_estate_bg[k,:] = mean_spec_estate_bg
            
            if ifplot:
                plt.figure()
                plt.plot(1e3*freq, mean_spec_gstate, 'b-')
                plt.plot(1e3*freq, mean_spec_estate, 'r-')
                plt.xlabel('Freq (MHz)')
                plt.xlim([-20,20])
                plt.title('#{} in {}: {} amp'.format(k+1, np.size(amp), amp[k]))
                plt.show()

    
    if ifprint == 0:
        sys.stdout = _original_stdout
    
    
    lim = np.max(np.abs(np.append(sig_estate, sig_gstate)))*1.2   
    binnum = 50
    
    
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    
    voltage = bins_y[0:-1]+binwidth/2
   
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplot:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
    
    if ifsave:
        data = {'totspectrum_gstate':totspectrum_gstate,
                'totspectrum_estate':totspectrum_estate,
                'totspectrum_gstate_bg':totspectrum_gstate_bg,
                'totspectrum_estate_bg':totspectrum_estate_bg,
                'sig_gstate':sig_gstate,
                'sig_estate':sig_estate,
                'amp':amp, 'totrepnum':totrepnum, 'freq':freq,
                'n':n, 'phase':phase}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
    
    return voltage, hist_g, hist_e
    
def AmpSweep_ge_ReadoutSNR_withave_monoThermalReadout(ave=5,photonnum = np.arange(0.2,8+0.3,0.5),
                            ifsave=1):
    
    totvoltage_amp = np.zeros([ave, np.size(photonnum), 49])
    tothist_g_amp = np.zeros([ave, np.size(photonnum), 49])
    tothist_e_amp = np.zeros([ave, np.size(photonnum), 49])
    
    totvoltage_correct = np.zeros([ave, np.size(photonnum), 49])
    tothist_g_correct = np.zeros([ave, np.size(photonnum), 49])
    tothist_e_correct = np.zeros([ave, np.size(photonnum), 49])
    
    
    
    
    

    for m in np.arange(ave):
        for k in np.arange(np.size(photonnum)):
                    
            
            [voltage_amp, hist_g_amp, hist_e_amp, voltage, hist_g, 
             hist_e] = single_readout_monoThermal(n = photonnum[k], totnum = 5000, ifplot=0, ifprint=0, 
                              phase = np.arange(0,2*np.pi,np.pi/4))
            
    
            
            totvoltage_amp[m, k, :] = voltage_amp
            tothist_g_amp[m, k, :] = hist_g_amp
            tothist_e_amp[m, k, :] = hist_e_amp
            
            totvoltage_correct[m, k, :] = voltage
            tothist_g_correct[m, k, :] = hist_g
            tothist_e_correct[m, k, :] = hist_e
            
            
            plt.figure()
            xdata = voltage
            ydata_g = hist_g
            ydata_e = hist_e
    
            index = np.argmax(ydata_g)
            startpoint = [np.max(ydata_g), xdata[index], np.ptp(xdata)/6]
            fitfunc = Gaussian;
            para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
            perr_g = np.sqrt(np.diag(pcov_g))
            
            index = np.argmax(ydata_e)
            startpoint = [np.max(ydata_e), xdata[index], np.ptp(xdata)/6]
            para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
            perr_e = np.sqrt(np.diag(pcov_e))
        
            SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
            plt.figure()
            plt.plot(xdata, ydata_g, 'b.')
            plt.plot(xdata, ydata_e, 'r.')
            plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
            plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
            
            plt.title('SNR is {}'.format(SNR))
            plt.show()

            
            plt.figure()
            plt.plot(voltage_amp, hist_g_amp, 'b.-')
            plt.plot(voltage_amp, hist_e_amp, 'r.-')
            plt.xlabel('voltage')
            plt.ylabel('Count')
            plt.title('# {} in {}: {} photonnum'.format(k+1+m*np.size(photonnum), 
                      np.size(photonnum)*ave, photonnum[k]))
            plt.show()
        
    if ifsave:
        data = {'totvoltage_amp':totvoltage_amp, 'tothist_g_amp':tothist_g_amp,
                'tothist_e_amp':tothist_e_amp, 
                'totvoltage_correct':totvoltage_correct, 'tothist_g_correct':tothist_g_correct,
                'tothist_e_correct':tothist_e_correct, 
                'photonnum':photonnum, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)    


def AmpSweep_ge_ReadoutSNR_withave_monoThermalReadout_FTmethod(ave=5,photonnum = np.arange(0.2,8+0.3,0.5),
                            ifsave=1):
    
    totvoltage = np.zeros([ave, np.size(photonnum), 49])
    tothist_g = np.zeros([ave, np.size(photonnum), 49])
    tothist_e = np.zeros([ave, np.size(photonnum), 49])
    
    
    
    

    for m in np.arange(ave):
        for k in np.arange(np.size(photonnum)):
                    
            
            [voltage, hist_g, hist_e] = single_readout_monoThermal_FTmethod(
                    n = photonnum[k], totnum = 5000, ifplot=0, ifprint=0, 
                              phase = np.arange(0,2*np.pi,np.pi/4), ifsave=0, 
                              fJPA = 5.6205)
            
    
            
            totvoltage[m, k, :] = voltage
            tothist_g[m, k, :] = hist_g
            tothist_e[m, k, :] = hist_e
            
            
            plt.figure()
            xdata = voltage
            ydata_g = hist_g
            ydata_e = hist_e
    
            index = np.argmax(ydata_g)
            startpoint = [np.max(ydata_g), xdata[index], np.ptp(xdata)/6]
            fitfunc = Gaussian;
            para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
            perr_g = np.sqrt(np.diag(pcov_g))
            
            index = np.argmax(ydata_e)
            startpoint = [np.max(ydata_e), xdata[index], np.ptp(xdata)/6]
            para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
            perr_e = np.sqrt(np.diag(pcov_e))
        
            SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
            plt.figure()
            plt.plot(xdata, ydata_g, 'b.')
            plt.plot(xdata, ydata_e, 'r.')
            plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
            plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
            
#            plt.title('SNR is {}'.format(SNR))

            plt.title('# {} in {}: {} n {:.2f} SNR'.format(k+1+m*np.size(photonnum), 
                      np.size(photonnum)*ave, photonnum[k], SNR))
            plt.show()
        
    if ifsave:
        data = {'totvoltage':totvoltage, 'tothist_g':tothist_g,
                'tothist_e':tothist_e, 
                'photonnum':photonnum, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)            
    
        
def AmpSweep_ge_ReadoutSNR_withave_coherentReadout(ave=5, amp = np.arange(0,1.41,0.05),
                            ifsave=1, ifprint=0, numrec = 3000):
    
    totvoltage = np.zeros([ave, np.size(amp), 49])
    tothist_g = np.zeros([ave, np.size(amp), 49])
    tothist_e = np.zeros([ave, np.size(amp), 49])
    tot_cangle = np.zeros([ave, np.size(amp)])
    
    
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(ave):
        for k in np.arange(np.size(amp)):
                    
            check_readout_multiplestate(measure = 'coherent', num_steps=5, 
                                        measTime = 2000, c_amp=amp[k])
            
            _, voltage, hist_g, hist_e, angle_c = xr.compare_ge_readout_simultaneously_noSequence(
                    num_steps=5, numrec=numrec, ifplot=0, ifplotSNR=0, 
                    ifplotamp=0,ifreturn=1)
            
    
            
            totvoltage[m, k, :] = voltage
            tothist_g[m, k, :] = hist_g
            tothist_e[m, k, :] = hist_e
            tot_cangle[m, k] = angle_c
            
            plt.figure()
            plt.plot(voltage, hist_g, 'b.-')
            plt.plot(voltage, hist_e, 'r.-')
            plt.xlabel('voltage')
            plt.ylabel('Count')
            plt.title('# {} in {}: {} amp'.format(k+1+m*np.size(amp), 
                      np.size(amp)*ave, amp[k]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'totvoltage':totvoltage, 'tothist_g':tothist_g,
                'tothist_e':tothist_e, 'tot_cangle':tot_cangle,
                'amp':amp, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        
        
def AmpSweep_ge_ReadoutSNR_withave_coherentReadout_FTmethod(ave=5, amp = np.arange(0,1.41,0.05),
                            ifsave=1, ifprint=0, numrec = 4000, weight='fixed'):
    
    totvoltage = np.zeros([ave, np.size(amp), 49])
    tothist_g = np.zeros([ave, np.size(amp), 49])
    tothist_e = np.zeros([ave, np.size(amp), 49])
    
    
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(ave):
        for k in np.arange(np.size(amp)):
                    
            check_readout_multiplestate(measure = 'coherent', num_steps=5, 
                                        measTime = 2000, c_amp=amp[k])
            
            if weight == 'fixed':
            
                _, voltage, hist_g, hist_e, _, _, _ = xr.compare_ge_readout_FT_simultaneously_nosequence_weight_fixedWeight(
                            num_steps=5, fJPA = 5.6205,numrec=numrec, 
                            ifplot=0, ifplotSNR=0,ifsave=0,weightType=1,width=0.5)
            elif weight == 'notfixed':
                _, voltage, hist_g, hist_e, _ = xr.compare_ge_readout_FT_simultaneously_nosequence_weight(
                        num_steps=5, fJPA = 5.6205, numrec=numrec, 
                        ifplot=0, ifplotSNR=0, ifsave=0)
                
            elif weight == '4peak':
                _, voltage, hist_g, hist_e, _, _, _ = xr.compare_ge_readout_FT_simultaneously_nosequence_weight_fixedWeight(
                            num_steps=5, fJPA = 5.6205,numrec=numrec, 
                            ifplot=0, ifplotSNR=0,ifsave=0,weightType=2,width=0.5)
            
    
            
            totvoltage[m, k, :] = voltage
            tothist_g[m, k, :] = hist_g
            tothist_e[m, k, :] = hist_e
            
            plt.figure()
            plt.plot(voltage, hist_g, 'b.-')
            plt.plot(voltage, hist_e, 'r.-')
            plt.xlabel('voltage')
            plt.ylabel('Count')
            plt.title('# {} in {}: {} amp'.format(k+1+m*np.size(amp), 
                      np.size(amp)*ave, amp[k]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'totvoltage':totvoltage, 'tothist_g':tothist_g,
                'tothist_e':tothist_e, 
                'amp':amp, 'ave':ave}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data) 
        

def ac_Stark_shift(start=1000,width=50,freqrange=[-0.07,0.01],
                   acStark_amp=0.1,amp=0.4,state=0, readout='coherent'):
    file_length = 8000
    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    thermal_amp = 2
    
    if state == 1:
    
        pi_ge = Pulse(start=file_length-6050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
        pi_ge.phase = 90
        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    if readout == 'coherent':
        main_pulse = Pulse(start = file_length-6050,duration = 2000, amplitude= amp )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    elif readout == 'thermal':
        main_pulse = Pulse(start = file_length-6050,duration = 2000, amplitude= thermal_amp )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
       
    acStark_ge = Pulse(start=file_length-6050+start, duration=width, amplitude=acStark_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='freq', start=freqrange[0],stop=freqrange[1], initial_pulse=acStark_ge)
    acStark_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='freq', start=freqrange[0],stop=freqrange[1], initial_pulse=acStark_ge)
    
    
    
    main_pulse = Pulse(start = file_length-3000,duration = 2000, amplitude= readout_amp )
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-3000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])
##END geom


def time_sweep_acStark(singleAve=10, amp = 0.4, acStark_amp = 0.05, width = 200,
                       start=np.arange(-400,2801,200),freqrange=[-0.07,0.03],
                       ifprint=0, ifsave=1):
    num_steps = 101
       
    
    totprob = np.zeros([np.size(start), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    
    for k in np.arange(np.size(start)):
                
        ac_Stark_shift(start=start[k],width=width,freqrange=freqrange,
                       acStark_amp=acStark_amp,amp=amp)
        
        sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                       sweepparam = np.linspace(freqrange[0],freqrange[1],num_steps), ave = singleAve, 
                       ifprint = 0, ifsave = 0, ifplot=0)
        

        
        totprob[k,:] = prob
        
        plt.figure()
        plt.plot(sweepparam, prob, 'k.-')
        plt.xlabel('Detuning (GHz)')
        plt.ylabel('$P_e$')
        plt.title('# {} in {}: {} ns Start'.format(k+1, np.size(start), start[k]))
        plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'start':start}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)  
        
def time_sweep_acStark_thermal(singleAve=10, acStark_amp = 0.05, width = 200,
                       start=np.arange(-400,3000,200),freqrange=[-0.07,0.03],
                       ifprint=0, ifsave=1):
    num_steps = 101
       
    
    totprob = np.zeros([np.size(start), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    
    for k in np.arange(np.size(start)):
                
        ac_Stark_shift(start=start[k],width=width,freqrange=freqrange,
                       acStark_amp=acStark_amp, readout='thermal')
        
        sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                       sweepparam = np.linspace(freqrange[0],freqrange[1],num_steps), ave = singleAve, 
                       ifprint = 0, ifsave = 0, ifplot=0)
        

        
        totprob[k,:] = prob
        
        plt.figure()
        plt.plot(sweepparam, prob, 'k.-')
        plt.xlabel('Detuning (GHz)')
        plt.ylabel('$P_e$')
        plt.title('# {} in {}: {} ns Start'.format(k+1, np.size(start), start[k]))
        plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'start':start}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)  
        
def timeAndPower_sweep_acStark_thermal(singleAve=10, acStark_amp = 0.05, width = 200,
                       start=np.arange(-400,3000,200),freqrange=[-0.07,0.03],
                       att=np.arange(40,5,-2),ifprint=0, ifsave=1):
    num_steps = 101
       
    
    totprob = np.zeros([np.size(att), np.size(start), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    
    for k in np.arange(np.size(start)):
        
        ac_Stark_shift(start=start[k],width=width,freqrange=freqrange,
                           acStark_amp=acStark_amp, readout='thermal')
        
        for m in np.arange(np.size(att)):
                    
            vc.set_attenuation(int(att[m]))
            
            
            
            sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                           sweepparam = np.linspace(freqrange[0],freqrange[1],num_steps), ave = singleAve, 
                           ifprint = 0, ifsave = 0, ifplot=0)
            
    
            
            totprob[m, k,:] = prob
            
            plt.figure()
            plt.plot(sweepparam, prob, 'k.-')
            plt.xlabel('Detuning (GHz)')
            plt.ylabel('$P_e$')
            plt.title('# {} in {}: {} ns Start {} dB att'.format(m+1+k*np.size(att), 
                      np.size(start)*np.size(att), start[k], att[m]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'start':start, 'att':att}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)  
        
def ampAndTime_sweep_acStark(singleAve=10, amp = np.arange(0,1.41,0.2), acStark_amp = 0.03, width = 300,
                   start=np.arange(-600,2701,300),freqrange=[-0.07,0.03],ifprint=1,ifsave=1):
    num_steps = 101
       
    
    totprob = np.zeros([np.size(amp), np.size(start), num_steps])
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    for m in np.arange(np.size(amp)):
        for k in np.arange(np.size(start)):
        
                
            ac_Stark_shift(start=start[k],width=width,freqrange=freqrange,
                           acStark_amp=acStark_amp,amp=amp[m])
            
            sweepparam, prob = readout_python(num_patterns = num_steps, num_records = 200, 
                           sweepparam = np.linspace(freqrange[0],freqrange[1],num_steps), ave = singleAve, 
                           ifprint = 0, ifsave = 0, ifplot=0)
            
    
            
            totprob[m,k,:] = prob
            
            plt.figure()
            plt.plot(sweepparam, prob, 'k.-')
            plt.xlabel('Detuning (GHz)')
            plt.ylabel('$P_e$')
            plt.title('# {} in {}: {} ns Start {} amp'.format(k+1+m*np.size(start),
                      np.size(start)*np.size(amp), start[k], amp[m]))
            plt.show()
        
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'totprob':totprob,
                'start':start, 'amp':amp}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)  


def check_readout_multiplestate(measure = 'highpower', num_steps=5, 
                                measTime = 2000, c_amp = 1, phase=0,
                                ifplotsequence=1, delay=0): #this is pulsed readout to ring up and ring down cavity dfor e state
# ===========================
# g state

    file_length = 8000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    noise_amp = 2
#    pi_ge = Pulse(start=file_length-3050-pi_ge_time, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='start',start=0,stop=-sweep_time, initial_pulse=pi_ge)
#    pi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='start',start=0,stop=-sweep_time ,initial_pulse=pi_ge)
    
    if measure == 'highpower':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= readout_amp )
        ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    elif measure == 'thermal':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= noise_amp )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    elif measure == 'coherent':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= c_amp )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    elif measure == 'coherent with phase':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= c_amp*np.cos(phase) )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= c_amp*np.sin(-phase) )
        ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=main_pulse)
        
    
    ## markers
    alazar_trigger = Pulse(start=file_length-5000-delay, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-5000-300:file_length-5000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
    
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])


# ===========================
# e state    

    file_length = 8000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    pi_ge = Pulse(start=file_length-5050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
    pi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    if measure == 'highpower':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= readout_amp )
        ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    elif measure == 'thermal':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= noise_amp )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    elif measure == 'coherent':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= c_amp )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
    elif measure == 'coherent with phase':
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= c_amp*np.cos(phase) )
        ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=main_pulse)
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= c_amp*np.sin(-phase) )
        ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-5000-delay, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-5000-300:file_length-5000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
    
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=num_steps, write_binary=True)
    ringupdown_seq.num_steps = num_steps * 2
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])


def check_readout_multiplestate_differentPhase(  
                                measTime = 2000, c_amp = 1, 
                                phase=np.arange(0,2*np.pi,np.pi/4),
                                ifplotsequence=1): #this is pulsed readout to ring up and ring down cavity dfor e state
# ===========================
# g state
    num_steps = np.size(phase)
    file_length = 8000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    noise_amp = 2
#    pi_ge = Pulse(start=file_length-3050-pi_ge_time, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='start',start=0,stop=-sweep_time, initial_pulse=pi_ge)
#    pi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='start',start=0,stop=-sweep_time ,initial_pulse=pi_ge)
    
    for k in np.arange(num_steps):
        pulseLength = measTime
        pulseStart = file_length-5000

        ringupdown_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] = c_amp*np.cos(phase[k])
        ringupdown_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] = c_amp*np.sin(-phase[k])
         
    
    ## markers
    alazar_trigger = Pulse(start=file_length-5000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-5000-300:file_length-5000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
    
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])


# ===========================
# e state    

    file_length = 8000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    pi_ge = Pulse(start=file_length-5050, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
    pi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    for k in np.arange(num_steps):
        pulseLength = measTime
        pulseStart = file_length-5000

        ringupdown_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] = c_amp*np.cos(phase[k])
        ringupdown_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] = c_amp*np.sin(-phase[k])
    
    ## markers
    alazar_trigger = Pulse(start=file_length-5000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if ifplotsequence:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-5000-300:file_length-5000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
    
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=num_steps, write_binary=True)
    ringupdown_seq.num_steps = num_steps * 2
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])


def check_readout_multiplestate_ef(measure = 1, num_steps=5, measTime = 2000): #this is pulsed readout to ring up and ring down cavity dfor e state
# ===========================
# g state

    file_length = 8000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    pi_ge = Pulse(start=file_length-5050-pi_ge_time, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
    pi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    if measure == 1:
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= readout_amp )
        ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-5000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-5000-300:file_length-5000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
    
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])


# ===========================
# e state    

    file_length = 8000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = 1
    pi_ge_time = pi_ge_time_setting
    pi2_ge_time = pi2_ge_time_setting
    pi_ef_time = pi_ef_time_setting
    pi2_ef_time = pi2_ef_time_setting
#    ssm_ge = ssm_ge_setting
    ssm_ef = ssm_ef_setting
    readout_amp = 1 
    
    pi_ge = Pulse(start=file_length-5050-pi_ef_time, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
    pi_ge.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    
    pi_ef = Pulse(start=file_length-5050, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ef)
    pi_ef.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ef)
    
    
    if measure == 1:
        main_pulse = Pulse(start = file_length-5000,duration = measTime, amplitude= readout_amp )
        ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-5000, duration=500, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-5000-300:file_length-5000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
    
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=num_steps, write_binary=True)
    ringupdown_seq.num_steps = num_steps * 2
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1.5,1.5])



    
    
def readout_python_withInfo(num_records = 200, ifplot=0,
                   ifsave=0, ave = 1, ifprint=1):    
    
    
    global seqInfo
    num_patterns = seqInfo['num_steps']
    sweepparam = seqInfo['sweepparam']
    sweepparam = sweepparam.flatten()
    sweeplabel = seqInfo['sweeplabel']
    
    prob_e = np.zeros(num_patterns)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
        
    for k in np.arange(ave):
        daq_params, rec_all_raw_ave = (
                dp.run_daq_rawdata(num_patterns=num_patterns,
                                   num_records_per_pattern=num_records))
        
        _, currentprob = xr.process_readout_pattern(rec_all_raw_ave, daq_params, ifplot=0)
        
        prob_e = prob_e + currentprob/ave
        
        if ifplot:
            plt.figure()
            plt.plot(sweepparam, prob_e, 'k.-')
            plt.xlabel(sweeplabel)
            plt.ylabel('P$_e$')
    
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'prob_e':prob_e}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
        
    return sweepparam, prob_e

def readout_python(sweepparam = np.linspace(0,50,101), 
                   num_patterns = 101, num_records = 200, ifplot=0,
                   ifsave=0, ave = 1, ifprint=1):    
    
    prob_e = np.zeros(num_patterns)
    
    if ifprint == 0:
        _original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
        
    for k in np.arange(ave):
        daq_params, rec_all_raw_ave = (
                dp.run_daq_rawdata(num_patterns=num_patterns,
                                   num_records_per_pattern=num_records))
        
        _, currentprob = xr.process_readout_pattern(rec_all_raw_ave, daq_params, ifplot=0)
        
        prob_e = prob_e + currentprob
        
        if ifplot:
            plt.figure()
            plt.plot(sweepparam, prob_e/(k+1), 'k.-')
            plt.title('# {} in {}'.format(k+1, ave))
            plt.xlabel('Sweepparam')
            plt.ylabel('P$_e$')
            plt.show()
    
    prob_e = prob_e/ave
    
    if ifprint == 0:
        sys.stdout = _original_stdout
        
    if ifsave:
        data = {'sweepparam':sweepparam, 'prob_e':prob_e}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\03_single_photon',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
        
    return sweepparam, prob_e
    

def freqset_f0toe1_transition(df=0.5, offset = 0, power=5):
    rm = pyvisa.ResourceManager() 
    
    f1 = 4.806 + df
    f2 = 5.606 + df + offset

    inst = rm.open_resource('GPIB0::19::INSTR')
    print(inst.query('*IDN?'))
    inst.write('FREQ {} GHz'.format(f2))
    inst.write('POW {} dBm'.format(power))
    currentFreq = float(inst.query('FREQ?'))
    currentPow = float(inst.query('POW?'))
    
    print('Current freq is {} GHz'.format(currentFreq*1e-9))
    print('Current power is {} dBm \n'.format(currentPow))
    
    inst.close()
    
    inst2 = rm.open_resource('USB0::0x03EB::0xAFFF::471-43A6D0600-1548::INSTR')
    print(inst2.query('*IDN?'))
    inst2.write('FREQ {} GHz'.format(f1))
    inst2.write('POW {} dBm'.format(power))
    currentFreq = float(inst2.query('FREQ?'))
    currentPow = float(inst2.query('POW?'))
    
    print('Current freq is {} GHz'.format(currentFreq*1e-9))
    print('Current power is {} dBm '.format(currentPow))
    
    inst2.close()
    
def test_Agilent():
    rm = pyvisa.ResourceManager()   
    inst = rm.open_resource('GPIB0::19::INSTR')
    print(inst.query('*IDN?'))
#    inst.write('FREQ 6.106 GHz')
#    inst.write('POW 5 dBm')
    
#    inst.write('FREQ 5.60465 GHz')
#    inst.write('POW 15 dBm')
    
#    inst.write('FREQ 5.6205 GHz')
#    inst.write('POW 19 dBm')
    
    inst.write('FREQ 5.8389 GHz')
    inst.write('POW -10 dBm')
    
    currentFreq = float(inst.query('FREQ?'))
    currentPow = float(inst.query('POW?'))
    
    print('Current freq is {} GHz'.format(currentFreq*1e-9))
    print('Current power is {} dBm'.format(currentPow))
    
    inst.close()
    
def test_BNC_readout():
    rm = pyvisa.ResourceManager()   
    inst = rm.open_resource('USB0::0x03EB::0xAFFF::471-43A6D0600-1548::INSTR')
    print(inst.query('*IDN?'))
#    inst.write('FREQ 5.306 GHz')
#    inst.write('POW 7 dBm')
    
#    inst.write('FREQ 0.36 GHz')
#    inst.write('POW 4.3 dBm')
    
    inst.write('FREQ 5.6185 GHz')
    inst.write('POW 15 dBm')
    
    currentFreq = float(inst.query('FREQ?'))
    currentPow = float(inst.query('POW?'))
    
    print('Current freq is {} GHz'.format(currentFreq*1e-9))
    print('Current power is {} dBm'.format(currentPow))
    
    inst.close()
    
def test_BNC_readout_2():
    rm = pyvisa.ResourceManager()   
    inst = rm.open_resource('USB0::0x03EB::0xAFFF::421-4385A0002-0619::INSTR')
    print(inst.query('*IDN?'))
#    inst.write('FREQ 6.106 GHz')
#    inst.write('POW 7 dBm')
    
#    inst.write('FREQ 5.60465 GHz')
#    inst.write('POW 15 dBm')
    
    inst.write('FREQ 5.6185 GHz')
    inst.write('POW 15 dBm')
    
    
    
    currentFreq = float(inst.query('FREQ?'))
    currentPow = float(inst.query('POW?'))
    
    print('Current freq is {} GHz'.format(currentFreq*1e-9))
    print('Current power is {} dBm'.format(currentPow))
    
    inst.close()
    

def set_qubit_BNC_freq(setfreq=4.172):
    rm = pyvisa.ResourceManager()   
    inst = rm.open_resource('USB0::0x03EB::0xAFFF::421-4385A0002-0683::INSTR')
    print(inst.query('*IDN?'))
    
    inst.write('FREQ {} GHz'.format(setfreq))

    currentFreq = float(inst.query('FREQ?'))
    
    print('Current freq is {} GHz'.format(currentFreq*1e-9))

    inst.close()

def program_stop_button():
    def start():
        """Enable scanning by setting the global flag to True."""
        global running
        running = True

    def stop():
        """Stop scanning by setting the global flag to False."""
        global running
        running = False
        
    def scanning():
        if running:  # Only do this if the Stop button has not been clicked
            print('start scan')
            xr.compare_ge_readout_simultaneously_noSequence(num_steps=5, 
                                                            numrec=1000,ifplotSNR=0,ifplotamp=0)
            
#            readout_singleRun(ifWX = 0)
            root.after(10, scanning)
            
        # After 1 second, call scanning again (create a recursive loop)
        else:
            print('stop')
        
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    global running
    running = True    
    root = Tk()
    root.wm_attributes('-topmost', 1)
    
    root.title("Title")
    root.geometry("500x500")
    
    app = tki.Frame(root)
    app.grid()
    
    start = tki.Button(app, text="Start Scan", command=scanning)
    stop = tki.Button(app, text="Stop", command=stop)
    
    
    start.grid()
    stop.grid()

    root.mainloop()

    
            
   
#    root.destroy() 

def Gaussian(x, A, c, sigma):
    return A * np.exp(-(x-c)**2/(2*sigma**2))

# Disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore print
def enablePrint():
    sys.stdout = sys.__stdout__   


def autoscale_y(ax,margin=0.05):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    import numpy as np

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>=lo) & (xd<=hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)
    
if __name__ == '__main__':
    pass
  #  pulsed_readout_ramsey_GE()
    #rabi_gaussian()
    #Ramsey()
    #T2_echo()
    #T1()
  #  ramsey_GE()
    #T1_constant()
   # rabi_constant()
    #ramsey_constant()
    #ramseyef_constant()
    #ramseygf_constant()
    #echogf_constant()
    #rabi()
    #ramsey_GE()
    