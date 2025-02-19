# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:07:08 2024

@author: A. Udenkwo
"""

import sys
file_dir= r"C:\Users\quantum1\Documents\Python Scripts\data_acquisition_scripts"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\quantum1\Documents\Python Scripts\sequence_generator"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Data\2020\201021_thermal_readout\Ramsey vs thermal noise"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\quantum1\Documents\Python Scripts\sequence_generator\py_sequences"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"Z:\Mark"
if file_dir not in sys.path: sys.path.append(file_dir)


from generator import *
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import traceback
import os
import old_references.nonlinear_QM_FPJPA as Fs
import wx_programs as wx
import daq.daq_programs_homo as daq_programs_homo
import analysis
import chevron
import bnc

pi = np.pi
instrument_path = r"C:\Users\quantum1\Documents\Python\instr\analyzer"
if instrument_path not in sys.path: sys.path.append(instrument_path )

analyzer_path = r"C:\Users\quantum1\Documents\Python\instr\python_interface\python_without_WX2184C"
if analyzer_path not in sys.path: sys.path.append(analyzer_path )

save_dir = r"C:\Data\2024\5plus4"
target_bnc_black_address = 'GPIB0::19::INSTR'
big_agilent_address = 'GPIB0::30::INSTR'
agilent_function_generator = 'GPIB0::10::INSTR'


#qubit bnc (BNC 1)
target_bnc_address_6 = 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR' #835 with touchscreen 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR'

#FFL1 Q1 Signal Hound


reps =500#1000#300
ro_dur = 5000#5000
## [QC,RO,pump,Qubit]

#1.5 to 0.0 for Q1, determines amplitude of pi pulse
wx_amps =[1.0,.975,.5,.5]#q1 ro wx.975,q2 sweep .1#[1.0,.975,.5,.5]##[1.0,0.85,1.5,.5]# wx_amps =[1,0.85,1.9,.5]#wx_amps =[.6,1.2,.5,.5]
#wx_amps =[1.99,.9,.5,1.99] #required for doing raman transitions, need to adjust the reaodut and qubit drive amplitudes

coupler_off_value = 0.7#.45
wx_offs = [coupler_off_value,0,0,0]


#keithley2401.set_current(0 ,step_size_mA=0.01)#Decreased by 5 mA since we started exploring the noisy qubit.


qubit_bnc =  4.600
RO_LO = 6.6247#6.575#6.6247
RO_LO_pwr = 16
ROq1 = 6.7275#6.727280##6.7278575#6.72765
ROq2 = 6.65555#6.655025#6.6557#6.65585 
ROq3 = 6.583064
ssm_geq1 =-0.111#-0.110 #GHz
ssm_geq2 =-0.152 #Ghz
ssm_efq1 =-0.2525 #-.2375
ssm_efq2 =-0.224
pi_ge_time_q1 = 51#21
pi_ge_time_q2 = 40
pi_ef_time_q1 = 77

which_qubit = 1 #qubit 1 = 1; qubit 2 = 0; both qubits = -1

if which_qubit == 1:
    pi_ge_time = pi_ge_time_q1
if which_qubit == 0:
    pi_ge_time = pi_ge_time_q2
    
    
    
qubit_1_thr=[-5000,5000] 
qubit_2_thr=[-5000,5000]#[-5000,5000]
IQ_angle_q1 = 40
IQ_angle_q2 = 95
if which_qubit==1:
    IQ_angle = IQ_angle_q1
if which_qubit==0:
    IQ_angle = IQ_angle_q2 
show_IQ=1#increasing rotates CCW
scale_Q1 = 1
#pi_ef_time = 40.46

if which_qubit == 1:
    ssm_ge = ssm_geq1
    ssm_ef = ssm_efq1
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
   # ssm_ef = ssm_ge - 0.15
    #RO_freq = ROq1
    RO_freq = RO_LO
    bnc.set_bnc_output(RO_freq,power_dBm=RO_LO_pwr, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc,power_dBm=13,bnc_addr=target_bnc_address_6)
if which_qubit == 0:
    ssm_ge = ssm_geq2
    ssm_ef = ssm_efq2
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    #RO_freq = ROq2
   # qubit_bnc =4.3#4.256
    #bnc.set_bnc_output(RO_freq,power_dBm=RO_LO_pwr, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc,power_dBm=13, bnc_addr=target_bnc_address_6)
    
if which_qubit == -1:
    #Do this for heterodyne
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    #bnc.set_bnc_output(RO_freq,power_dBm=RO_LO_pwr, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc, power_dBm=13,bnc_addr=target_bnc_address_6)

spec_ge= 0
spec_ef= 0
pis_nopi = 1
run_sweep = 0
run_sweep_AWG =0
run_rabi = 0
run_rabi_ef = 0
run_T1 = 0
run_T1_ef =0
run_ramsey =0
run_qubit_bnc_sweep = 0
run_ramsey_quantum_effiency = 0
run_ramsey_quantum_effiency_sweep = 0
run_ro_freq_sweep = 0
run_ro_amp_sweep = 0
parametric_coupling=0
parametric_coupling_time_domain=0
run_TWPA_snr_sweep = 0
run_TWPA_snr_sweep_corrected = 1


if spec_ge:
    num_steps =101
    f1=-.20
    f2=-.1
    if which_qubit == 0:  #readout on both qubits, excited q2
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=.7,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
    if which_qubit == 1: ##readout on both qubits, excited q1
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=.5,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)#.08

#spectroscopy_ef(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0)
if spec_ef:
    num_steps =201
    f1=-.3
    f2=-.2
    if which_qubit == 0:  #readout on both qubits, excited q2
        Fs.spectroscopy_ef(num_steps,ssm_ge = ssm_geq2,pi_ge = pi_ge_time_q2,ssm_start=f1,ssm_stop=f2,spec_amp= 1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
    if which_qubit == 1: ##readout on both qubits, excited q1
        Fs.spectroscopy_ef(num_steps,ssm_ge = ssm_geq1,pi_ge = pi_ge_time_q1,ssm_start=f1,ssm_stop=f2,spec_amp= .1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
 
if run_rabi_ef:#rabi_ef(num_steps=51,sweep_time=200,pi_ge = 20,ssm_ge = -0.150,ssm_ef=-0.300,ef_amp=1,ROIF1= 0,ROI2= 0,q=0)
    num_steps = 101
    sweep_time = 200 #ns
    ef_amp = 1
    if which_qubit == 0:
        Fs.rabi_ef(num_steps,sweep_time,pi_ge_time_q2,ssm_geq2,ssm_efq2,ef_amp,ROIF1,ROIF2,which_qubit)
    if which_qubit == 1:
        Fs.rabi_ef(num_steps,sweep_time,pi_ge_time_q1,ssm_geq1,ssm_efq1,ef_amp,ROIF1,ROIF2,which_qubit)
        
if pis_nopi: #run to correct two qubit rabi tomography
    num_steps = 3
    reps = 15000
    if which_qubit == 0:
        
        Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
        # Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|e>, sends population to |e>
        #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
        
    if which_qubit == 1:
        Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|g>, keeps population in ground state
        Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|e>, sends population to |e>
        #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge,ssm_ef=ssm_ef);

if run_rabi:
    num_steps = 101
    sweep_time = 200 #ns
    if which_qubit == 0:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit)
    if which_qubit == 1:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit)  

if run_T1:
    num_steps = 101
    sweep_time =100000#100000#ns
    if which_qubit == 0:
        Fs.T1_ge(num_steps,sweep_time,ssm_ge,pi_ge_time,q=which_qubit,ifload = 1,ROIF=ROIF2)
    if which_qubit == 1:
        Fs.T1_ge(num_steps,sweep_time,ssm_ge,pi_ge_time,q=which_qubit,ifload = 1,ROIF=ROIF1) 

if run_T1_ef:
    num_steps = 101
    sweep_time =100000#100000#ns
    # if which_qubit == 0:
    #     Fs.T1_ge(num_steps,sweep_time,ssm_ge,ssm_ef_q1,pi_ge_time,pi_ef_time_q1,q=which_qubit,ROIF=ROIF2)
        #T1_ef(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,ssm_ef=-0.30,pi_ge_time=20,pi_ef_time=20,q=0)
    if which_qubit == 1:
        Fs.T1_ef(num_steps,sweep_time,ssm_ge,ssm_efq1,pi_ge_time,pi_ef_time_q1,q=which_qubit,ifload = 1,ROIF=ROIF1)
        
if run_ramsey:
    piecho = 0
    num_steps = 201
    sweep_time =20000#20000 #ns
    osc_num = 0#-0.25*sweep_time/1000+6
    if which_qubit == 0:
        Fs.ramsey(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ROIF2,q=which_qubit)
    if which_qubit == 1:
        Fs.ramsey(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ROIF1,q=which_qubit)
        
if run_ramsey_quantum_effiency:
    num_steps = 101
    RO_ram_amp = 0
    ro_dur = 5000#2000
    if which_qubit == 0:
        Fs.ramsey_quantum_effiency(num_steps,ssm_ge,pi_ge_time,ROIF2,q=which_qubit, RO_ram_amp = RO_ram_amp )
    if which_qubit == 1:
        Fs.ramsey_quantum_effiency(num_steps,ssm_ge,pi_ge_time,ROIF1,q=which_qubit, RO_ram_amp = RO_ram_amp)

        
if parametric_coupling:
    num_steps = 101
    f1=-.005
    f2=-.03
    if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=.5,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
        #parametric_coupling(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0)
    if which_qubit == 1:
        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=.5,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
        
if parametric_coupling_time_domain:
    num_steps = 51
    f_parametric=-0.0705#-0.0337#-.0342
    sweep_time = 5000
    #if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
    Fs.parametric_coupling_time_domain(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_para=f_parametric,spec_amp=.5,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit,sweep_time=sweep_time)
                  
           
##Do actual heterodyne
wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
#acquire the data...
n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2,rec_readout_vs_pats_1,rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)



def make_hist(IQdata, ax_hist,labelstr):
    binwidth = 200
    lim_IQ_min = np.min([IQdata])
    lim_IQ_max = np.max([IQdata])
    
    bins_IQ = np.arange(lim_IQ_min, lim_IQ_max, binwidth)
    
    counts_IQ = ax_hist.hist(IQdata, bins=bins_IQ, histtype='step', orientation='vertical', label = labelstr)[0]
    
    gauss_xaxis = np.linspace(lim_IQ_min,lim_IQ_max, len(counts_IQ))
    
    
    return bins_IQ,counts_IQ, gauss_xaxis

###0 is I, 1 is Q
if which_qubit == 0:
    axI2= plt.axes()


    print("I histogram")
    bins_Ig2, counts_Ig2, gauss_xaxis_Ig2 = make_hist(rec_readout_vs_pats_2[0,:,0],axI2, 'I2g')
    bins_Ie2, counts_Ie2, gauss_xaxis_Ie2 = make_hist(rec_readout_vs_pats_2[0,:,1],axI2, 'I2e')

    _2, mu_Ig2, std_Ig2 =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,0], gauss_xaxis_Ig2, counts_Ig2,axI2)

    _2, mu_Ie2, std_Ie2 =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,1], gauss_xaxis_Ie2, counts_Ie2,axI2)

    axI2.legend()
    plt.show()

    axQ2= plt.axes()

    print("Q histogram")
    bins_Qg2, counts_Qg2, gauss_xaxis_Qg2 = make_hist(rec_readout_vs_pats_2[1,:,0],axQ2,'Q2g')
    bins_Qe2, counts_Qe2, gauss_xaxis_Qe2 = make_hist(rec_readout_vs_pats_2[1,:,1],axQ2, 'Q2e')

    _2, mu_Qg2, std_Qg2 =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,0], gauss_xaxis_Qg2, counts_Qg2,axQ2)

    _2, mu_Qe2, std_Qe2 =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,1], gauss_xaxis_Qe2, counts_Qe2,axQ2)

    axQ2.legend()
    plt.show()


    #Calculate SNR

    signal = np.sqrt( (mu_Ig2 - mu_Ie2)**2 + (mu_Qg2 - mu_Qe2)**2 )
    noise = (np.abs(std_Ig2) + np.abs(std_Qg2) + np.abs(std_Ie2) + np.abs(std_Qg2))/4


    SNR2 = signal/noise

    print("SNR of Qubit2:",SNR2)
    

    
    


    
#if which_qubit==1:                                                                                                                                                                                            
print("Qubit 1")
P_Q1 = prob_vs_pats_1[0];
#plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
I_Q1 = rec_avg_vs_pats_1[0];
Q_Q1 = rec_avg_vs_pats_1[1];


if show_IQ:
    if spec_ge ==1:
        freq = np.linspace(f1,f2,num_steps)
        plt.plot(freq,I_Q1);plt.title('I Q1');plt.show()
        plt.plot(freq,Q_Q1);plt.title('Q Q1');plt.show()
    else:
        plt.plot(I_Q1);plt.title('I Q1');plt.show()
        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()

#Scaled Q1

#if which_qubit==0:
print("Qubit 2")
P_Q2 = prob_vs_pats_2[0];
#    plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
I_Q2 = rec_avg_vs_pats_2[0];
Q_Q2 = rec_avg_vs_pats_2[1];
if show_IQ:
    if spec_ge ==1:
        freq = np.linspace(f1,f2,num_steps)
        plt.plot(freq,I_Q2);plt.title('I Q2');plt.show()
        plt.plot(freq,Q_Q2);plt.title('Q Q2');plt.show()
    else:
        plt.plot(I_Q2);plt.title('I Q2');plt.show()
        plt.plot(Q_Q2);plt.title('Q Q2');plt.show()

if which_qubit == 0:
    I = I_Q2; Q = Q_Q2

if which_qubit == 1:
    I = I_Q1; Q = Q_Q1



if spec_ge or spec_ef:
    freq = np.linspace(f1,f2,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        a=np.median(Q)
        if(np.max(Q)-a>a-np.min(Q)):
            freq_index = np.where(Q == np.max(Q))
        else:
            freq_index=np.where(Q==np.min(Q))
        print("Q")
        plt.plot(freq,Q)
                

    else: #    if Irange>Qrange:
        a=np.median(I)
        if(np.max(I)-a>a-np.min(I)):
                freq_index = np.where(I == np.max(I))
        else:
            freq_index=np.where(I==np.min(I))
       
        print("I")
        plt.plot(freq,I)

    ssm_ge = freq[freq_index]
    print(ssm_ge, freq_index)
    
if run_rabi or run_rabi_ef:
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        times = np.linspace(0,sweep_time/1000,num_steps)
        pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[11,0.3,np.abs(np.max(Q)-np.min(Q)),38,Q[0]])
        pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
        print("\u03C0_ge time = {} ns".format(pi_ge))
    else:    
        times = np.linspace(0,sweep_time/1000,num_steps)
        pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[11,0.3,np.abs(np.max(I)-np.min(I)),38,I[0]])
        pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
        print("\u03C0_ge time = {} ns".format(pi_ge))
    
if run_sweep:
    current_start=-20
    current_stop=-50
    current_steps = 31
    address = "GPIB0::26::INSTR"
    out_Q1, out_I1,  out_Q2, out_I2 = chevron.sweep_keithley(address , current_start,current_stop,current_steps,num_steps,reps,ro_dur,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,qubit_1_thr,qubit_2_thr)
    # save_basename = '\\keithley_sweep6_-20to-50_q2'
    # np.savetxt(save_dir+save_basename+"_I1_1000reps.txt",out_I1)
    # np.savetxt(save_dir+save_basename+"_Q1_1000reps.txt",out_Q1)
    # np.savetxt(save_dir+save_basename+"_I2_1000reps.txt",out_I2)
    # np.savetxt(save_dir+save_basename+"_Q2_1000reps.txt",out_Q2)


if run_sweep_AWG:
    voltage_start=0
    voltage_stop=-6#resistor is 300 ohm
    voltage_steps = 21
    address = agilent_function_generator
    out_Q1, out_I1,  out_Q2, out_I2 = chevron.sweep_DC_AWG(address , voltage_start,voltage_stop,voltage_steps,num_steps,reps,ro_dur,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,qubit_1_thr,qubit_2_thr)
    
if run_qubit_bnc_sweep:
    
    start_freq= 4.600 + .02
    stop_freq= 4.600 - .02
    freq_steps = 21
    pow_setting = 13
    out_Q1, out_I1,  out_Q2, out_I2 = chevron.sweep_bnc_qubit(start_freq, stop_freq, freq_steps,num_steps,reps,ro_dur,pow_setting,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,qubit_1_thr,qubit_2_thr,verbose=True)

    save_basename = '\\chevron_sweep1_q2_amp_.4'
    np.savetxt(save_dir+save_basename+"_I1reps.txt",out_I1)
    np.savetxt(save_dir+save_basename+"_Q1reps.txt",out_Q1)
    np.savetxt(save_dir+save_basename+"_I2reps.txt",out_I2)
    np.savetxt(save_dir+save_basename+"_Q2reps.txt",out_Q2)



if run_T1 or run_T1_ef:
    times = np.linspace(0,sweep_time/1000,num_steps)
    Qrange = np.max(Q)-np.min(Q)
    Irange = np.max(I)-np.min(I)
    Qrange_abs = abs(Qrange)
    Irange_abs = abs(Irange)
    
    if Qrange_abs>Irange_abs:
        T1_ge_fit_vals,error_vals,T1_fit_data,_ = analysis.fit_exp_decay(times,Q,guess_vals=[-Qrange,0.5,Q[-1]])
    else:
        T1_ge_fit_vals,error_vals,T1_fit_data,_ = analysis.fit_exp_decay(times,I,guess_vals=[-Irange,0.5,I[-1]])

    print(T1_ge_fit_vals)
    T1_ge = 1/T1_ge_fit_vals[1]
    print("T1_ge = {} \u03BCs".format(T1_ge))
    
if run_ramsey:
    times = np.linspace(0,sweep_time/1000,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[0.389417674062825 ,0.04600696770461763,-67.76128980276253,-263.4006117568909,322.4143200721203])
    if Irange>Qrange:
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[0.389417674062825 ,0.04600696770461763,-67.76128980276253,-263.4006117568909,322.4143200721203])
    T2 = 1/T2_fit_vals[1]
    print("T2* = {} \u03BCs".format(T2))
    
if run_ramsey_quantum_effiency:
    times = np.linspace(0,2*np.pi,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        fit_vals,_,_,_ = analysis.fit_sine(times,Q,guess_vals=[.16 ,100,-67.76128980276253,-263.4006117568909])
    if Irange>Qrange:
        fit_vals,_,_,_ = analysis.fit_sine(times,I,guess_vals=[.16,100,-67.76128980276253,-263.4006117568909])
    
    
    print('Angular freq:',fit_vals[0]*2*np.pi)
    print('Amp:',fit_vals[1])
    # if Qrange>Irange:
    #     fit_vals,_,_,_ = analysis.fit_sine_fix_freq(times,Q,guess_vals=[100,-67.76128980276253,1500],fixed_freq=0.15652176278260993)
    # if Irange>Qrange:
    #     fit_vals,_,_,_ = analysis.fit_sine_fix_freq(times,I,guess_vals=[100,-67.76128980276253,1500],fixed_freq=0.15652176278260993)

      
if run_ro_freq_sweep:
    ro_start =6.65450 #6.7265#6.6540
    ro_stop = 6.6565#6.7285#6.6550
    ro_steps = 21
    ro_freq = np.linspace(ro_start, ro_stop,ro_steps)
    num_steps = 3
    reps = 15000
    data_snr=np.zeros(ro_steps)
    if which_qubit == 0:  
        for i in range(len(ro_freq)):
            ROq2 = ro_freq[i]
            ROIF2 = ROq2 - RO_LO
            Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
            # Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
            # Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|e>, sends population to |e>
            ##Do actual heterodyne
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
            #acquire the data...
            n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2,rec_readout_vs_pats_1,rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
            axI= plt.axes()


            print("I histogram")
            bins_Ig, counts_Ig, gauss_xaxis_Ig = make_hist(rec_readout_vs_pats_2[0,:,0],axI, 'Ig')
            bins_Ie, counts_Ie, gauss_xaxis_Ie = make_hist(rec_readout_vs_pats_2[0,:,1],axI, 'Ie')

            _, mu_Ig, std_Ig =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,0], gauss_xaxis_Ig, counts_Ig,axI)

            _, mu_Ie, std_Ie =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,1], gauss_xaxis_Ie, counts_Ie,axI)

            axI.legend()
            plt.show()

            axQ= plt.axes()

            print("Q histogram")
            bins_Qg, counts_Qg, gauss_xaxis_Qg = make_hist(rec_readout_vs_pats_2[1,:,0],axQ,'Qg')
            bins_Qe, counts_Qe, gauss_xaxis_Qe = make_hist(rec_readout_vs_pats_2[1,:,1],axQ, 'Qe')

            _, mu_Qg, std_Qg =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,0], gauss_xaxis_Qg, counts_Qg,axQ)

            _, mu_Qe, std_Qe =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,1], gauss_xaxis_Qe, counts_Qe,axQ)

            axQ.legend()
            plt.show()


            #Calculate SNR

            signal = np.sqrt( (mu_Ig - mu_Ie)**2 + (mu_Qg - mu_Qe)**2 )
            noise = (np.abs(std_Ig) + np.abs(std_Qg) + np.abs(std_Ie) + np.abs(std_Qg))/4

            SNR = signal/noise
            print("SNR of Qubit2:",SNR)
            data_snr[i] = SNR
        plt.plot(ro_freq,data_snr)
          
 


if run_ro_amp_sweep:
    amp_start =.1
    amp_stop = .6
    amp_steps = 11
    ro_amp = np.linspace(amp_start, amp_stop,amp_steps)
    num_steps = 3
    data_snr=np.zeros(amp_steps)
    reps = 15000
    if which_qubit == 0: 
   
       for i in range(len(ro_amp)):
           wx_amps =[1.0,ro_amp[i],.5,.5]
           # Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
           # bnc.set_bnc_output(4.8,power_dBm=ro_amp[i], bnc_addr=big_agilent_address)   
           # Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
           # Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|e>, sends population to |e>
           ##Do actual heterodyne
           wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
           #acquire the data...
           n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2,rec_readout_vs_pats_1,rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           axI= plt.axes()


           print("I histogram")
           bins_Ig, counts_Ig, gauss_xaxis_Ig = make_hist(rec_readout_vs_pats_2[0,:,0],axI, 'Ig')
           bins_Ie, counts_Ie, gauss_xaxis_Ie = make_hist(rec_readout_vs_pats_2[0,:,1],axI, 'Ie')

           _, mu_Ig, std_Ig =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,0], gauss_xaxis_Ig, counts_Ig,axI)

           _, mu_Ie, std_Ie =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,1], gauss_xaxis_Ie, counts_Ie,axI)

           axI.legend()
           plt.show()

           axQ= plt.axes()

           print("Q histogram")
           bins_Qg, counts_Qg, gauss_xaxis_Qg = make_hist(rec_readout_vs_pats_2[1,:,0],axQ,'Qg')
           bins_Qe, counts_Qe, gauss_xaxis_Qe = make_hist(rec_readout_vs_pats_2[1,:,1],axQ, 'Qe')

           _, mu_Qg, std_Qg =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,0], gauss_xaxis_Qg, counts_Qg,axQ)

           _, mu_Qe, std_Qe =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,1], gauss_xaxis_Qe, counts_Qe,axQ)

           axQ.legend()
           plt.show()


           #Calculate SNR

           signal = np.sqrt( (mu_Ig - mu_Ie)**2 + (mu_Qg - mu_Qe)**2 )
           noise = (np.abs(std_Ig) + np.abs(std_Qg) + np.abs(std_Ie) + np.abs(std_Qg))/4

           SNR = signal/noise
           print("SNR of Qubit:",SNR)
           data_snr[i] = SNR
       plt.plot(ro_amp,data_snr) 

       plt.plot(ro_amp,data_snr)
       plt.show()
    analysis.fit_line(ro_amp,data_snr)    

if run_ramsey_quantum_effiency_sweep:
    RO_ram_amp_start =0
    RO_ram_amp_stop = .1
    RO_ram_amp_steps = 21
    RO_ram_ro_amp = np.linspace(RO_ram_amp_start, RO_ram_amp_stop,RO_ram_amp_steps)
    ramsey_amp=np.zeros(RO_ram_amp_steps)
    num_steps = 101
    reps= 5000
    ro_dur = 5000#2000
    ramsey_I,ramsey_Q=np.zeros((2,RO_ram_amp_steps,num_steps))
    fixed_freq=0.15680202279076202
    if which_qubit == 0:
        for i in range(len(RO_ram_ro_amp)):
            RO_ram_amp = RO_ram_ro_amp[i]
            Fs.ramsey_quantum_effiency(num_steps,ssm_ge,pi_ge_time,ROIF2,q=which_qubit, RO_ram_amp = RO_ram_amp )
            ##Do actual heterodyne
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
            #acquire the data...
            n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2,rec_readout_vs_pats_1,rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
            print("Qubit 1")
            P_Q1 = prob_vs_pats_1[0];
            #plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
            I_Q1 = rec_avg_vs_pats_1[0];
            Q_Q1 = rec_avg_vs_pats_1[1];
            plt.plot(I_Q1);plt.title('I Q1');plt.show()
            plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
            print("Qubit 2")
            P_Q2 = prob_vs_pats_2[0];
            #    plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
            I_Q2 = rec_avg_vs_pats_2[0];
            Q_Q2 = rec_avg_vs_pats_2[1];
            plt.plot(I_Q2);plt.title('I Q2');plt.show()
            plt.plot(Q_Q2);plt.title('Q Q2');plt.show()
            I = I_Q2; Q = Q_Q2
            times = np.linspace(0,2*np.pi,num_steps)
            Qrange = abs(np.max(Q)-np.min(Q))
            Irange = abs(np.max(I)-np.min(I))
            if Qrange>Irange:
                fit_vals,_,_,_ = analysis.fit_sine_fix_freq(times,Q,guess_vals=[100,-67.76128980276253,1500],fixed_freq=fixed_freq)
            if Irange>Qrange:
                fit_vals,_,_,_ = analysis.fit_sine_fix_freq(times,I,guess_vals=[100,-67.76128980276253,1500],fixed_freq=fixed_freq)
            
            print('Angular freq:',fixed_freq*2*np.pi)
            print('Amp:',fit_vals[0])
            ramsey_Q[i]=Q
            ramsey_I[i]=I
            ramsey_amp[i] = abs(fit_vals[0])  
        plt.plot(RO_ram_ro_amp,ramsey_amp)
        plt.show()
        for i in range(len(RO_ram_ro_amp)):
            plt.plot(times,ramsey_I[i])
        plt.title('I Q2')
        plt.show()
        for i in range(len(RO_ram_ro_amp)):
            plt.plot(times,ramsey_Q[i])
        plt.title('Q Q2')
        plt.show()
        save_basename = '\\ramsey_quantum_effiency'
        np.savetxt(save_dir+save_basename+"_ramsey_amp.txt",ramsey_amp)
        np.savetxt(save_dir+save_basename+"RO_ram_ro_amp.txt",RO_ram_ro_amp)
        Gaussian_fit_vals,_,_,_ = analysis.fit_gaussian_points(RO_ram_ro_amp,ramsey_amp,guess_vals=[300,0,1])
        

            


                
        
        
    if which_qubit == 1:
        for i in range(len(RO_ram_ro_amp)):
            RO_ram_amp = RO_ram_ro_amp[i]
            Fs.ramsey_quantum_effiency(num_steps,ssm_ge,pi_ge_time,ROIF1,q=which_qubit, RO_ram_amp = RO_ram_amp) 
            ##Do actual heterodyne
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
            #acquire the data...
            n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2,rec_readout_vs_pats_1,rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
            print("Qubit 1")
            P_Q1 = prob_vs_pats_1[0];
            #plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
            I_Q1 = rec_avg_vs_pats_1[0];
            Q_Q1 = rec_avg_vs_pats_1[1];
            plt.plot(I_Q1);plt.title('I Q1');plt.show()
            plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
            print("Qubit 2")
            P_Q2 = prob_vs_pats_2[0];
            #    plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
            I_Q2 = rec_avg_vs_pats_2[0];
            Q_Q2 = rec_avg_vs_pats_2[1];
            plt.plot(I_Q2);plt.title('I Q2');plt.show()
            plt.plot(Q_Q2);plt.title('Q Q2');plt.show()
            I = I_Q1; Q = Q_Q1
            times = np.linspace(0,2*np.pi,num_steps)
            Qrange = abs(np.max(Q)-np.min(Q))
            Irange = abs(np.max(I)-np.min(I))
            if Qrange>Irange:
                fit_vals,_,_,_ = analysis.fit_sine_fix_freq(times,Q,guess_vals=[100,-67.76128980276253,1500],fixed_freq=fixed_freq)
            if Irange>Qrange:
                fit_vals,_,_,_ = analysis.fit_sine_fix_freq(times,I,guess_vals=[100,-67.76128980276253,1500],fixed_freq=fixed_freq)
            
            print('Angular freq:',fixed_freq*2*np.pi)
            print('Amp:',fit_vals[0])
            ramsey_Q[i]=Q
            ramsey_I[i]=I
            ramsey_amp[i] = abs(fit_vals[0])  
        plt.plot(RO_ram_ro_amp,ramsey_amp)
        plt.show()
        for i in range(len(RO_ram_ro_amp)):
            plt.plot(times,ramsey_I[i])
        plt.title('I Q1')
        plt.show()
        for i in range(len(RO_ram_ro_amp)):
            plt.plot(times,ramsey_Q[i])
        plt.title('Q Q1')
        plt.show()
       
# if run_TWPA_snr_sweep:
#     amp_start =-10
#     amp_stop = 2
#     amp_steps = 50
#     freq_start = 8
#     freq_stop = 9
#     freq_steps = 101
#     twpa_amp = np.linspace(amp_start, amp_stop,amp_steps)
#     twpa_freq = np.linspace(freq_start, freq_stop,freq_steps)
#     num_steps = 3
#     data_snr=np.zeros(amp_steps)
#     reps = 15000
#     if which_qubit == 0: 
   
#        for i in range(len(twpa_amp)):
#            for j in range(len(twpa_freq)):
#                bnc.set_bnc_output(twpa_freq[j],power_dBm=twpa_amp[i], bnc_addr=big_agilent_address)   
#                wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
#                n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2,rec_readout_vs_pats_1,rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
#                axI= plt.axes()


#                print("I histogram")
#                bins_Ig, counts_Ig, gauss_xaxis_Ig = make_hist(rec_readout_vs_pats_2[0,:,0],axI, 'Ig')
#                bins_Ie, counts_Ie, gauss_xaxis_Ie = make_hist(rec_readout_vs_pats_2[0,:,1],axI, 'Ie')

#                _, mu_Ig, std_Ig =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,0], gauss_xaxis_Ig, counts_Ig,axI)

#                _, mu_Ie, std_Ie =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0,:,1], gauss_xaxis_Ie, counts_Ie,axI)

#                axI.legend()
#                plt.show()

#                axQ= plt.axes()

#                print("Q histogram")
#                bins_Qg, counts_Qg, gauss_xaxis_Qg = make_hist(rec_readout_vs_pats_2[1,:,0],axQ,'Qg')
#                bins_Qe, counts_Qe, gauss_xaxis_Qe = make_hist(rec_readout_vs_pats_2[1,:,1],axQ, 'Qe')

#                _, mu_Qg, std_Qg =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,0], gauss_xaxis_Qg, counts_Qg,axQ)

#                _, mu_Qe, std_Qe =analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1,:,1], gauss_xaxis_Qe, counts_Qe,axQ)

#                axQ.legend()
#                plt.show()


#                #Calculate SNR
    
#                signal = np.sqrt( (mu_Ig - mu_Ie)**2 + (mu_Qg - mu_Qe)**2 )
#                noise = (np.abs(std_Ig) + np.abs(std_Qg) + np.abs(std_Ie) + np.abs(std_Qg))/4
    
#                SNR = signal/noise
#                print("SNR of Qubit2:",SNR)
#                data_snr[i] = SNR
#        plt.plot(ro_amp,data_snr) 

#        plt.plot(ro_amp,data_snr)
#        plt.show()
#     analysis.fit_line(ro_amp,data_snr)                   
    
if run_TWPA_snr_sweep:
    amp_start = -6
    amp_stop = 2
    amp_steps = 21
    freq_start = 7
    freq_stop = 9
    freq_steps = 201
    twpa_amp = np.linspace(amp_start, amp_stop, amp_steps)
    twpa_freq = np.linspace(freq_start, freq_stop, freq_steps)
    num_steps = 3
    data_snr = np.zeros((amp_steps, freq_steps))  # 2D array for SNR
    reps = 15000
    if which_qubit == 0:
   
        for j in range(len(twpa_freq)):
            for i in range(len(twpa_amp)):
                bnc.set_bnc_output(twpa_freq[j], power_dBm=twpa_amp[i], bnc_addr=big_agilent_address)   
                wx.wx_set_and_amplitude_and_offset(amp=wx_amps, offset=wx_offs)
                n_vs_pats_1, n_vs_pats_2, rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2, rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2, prob_vs_pats_1, prob_vs_pats_2, n_readout_1, n_readout_2, rec_readout_vs_pats_1, rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1, ROIF2, deg_1=IQ_angle_q1, deg_2=IQ_angle_q2, num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, qubit_1_thr=qubit_1_thr, qubit_2_thr=qubit_2_thr, verbose=True)
                axI = plt.axes()
                print("I histogram")
                bins_Ig, counts_Ig, gauss_xaxis_Ig = make_hist(rec_readout_vs_pats_2[0, :, 0], axI, 'Ig')
                bins_Ie, counts_Ie, gauss_xaxis_Ie = make_hist(rec_readout_vs_pats_2[0, :, 1], axI, 'Ie')

                _, mu_Ig, std_Ig = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0, :, 0], gauss_xaxis_Ig, counts_Ig, axI)
                _, mu_Ie, std_Ie = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0, :, 1], gauss_xaxis_Ie, counts_Ie, axI)

                axI.legend()
                plt.show()

                axQ = plt.axes()
                print("Q histogram")
                bins_Qg, counts_Qg, gauss_xaxis_Qg = make_hist(rec_readout_vs_pats_2[1, :, 0], axQ, 'Qg')
                bins_Qe, counts_Qe, gauss_xaxis_Qe = make_hist(rec_readout_vs_pats_2[1, :, 1], axQ, 'Qe')

                _, mu_Qg, std_Qg = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1, :, 0], gauss_xaxis_Qg, counts_Qg, axQ)
                _, mu_Qe, std_Qe = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1, :, 1], gauss_xaxis_Qe, counts_Qe, axQ)

                axQ.legend()
                plt.show()

                # Calculate SNR
                signal = np.sqrt((mu_Ig - mu_Ie)**2 + (mu_Qg - mu_Qe)**2)
                noise = (np.abs(std_Ig) + np.abs(std_Qg) + np.abs(std_Ie) + np.abs(std_Qe)) / 4

                SNR = signal / noise
                print(f"SNR of Qubit2 at Amp={twpa_amp[i]}, Freq={twpa_freq[j]}: {SNR}")

                data_snr[i, j] = SNR  # Store SNR in the 2D array

        # Generate 2D plot for SNR values
        plt.imshow(data_snr, aspect='auto', extent=[freq_start, freq_stop, amp_start, amp_stop], origin='lower')
        plt.colorbar(label='SNR')
        plt.xlabel('TWPA Frequency (GHz)')
        plt.ylabel('TWPA Amplitude (dBm)')
        plt.title('SNR vs TWPA Frequency and Amplitude')
        plt.show()

        # Optionally fit and plot line for SNR if needed (you might want to adjust this part)
        # For example, you can fit a line to slices in the SNR matrix
        
if run_TWPA_snr_sweep_corrected:
    # Experiment parameters
    amp_start = -6
    amp_stop = 2
    amp_steps = 41  # Number of amp points
    freq_start = 7
    freq_stop = 9
    freq_steps = 41  # Number of freq points
    twpa_amp = np.linspace(amp_start, amp_stop, amp_steps)
    twpa_freq = np.linspace(freq_start, freq_stop, freq_steps)
    # rem_twp=twpa_freq[twpa_freq>7.84]
    num_steps = 3
    reps = 15000
    # File to save the results
    # snr_output_folder = "ramsey_quantum_effiency"
    output_folder=r"C:/Data/2024/5plus4/ramsey_quantum_effiency"
    snr_output_file = os.path.join(output_folder,"snr_experiment_data7,9.csv")

    columns = ["freq", "amp", "SNR"]

    # Generate the freq and amp arrays
    twpa_amp = np.linspace(amp_start, amp_stop, amp_steps)
    twpa_freq = np.linspace(freq_start, freq_stop, freq_steps)

    # Attempt to load existing data
    try:
        data_df = pd.read_csv(snr_output_file)
        print(f"Loaded existing data from {snr_output_file}")
    except FileNotFoundError:
        data_df = pd.DataFrame(columns=columns)

    # Main experiment loop
    try:
        for freq in twpa_freq:
            # Check if the current freq is already completed
            if freq in data_df["freq"].unique():
                print(f"Freq {freq} already completed, skipping...")
                continue

            # Store the results for the current freq
            freq_data = []

            try:
                for amp in twpa_amp:
                    bnc.set_bnc_output(freq, power_dBm=amp, bnc_addr=big_agilent_address)   
                    wx.wx_set_and_amplitude_and_offset(amp=wx_amps, offset=wx_offs)
                    n_vs_pats_1, n_vs_pats_2, rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2, rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2, prob_vs_pats_1, prob_vs_pats_2, n_readout_1, n_readout_2, rec_readout_vs_pats_1, rec_readout_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1, ROIF2, deg_1=IQ_angle_q1, deg_2=IQ_angle_q2, num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, qubit_1_thr=qubit_1_thr, qubit_2_thr=qubit_2_thr, verbose=True)
                    axI = plt.axes()
                    print("I histogram")
                    bins_Ig, counts_Ig, gauss_xaxis_Ig = make_hist(rec_readout_vs_pats_2[0, :, 0], axI, 'Ig')
                    bins_Ie, counts_Ie, gauss_xaxis_Ie = make_hist(rec_readout_vs_pats_2[0, :, 1], axI, 'Ie')

                    _, mu_Ig, std_Ig = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0, :, 0], gauss_xaxis_Ig, counts_Ig, axI)
                    _, mu_Ie, std_Ie = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[0, :, 1], gauss_xaxis_Ie, counts_Ie, axI)

                    axI.legend()
                    plt.show()

                    axQ = plt.axes()
                    print("Q histogram")
                    bins_Qg, counts_Qg, gauss_xaxis_Qg = make_hist(rec_readout_vs_pats_2[1, :, 0], axQ, 'Qg')
                    bins_Qe, counts_Qe, gauss_xaxis_Qe = make_hist(rec_readout_vs_pats_2[1, :, 1], axQ, 'Qe')

                    _, mu_Qg, std_Qg = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1, :, 0], gauss_xaxis_Qg, counts_Qg, axQ)
                    _, mu_Qe, std_Qe = analysis.fit_gaussian_no_plot(rec_readout_vs_pats_2[1, :, 1], gauss_xaxis_Qe, counts_Qe, axQ)

                    axQ.legend()
                    plt.show()

                    # Calculate SNR
                    signal = np.sqrt((mu_Ig - mu_Ie)**2 + (mu_Qg - mu_Qe)**2)
                    noise = (np.abs(std_Ig) + np.abs(std_Qg) + np.abs(std_Ie) + np.abs(std_Qe)) / 4

                    SNR = signal / noise
                    print(f"SNR of Qubit2 at Amp={amp}, Freq={freq}: {SNR}")
                   
                    freq_data.append({"freq": freq, "amp": amp, "SNR": SNR})

                # Save the current freq's results
                freq_df = pd.DataFrame(freq_data)
                data_df = pd.concat([data_df, freq_df], ignore_index=True)
                data_df.to_csv(snr_output_file, index=False)  # Save to file after each freq
                print(f"Saved data for freq={freq}")

            except Exception as e:
                # Handle errors for a single freq
                print(f"Error at freq={freq}: {e}")
                print(traceback.format_exc())
                raise e  # Stop execution and save progress

    except Exception as e:
        # Handle errors in the main loop
        print(f"Experiment interrupted: {e}")
        print("Saving partial data...")
        data_df.to_csv(snr_output_file, index=False)
        print(f"Partial data saved to {snr_output_file}")
        raise e

    print("Experiment completed successfully.")
    pivot_data = data_df.pivot(index="amp", columns="freq", values="SNR")

    # Convert the pivot table to a NumPy array for plotting
    data_snr = pivot_data.values
    plt.imshow(data_snr, aspect='auto', extent=[7, 9, amp_start, amp_stop], origin='lower')
    plt.colorbar(label='SNR')
    plt.xlabel('TWPA Frequency (GHz)')
    plt.ylabel('TWPA Amplitude (dBm)')
    plt.title('SNR vs TWPA Frequency and Amplitude')
    plt.show()
# data_df = data_df.drop_duplicates(subset=["freq","amp"],keep="first")
# data_df = data_df.sort_values(by=["freq","amp"]).reset_index(drop=True)

# print(data_df["freq"].unique())


# duplicates=data_df.duplicated(subset=["freq","amp"])
# print(data_df[duplicates])
# snr_output_file = "new_snr_experiment_data.csv"
# data_df.to_csv(snr_output_file, index=False)