# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:56:59 2024

@author: quantum1
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

import vna_analysis as vna
import wx_programs as wx
import daq_programs_homo
import analysis
import chevron
import bnc
import keithley2401


keithley_current = 2.1444 #at half a magnetic flux quantum 
keithley_address = 'GPIB0::24::INSTR'
big_agilent_address = 'GPIB0::30::INSTR'

run_initial_trace = 1
run_current_sweep = 0
run_pump_sweep = 0
run_gain_profile = 0

if run_initial_trace:
    keithley2401.set_current(keithley_address, 2.1444, step_size_mA=0.01)
    time.sleep(5)
    sdata = vna.get_smith_data(restart_average=True)
    freqs = sdata.freqs
    real_part = sdata.real
    imag_part = sdata.imag
    lin_mag = sdata.lin_mag
    log_mag = sdata.log_mag
    plt.plot(freqs,log_mag)
    plt.show()
    save_basename = '\\data_initial_trace_02_18'
    save_dir_NRC = r'C:\\Data\\2025\\TWPA_Calibration'
    np.savetxt(save_dir_NRC+save_basename+".txt",log_mag)
    

if run_current_sweep:
    start_current = 1.8
    stop_current = 2.5
    num_points = 41
    
    sweep_vals = np.linspace(start_current, stop_current, num_points)
    log_mag=np.zeros( ( num_points, 501))
    phase_deg=np.zeros( (num_points, 501))
    for i,ss in enumerate(sweep_vals):
        keithley2401.set_current(keithley_address, ss, step_size_mA=0.01)
        time.sleep(10)
        sdata = vna.get_smith_data()
        freqs = sdata.freqs
        log_mag[i] = sdata.log_mag
        phase_deg[i] = sdata.phase_deg
        # plt.plot(log_mag[i])
        # plt.show()
        # plt.plot(phase_deg[i])
        # plt.show()
    plt.imshow(log_mag, extent=[4,8,stop_current,start_current],aspect='auto' )
    plt.show()
    plt.imshow(phase_deg, extent=[4,8,stop_current,start_current],aspect='auto' )
    plt.show()
    save_basename = '\\current_sweep_TWPA_trans_current_1.8,2.5_mA,4,8GHz'
    save_dir_NRC = r'C:\\Data\\2025\\TWPA_Calibration'
    np.savetxt(save_dir_NRC+save_basename+"_mag.txt",log_mag)
    np.savetxt(save_dir_NRC+save_basename+"_deg.txt",phase_deg)
    
if run_pump_sweep:
    start_freq =7#4.0
    stop_freq = 9#5.5
    freq_step = .05
    num_freq = abs(int((stop_freq-start_freq)/freq_step))+1
    start_pow = 2#-2
    stop_pow = -8
    pow_step =.2
    num_pow = abs(int((stop_pow-start_pow)/pow_step))+1
    pow_list = np.linspace(start_pow, stop_pow, num_pow)
    freq_list = np.linspace(start_freq, stop_freq, num_freq)
    average_mag = np.zeros( ( num_freq, num_pow))
    keithley2401.set_current(keithley_address, keithley_current, step_size_mA=0.01)
    time.sleep(1)
    save_basename = '\\pump_sweep_TWPA_trans_freq7,9GHz,amp2,-8dBm'
    save_dir_NRC = r'C:\\Data\\2025\\TWPA_Calibration'
    for i,freq in enumerate(freq_list):
        for j,power in enumerate(pow_list):
            target_bnc_address =  big_agilent_address
            bnc.set_bnc_output(freq,power, bnc_addr=target_bnc_address)
            time.sleep(1.5)
            sdata = vna.get_smith_data()
            freqs = sdata.freqs
            log_mag = sdata.log_mag
            average_mag[i,j]=np.mean(log_mag)      
            plt.plot(freqs,log_mag)
            plt.show()
    plt.imshow(average_mag, extent=[start_pow,stop_pow,stop_freq,start_freq],aspect='auto' )
    plt.colorbar(label='Color Scale')
    plt.ylabel("Pump Frequency")
    plt.xlabel("Pump Amplitude")
    plt.show()
    max_index_mag = np.unravel_index(np.argmax(average_mag),average_mag.shape)
    print(freq_list[max_index_mag[0]])
    print(pow_list[max_index_mag[1]])
    print(average_mag[max_index_mag[0],max_index_mag[1]])
    np.savetxt(save_dir_NRC+save_basename+".txt",average_mag)

      
if run_gain_profile:
    freq = 4.5
    amp = -5.4
    target_bnc_address =  big_agilent_address
    bnc.set_bnc_output(freq,amp, bnc_addr=target_bnc_address)
    time.sleep(10)
    sdata = vna.get_smith_data()
    freqs = sdata.freqs
    log_mag_pump = sdata.log_mag
    plt.plot(freqs,log_mag_pump)
    plt.show()
    

