import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
from analysis import *


def get_IQ_averages(values):
    '''
    This function takes the values input, returns your IQ averages
    '''
    I1 = values.rec_avg_vs_pats_1[0]
    Q1 = values.rec_avg_vs_pats_1[1]
    I2 = values.rec_avg_vs_pats_2[0]
    Q2 = values.rec_avg_vs_pats_1[1]

    return pd.DataFrame({'I1': I1, 'Q1':Q1, 'I2': I2, 'Q2':Q2})


def spectroscopy_plot(freq_list: list, values: dict, vert_line_value: list = None, qubit_num=1):
    """
    For a given qubit spectroscopy, plot the IQ values on side-by-side subplots with optional vertical lines.

    Args:
        freq_list: list of frequencies used in spectroscopy.
        values: dictionary from readout.
        vert_line_value: optional list with two values. The first is the vertical line position for the I channel,
                         and the second is for the Q channel.
        qubit_num: qubit number to run spectroscopy on.
    """
    IQ_vals = get_IQ_averages(values)
    IQ_vals['freq'] = freq_list

    # Create a figure with two side-by-side subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 4), dpi=150)
    
    # Plot the I channel on the first subplot
    sns.lineplot(data=IQ_vals, x='freq', y=f'I{qubit_num}', ax=axs[0])
    axs[0].set_title(f'Qubit {qubit_num} I-channel')
    axs[0].set_xlabel('Frequency')
    axs[0].set_ylabel('I Value')
    
    # Optionally add a vertical line to the I channel plot
    if vert_line_value[0] is not None and len(vert_line_value) > 0:
        axs[0].axvline(x=vert_line_value[0], color='red', linestyle='--')
    
    # Plot the Q channel on the second subplot
    sns.lineplot(data=IQ_vals, x='freq', y=f'Q{qubit_num}', ax=axs[1])
    axs[1].set_title(f'Qubit {qubit_num} Q-channel')
    axs[1].set_xlabel('Frequency')
    axs[1].set_ylabel('Q Value')
    
    # Optionally add a vertical line to the Q channel plot
    if vert_line_value[1] is not None and len(vert_line_value) > 1:
        axs[1].axvline(x=vert_line_value[1], color='red', linestyle='--')
    
    plt.tight_layout()
    plt.show()


def rabi_plot(sweep_time, num_steps, values, qubit_num=1):
    IQ_data = get_IQ_averages(values)
    Q = IQ_data[f'Q{qubit_num}']
    I = IQ_data[f'I{qubit_num}']
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        times = np.linspace(0,sweep_time/1000,num_steps)
        pi_ge_fit_vals,_,_,_ = fit_sine_decay(times,Q,guess_vals=[11,0.3,np.abs(np.max(Q)-np.min(Q)),38,Q[0]])
        pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
        print("\u03C0_ge time = {} ns".format(pi_ge))
    else:    
            times = np.linspace(0,sweep_time/1000,num_steps)
            pi_ge_fit_vals,_,_,_ = fit_sine_decay(times,I,guess_vals=[11,0.3,np.abs(np.max(I)-np.min(I)),38,I[0]])
            pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
            print("\u03C0_ge time = {} ns".format(pi_ge))
    