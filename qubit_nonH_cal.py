import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import sys
import hardware_control.wx_programs as wx
import hardware_control.bnc as bnc
from classes.generator import *
from hardware_control.hardware_config import *
from experiment_configuration.values import *
from classes.qubit_class import *
import daq.daq_programs_homo as daq
import seaborn as sns
import standard_sequences.rabi as rabi
import analysis.plotting as plotting
import analysis.analysis as analysis
import standard_sequences.sweeppiamp as sweeppiamp
from skopt import gp_minimize
from skopt.space import Categorical
from skopt.callbacks import CheckpointSaver
import classifiers.classifier as classifier
import pickle
import standard_sequences.parametric_coupling as parametric_coupling
from datetime import datetime

def save_optimization_results(params_dict, experiment_name, save_dir=None):
    """Save optimization results to a dated text file"""
    if save_dir is None:
        save_dir = os.getcwd()
    
    log_dir = os.path.join(save_dir, 'optimization_logs')
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f"optimization_parameters_{timestamp}.txt"
    filepath = os.path.join(log_dir, filename)
    
    # Check if file exists and read existing content
    existing_params = {}
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if ':' in line:
                    key, value = line.split(':', 1)
                    existing_params[key.strip()] = value.strip()
    
    # Update with new parameters
    for key, value in params_dict.items():
        if experiment_name:
            param_name = f"{experiment_name}_{key}"
        else:
            param_name = key
        existing_params[param_name] = value
    
    # Write all parameters to file
    with open(filepath, 'w') as f:
        f.write(f"Optimization Parameters\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        # Write parameters in sorted order
        for key in sorted(existing_params.keys()):
            f.write(f"{key}: {existing_params[key]}\n")
    
    print(f"Results saved to {filepath}")
    return filepath

def run_rabi(
    q1,
    q2,
    general_vals_dict,
    num_steps: int,
    sweep_time: float,
    reps: int,
    manifold: str,
    save_dir: str = None
):
    """
    Runs a single instance of the geRabi experiment
    and processes the resulting IQ data to compute probabilities vs. time.

    Returns:
        tuple: A tuple containing:
            - values: Raw DAQ values
            - IQ_df: DataFrame with IQ data
            - rabi_amp: Fitted Rabi amplitude
            - pi_time: Fitted pi pulse time in ns
    """
    if manifold == "ge":
        rabi.rabi_ge(q1, q2, general_vals_dict, num_steps, sweep_time)
    elif manifold == "ef":
        rabi.rabi_ef(q1, q2, general_vals_dict, num_steps, sweep_time)
    else:
        raise ValueError(f"Manifold {manifold} not supported")
    
    # Run the experiment
    wx.wx_set_and_amplitude_and_offset(
        amp=general_vals_dict["wx_amps"], offset=general_vals_dict["wx_offs"]
    )
    # Acquire the raw IQ data
    values = daq.run_daq_het_2q(
        q1, q2, num_patterns=num_steps, num_records_per_pattern=reps, verbose=True
    )

    # Retrieve raw IQ data from the acquired values
    I1_raw = values.rec_readout_1[0]
    Q1_raw = values.rec_readout_1[1]
    I2_raw = values.rec_readout_2[0]
    Q2_raw = values.rec_readout_2[1]

    # Build a DataFrame from the IQ data for classification
    IQ_df = pd.DataFrame({"I1": I1_raw, "Q1": Q1_raw, "I2": I2_raw, "Q2": Q2_raw})
    
    # Get average IQ values for fitting
    IQ_data = plotting.get_IQ_averages(values)
    Q = IQ_data[f"Q1"]  # Using Q1 for fitting
    I = IQ_data[f"I1"]  # Using I1 for fitting
    
    # Determine which channel has larger range for fitting
    Qrange = abs(np.max(Q) - np.min(Q))
    Irange = abs(np.max(I) - np.min(I))
    
    # Create time array for fitting
    times = np.linspace(0, sweep_time / 1000, num_steps)
    
    # Fit the data using the channel with larger range
    if Qrange > Irange:
        fit_data = Q
    else:
        fit_data = I
        
    # Perform the Rabi fit
    pi_ge_fit_vals, _, _, _ = analysis.fit_sine_decay(
        times, fit_data, guess_vals=[11, 0.3, np.abs(np.max(fit_data) - np.min(fit_data)), 38, fit_data[0]]
    )
    
    # Extract Rabi amplitude and pi time
    rabi_amp = pi_ge_fit_vals[2]  # Amplitude from fit
    pi_time = abs((1 / 2 / pi_ge_fit_vals[0]) * 1000)  # Convert to ns
    
    # Save results with clearer parameter names
    rabi_results = {
        f'rabi_{manifold}_amplitude': rabi_amp,
        f'rabi_{manifold}_pi_time': pi_time,
        f'rabi_{manifold}_num_steps': num_steps,
        f'rabi_{manifold}_sweep_time': sweep_time,
        f'rabi_{manifold}_reps': reps,
        f'rabi_{manifold}_fit_frequency': pi_ge_fit_vals[0],
        f'rabi_{manifold}_fit_decay': pi_ge_fit_vals[1],
        f'rabi_{manifold}_fit_offset': pi_ge_fit_vals[4]
    }
    save_optimization_results(rabi_results, '', save_dir)
    
    return values, IQ_df, rabi_amp, pi_time

def pi_tune(
    q1: object,
    q2: object,
    general_vals_dict: dict,
    num_steps: int,
    reps: int,
    amp: float,
    pi_time: int,
    manifold: str = "ge",
):
    """
    Runs a single instance of the pi pulse tuning experiment for either ge or ef manifold
    and processes the resulting IQ data.

    Returns:
        tuple: (values, IQ_df, df_prob) containing raw values, IQ data, and state probabilities
    """
    # Run the experiment based on manifold
    if manifold == "ge":
        sweeppiamp.pi_ge_gaussian_tune_iq(
            q1,
            q2,
            general_vals_dict,
            num_steps=3,
            amp=amp,
            pi_ge_time=pi_time,
            file_length=16000,
            verbose=False,
        )
    elif manifold == "ef":
        sweeppiamp.pi_ef_gaussian_tune_iq(
            q1,
            q2,
            general_vals_dict,
            num_steps=3,
            amp=amp,
            pi_ef_time=pi_time,
            file_length=16000,
            verbose=False,
        )
    else:
        raise ValueError(f"Manifold {manifold} not supported")

    wx.wx_set_and_amplitude_and_offset(
        amp=general_vals_dict["wx_amps"], offset=general_vals_dict["wx_offs"]
    )
    
    # Acquire the raw IQ data
    values = daq.run_daq_het_2q(
        q1, q2, num_patterns=num_steps, num_records_per_pattern=reps, verbose=False
    )

    # Retrieve raw IQ data from the acquired values
    I1_raw = values.rec_readout_1[0]
    Q1_raw = values.rec_readout_1[1]
    I2_raw = values.rec_readout_2[0]
    Q2_raw = values.rec_readout_2[1]

    # Build a DataFrame from the IQ data for classification
    IQ_df = pd.DataFrame({"I1": I1_raw, "Q1": Q1_raw, "I2": I2_raw, "Q2": Q2_raw})
    
    # Classify the IQ data
    classified = classifier.classify(IQ_df)
    states = classified["predicted"]
    states_reshaped = classifier.reshape_for_exp(states, reps, num_steps)
    probabilties = classifier.probabilities(states_reshaped)
    population = classifier.population(states_reshaped)

    # Build probability DataFrame
    df_prob = pd.DataFrame(
        {
            "P_f": probabilties["P_f"],
            "P_e": probabilties["P_e"],
            "P_g": probabilties["P_g"],
        }
    )

    return values, IQ_df, df_prob

def optimize_pi_pulse(
    q1: object,
    q2: object,
    general_vals_dict: dict,
    reps: int,
    initial_amp: float,
    initial_time: float,
    manifold: str = "ge",
    n_calls: int = 35,
    save_dir: str = None,
):
    """
    Optimizes the pi pulse parameters using Gaussian Process optimization.
    
    Args:
        q1, q2: Qubit objects
        general_vals_dict: General values dictionary
        reps: Number of repetitions
        initial_amp: Initial amplitude guess (from Rabi)
        initial_time: Initial time guess (from Rabi)
        manifold: 'ge' or 'ef'
        n_calls: Number of optimization steps
        save_dir: Directory to save results. If None, uses current directory.
        
    Returns:
        tuple: (optimal_amp, optimal_time, result) containing the optimized parameters
               and the full optimization result
    """
    # Define search space around initial guesses
    amp_range = 0.5  # Search ±50% of initial amplitude
    time_range = 0.5  # Search ±50% of initial time
    
    amp_vals = np.round(
        np.linspace(
            initial_amp * (1 - amp_range),
            initial_amp * (1 + amp_range),
            50
        ),
        2
    )
    time_vals = np.round(
        np.linspace(
            initial_time * (1 - time_range),
            initial_time * (1 + time_range),
            50
        ),
        1
    )
    
    space = [
        Categorical(amp_vals.tolist(), name="amp"),
        Categorical(time_vals.tolist(), name="pi_time"),
    ]

    def objective(params):
        amp, pi_time = params
        _, _, df_prob = pi_tune(
            q1,
            q2,
            general_vals_dict,
            num_steps=3,
            reps=reps,
            amp=amp,
            pi_time=pi_time,
            manifold=manifold,
        )
        
        if manifold == "ge":
            # For ge, maximize e state population difference from g state
            score = np.mean(df_prob["P_g"]) - np.mean(df_prob["P_e"])
        else:
            # For ef, maximize f state population difference from e state
            score = np.mean(df_prob["P_e"]) - np.mean(df_prob["P_f"])
        
        return score

    # Run optimization
    print(f"Starting optimization for {manifold} manifold...")
    result = gp_minimize(
        func=objective,
        dimensions=space,
        acq_func="EI",
        n_calls=n_calls,
        n_initial_points=10,
        noise="gaussian",
        random_state=42
    )

    optimal_amp = result.x[0]
    optimal_time = result.x[1]
    
    # Save results if directory provided
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
        # Save optimization results
        optimization_results = {
            f'pi_pulse_{manifold}_optimal_amplitude': optimal_amp,
            f'pi_pulse_{manifold}_optimal_time': optimal_time,
            f'pi_pulse_{manifold}_best_score': result.fun,
            f'pi_pulse_{manifold}_initial_amplitude': initial_amp,
            f'pi_pulse_{manifold}_initial_time': initial_time,
            f'pi_pulse_{manifold}_n_calls': n_calls,
            f'pi_pulse_{manifold}_n_initial_points': 10,
            f'pi_pulse_{manifold}_amp_range': amp_range,
            f'pi_pulse_{manifold}_time_range': time_range
        }
        save_optimization_results(optimization_results, '', save_dir)
        
        # Save optimization trajectory
        trajectory_file = os.path.join(save_dir, f"pi_pulse_{manifold}_trajectory.csv")
        with open(trajectory_file, 'w') as f:
            f.write("Iteration,Amplitude,Time,Score\n")
            for i, (x, y) in enumerate(zip(result.x_iters, result.func_vals)):
                f.write(f"{i},{x[0]},{x[1]},{y}\n")
    
    print(f"\nOptimization results for {manifold} manifold:")
    print(f"Optimal amplitude: {optimal_amp}")
    print(f"Optimal pi time: {optimal_time}")
    print(f"Best score: {result.fun}")
    print(f"Number of function evaluations: {result.nfev}")
    
    return optimal_amp, optimal_time, result

def swap_gate_sweep(
    q1: object,
    q2: object,
    general_vals_dict: dict,
    num_steps: int,
    reps: int,
    sweep_time: int,
    sweep_amp: list,
    sweep_freq: list
):
    """
    This function should sweep the swap gate amplitude and time for a given qubit pair
    """

def optimize_swap_gate(
    q1: object,
    q2: object,
    general_vals_dict: dict,
    reps: int,
    initial_freq: float,
    initial_amp: float,
    n_calls: int = 35,
    save_dir: str = None
):
    """
    Optimizes swap gate parameters using Gaussian Process optimization.
    Tries to maximize I1 and I2 while minimizing Q1 and Q2.
    
    Args:
        q1, q2: Qubit objects
        general_vals_dict: General values dictionary
        reps: Number of repetitions
        initial_freq: Initial frequency guess
        initial_amp: Initial amplitude guess
        n_calls: Number of optimization steps
        save_dir: Directory to save checkpoint files
        
    Returns:
        tuple: (optimal_freq, optimal_time, optimal_amp, result)
    """
    def objective(params):
        freq, time, amp = params
        
        # Run swap gate with current parameters
        pmc.parametric_coupling_time_domain(
            q1,
            q2,
            general_vals_dict,
            num_steps=3,
            ssm_para=freq,
            spec_amp=amp,
            sweep_time=time,
            phase=0,
            verbose=False,
        )
        
        # Get IQ data
        values = daq.run_daq_het_2q(
            q1, q2, num_patterns=3, num_records_per_pattern=reps, verbose=False
        )
        IQ_df = plotting.get_IQ_averages(values)
        
        # Calculate score:
        # Want to maximize I1 and I2 (negative contribution to score)
        # Want to minimize Q1 and Q2 (positive contribution to score)
        score = (IQ_df["Q1"][2] + IQ_df["Q2"][2]) - (IQ_df["I1"][2] + IQ_df["I2"][2])
        return score

    # Define the parameter space
    space = [
        (initial_freq - 0.05, initial_freq + 0.05),  # Frequency range
        (40.0, 400.0),                              # Time range in ns
        (initial_amp * 0.5, initial_amp * 1.5),      # Amplitude range
    ]

    # Setup checkpoint saving
    checkpoint_path = None
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        checkpoint_path = os.path.join(save_dir, "swap_gate_checkpoint.pkl")
        callback = [CheckpointSaver(checkpoint_path)]
    else:
        callback = None

    # Run optimization
    result = gp_minimize(
        objective,
        space,
        n_calls=n_calls,
        n_random_starts=15,
        callback=callback,
        random_state=42
    )

    optimal_freq, optimal_time, optimal_amp = result.x
    
    # Save results with clearer parameter names
    optimization_results = {
        'swap_gate_optimal_frequency': optimal_freq,
        'swap_gate_optimal_time': optimal_time,
        'swap_gate_optimal_amplitude': optimal_amp,
        'swap_gate_best_score': result.fun,
        'swap_gate_initial_frequency': initial_freq,
        'swap_gate_initial_amplitude': initial_amp,
        'swap_gate_n_calls': n_calls,
        'swap_gate_n_random_starts': 15
    }
    save_optimization_results(optimization_results, '', save_dir)

    print(f"Optimization completed:")
    print(f"Best frequency: {optimal_freq:.6f}")
    print(f"Best time: {optimal_time:.2f}")
    print(f"Best amplitude: {optimal_amp:.6f}")
    print(f"Best score: {result.fun}")

    return optimal_freq, optimal_time, optimal_amp, result