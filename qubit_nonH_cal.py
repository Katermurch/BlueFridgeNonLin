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

def run_rabi(
    q1,
    q2,
    general_vals_dict,
    num_steps: int,
    sweep_time: float,
    reps: int,
    manifold: str
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
    IQ_data = analysis.get_IQ_averages(values)
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
        save_dir: Directory to save checkpoint files. If None, uses current directory.
        
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

    def minimization_function(params):
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
            min_val = np.mean(df_prob["P_g"]) - np.mean(df_prob["P_e"])
        else:
            # For ef, maximize f state population difference from e state
            min_val = np.mean(df_prob["P_e"]) - np.mean(df_prob["P_f"])
        
        return min_val

    # Setup checkpoint saving with proper directory
    if save_dir is None:
        save_dir = os.getcwd()  # Use current working directory if none specified
    
    # Create save directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    
    checkpoint_file = os.path.join(save_dir, f"gp_minimize_checkpoint_{manifold}.pkl")
    checkpoint_saver = CheckpointSaver(checkpoint_file, compress=9)

    # Try to load existing checkpoint
    checkpoint_result = None
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, "rb") as f:
                checkpoint_result = pickle.load(f)
            print(f"Checkpoint loaded successfully from {checkpoint_file}")
        except Exception as e:
            print(f"Checkpoint corrupted ({e}), deleting and starting over.")
            os.remove(checkpoint_file)
            checkpoint_result = None

    # Run or continue optimization
    if checkpoint_result is not None:
        x0 = checkpoint_result.x_iters
        y0 = checkpoint_result.func_vals
        already_done = len(x0)
        remaining = n_calls - already_done
        print(f"Already completed: {already_done} calls; Remaining: {remaining}")
        if remaining <= 0:
            result = checkpoint_result
        else:
            result = gp_minimize(
                func=minimization_function,
                dimensions=space,
                acq_func="EI",
                n_calls=remaining,
                n_initial_points=0,
                noise="gaussian",
                random_state=42,
                x0=x0,
                y0=y0,
                callback=[checkpoint_saver],
            )
    else:
        print(f"Starting new optimization. Checkpoints will be saved to {checkpoint_file}")
        result = gp_minimize(
            func=minimization_function,
            dimensions=space,
            acq_func="EI",
            n_calls=n_calls,
            n_initial_points=10,
            noise="gaussian",
            random_state=42,
            callback=[checkpoint_saver],
        )

    optimal_amp = result.x[0]
    optimal_time = result.x[1]
    
    print(f"\nOptimization results for {manifold} manifold:")
    print(f"Optimal amplitude: {optimal_amp}")
    print(f"Optimal pi time: {optimal_time}")
    print(f"Best score: {result.fun}")
    
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