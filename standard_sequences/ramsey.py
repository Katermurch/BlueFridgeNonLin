import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append(
    r"C:\Users\quantum1\Documents\Python Scripts\Important Blue Fridge Python Files\New\nonlinear_QM"
)
from classes.generator import *
from hardware_control.wx_programs import *


def ramsey(
    qubit1: object,
    gen_vals: dict,
    num_steps=101,
    t1_time=100000,
    pi_echo_coef=0,
    osc=0,
    verbose=True,
):
    """
    Performs a Ramsey experiment to measure qubit dephasing time (T2*).

    Args:
        qubit1 (object): The qubit under test.
        gen_vals (dict): General experiment parameters, including mixer offsets.
        num_steps (int, optional): Number of delay steps for T2* measurement. Defaults to 101.
        t1_time (int, optional): Maximum wait time before applying second π/2 pulse (in ns). Defaults to 100000.
        pi_echo_coef (float, optional): Echo pulse coefficient. Defaults to 0 (no echo).
        osc (int, optional): Oscillation frequency for the phase of the final π/2 pulse. Defaults to 0.
        verbose (bool, optional): If True, generates visualizations of the pulse sequence. Defaults to True.

    Returns:
        Sequence: The generated pulse sequence for the Ramsey experiment.
    """

    file_length = 16000  # Fixed file length
    ringupdown_seq = Sequence(file_length, num_steps)  # Initialize sequence

    # Extract qubit parameters
    ge_amp = qubit1.ge_amp
    pi_ge_time = qubit1.ge_time
    ssm_ge = qubit1.ge_ssm
    readout_amp1 = qubit1.ro_amp
    readout_dur = qubit1.ro_dur
    ROIF1 = qubit1.ROIF
    buffer = 0  # Buffer time before readout

    # Apply the first π/2 pulse at the start of the Ramsey sequence
    t2_ge = Pulse(
        start=file_length
        - readout_dur
        - buffer
        - pi_ge_time / 2
        - pi_echo_coef * pi_ge_time,
        duration=-pi_ge_time / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-t1_time, initial_pulse=t2_ge
    )

    # Apply an optional echo π pulse halfway through the wait time
    if pi_echo_coef != 0:
        t2_ge = Pulse(
            start=file_length - readout_dur - buffer - pi_ge_time / 2,
            duration=-pi_ge_time * pi_echo_coef,
            amplitude=ge_amp,
            ssm_freq=ssm_ge,
            phase=0,
        )
        ringupdown_seq.add_sweep(
            channel=4,
            sweep_name="start",
            start=0,
            stop=-t1_time / 2,
            initial_pulse=t2_ge,
        )

    # Apply the second π/2 pulse, with an optional phase sweep for oscillations
    final_pi2_pulse = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-pi_ge_time / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )

    if osc != 0:
        ringupdown_seq.add_sweep(
            channel=4,
            sweep_name="phase",
            start=0,
            stop=osc * 360,
            initial_pulse=final_pi2_pulse,
        )
    else:
        ringupdown_seq.add_sweep(
            channel=4, sweep_name="none", initial_pulse=final_pi2_pulse
        )

    # Readout pulse for the qubit
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Trigger for Alazar data acquisition
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=50, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )
    channel1_ch = ringupdown_seq.channel_list[0][0]  # Qubit control channel
    channel3_ch = ringupdown_seq.channel_list[2][0]  # Additional control channel
    marker1 = ringupdown_seq.channel_list[0][2]  # Marker signal
    # Plot pulse sequence if verbose mode is enabled
    if verbose:

        channel = channel1_ch + channel3_ch + marker1

        plt.figure()
        plt.imshow(
            channel[
                :, file_length - readout_dur - 1000 - 4000 : file_length - readout_dur
            ],
            aspect="auto",
        )
        plt.show()

    # Save and load the sequence for execution
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=0,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq


def ramsey_quantum_efficiency(
    qubit1: object,
    gen_vals: dict,
    num_steps=101,
    ssm_ge=-0.15,
    pi_ge=20,
    ROIF=0,
    RO_ram_amp=0.0,
    verbose=True,
):
    """
    Performs a Ramsey experiment with an additional readout pulse to measure quantum efficiency.

    Args:
        qubit1 (object): The qubit under test.
        gen_vals (dict): General experiment parameters, including mixer offsets.
        num_steps (int, optional): Number of steps in the Ramsey sequence. Defaults to 101.
        ssm_ge (float, optional): Single sideband modulation frequency for the π/2 pulses. Defaults to -0.15.
        pi_ge (int, optional): Duration of the π/2 pulses in ns. Defaults to 20.
        ROIF (float, optional): Readout intermediate frequency. Defaults to 0.
        RO_ram_amp (float, optional): Amplitude of additional readout pulse for efficiency measurement. Defaults to 0.0.
        verbose (bool, optional): If True, generates visualizations of the pulse sequence. Defaults to True.

    Returns:
        Sequence: The generated pulse sequence for Ramsey quantum efficiency measurement.
    """

    file_length = 16000  # Fixed file length
    ringupdown_seq = Sequence(file_length, num_steps)  # Initialize sequence

    # Extract qubit parameters
    ge_amp = qubit1.ge_amp
    readout_amp1 = qubit1.ro_amp
    readout_dur = 5000  # Adjusted readout duration for quantum efficiency measurement
    buffer = 1000  # Buffer time before readout

    # Apply the first π/2 pulse
    first_pi2_pulse = Pulse(
        start=file_length - 2 * readout_dur - buffer - pi_ge / 2 - 100,
        duration=-pi_ge / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="none", initial_pulse=first_pi2_pulse
    )

    # Apply additional readout pulse for quantum efficiency measurement
    readout_pulse_efficiency = Pulse(
        start=file_length - readout_dur - 50 - pi_ge / 2 - buffer,
        duration=-readout_dur,
        amplitude=RO_ram_amp,
        ssm_freq=ROIF,
        phase=-file_length * ROIF * 360,
    )
    ringupdown_seq.add_sweep(
        channel=2, sweep_name="none", initial_pulse=readout_pulse_efficiency
    )

    # Apply the second π/2 pulse with phase sweep for quantum efficiency measurement
    final_pi2_pulse = Pulse(
        start=file_length - readout_dur - 50,
        duration=-pi_ge / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="phase", start=0, stop=360, initial_pulse=final_pi2_pulse
    )

    # Readout pulse for the qubit
    readout_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp1,
        ssm_freq=ROIF,
        phase=-file_length * ROIF * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=readout_pulse)

    # Trigger for Alazar data acquisition
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=50, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )
    channel1_ch = ringupdown_seq.channel_list[0][0]  # Qubit control channel
    channel3_ch = ringupdown_seq.channel_list[2][0]  # Additional control channel
    marker1 = ringupdown_seq.channel_list[0][2]  # Marker signal

    channel = channel1_ch + channel3_ch + marker1

    # Plot pulse sequence if verbose mode is enabled
    if verbose:

        plt.figure()
        plt.imshow(
            channel[
                :, file_length - readout_dur - 1000 - 4000 : file_length - readout_dur
            ],
            aspect="auto",
        )
        plt.show()

    # Save and load the sequence for execution
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=0,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq
