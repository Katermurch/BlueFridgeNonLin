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
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    sweep_time=1000,
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
    mixer_offset_ge=qubit1.mixer_offset_ge
    mixer_offset_ef=qubit1.mixer_offset_ef
    buffer = 0 
    # Buffer time before readout

    # Apply the first π/2 pulse at the start of the Ramsey sequence
    t2_ge_Q = Pulse(
        start=file_length
        - readout_dur
        - buffer
        - pi_ge_time / 2,
        duration=-pi_ge_time / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=t2_ge_Q
    )

    t2_ge_I = Pulse(
        start=file_length
        - readout_dur
        - buffer
        - pi_ge_time / 2,
        duration=-pi_ge_time / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=t2_ge_I
    )


    # Apply the second π/2 pulse, with an optional phase sweep for oscillations
    final_pi2_pulse_Q = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-pi_ge_time / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )

    ringupdown_seq.add_sweep(
            channel=4, sweep_name="none", initial_pulse=final_pi2_pulse_Q)
    
    final_pi2_pulse_I = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-pi_ge_time / 2,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
            channel=1, sweep_name="none", initial_pulse=final_pi2_pulse_I)

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
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq


def ramsey_ef(
    qubit1: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    t1_time=100000,
    swap_time=1000,
    swap_freq=0,
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
    ef_amp = qubit1.ef_amp
    pi_ef_time = qubit1.ef_time
    pi_ge_time = qubit1.ge_time
    ssm_ge = qubit1.ge_ssm
    ssm_ef = qubit1.ef_ssm
    readout_amp1 = qubit1.ro_amp
    readout_dur = qubit1.ro_dur
    ef_half_amp= qubit1.ef_half_amp
    ef_half_time= qubit1.ef_half_time
    ROIF1 = qubit1.ROIF
    ROIF2 = qubit2.ROIF
    mixer_offset_ge=qubit1.mixer_offset_ge
    mixer_offset_ef=qubit1.mixer_offset_ef
    buffer = 0 
    # Buffer time before readout
    pi_ge_Q = Pulse(
        start=file_length
        - readout_dur
        - buffer
        -swap_time
        -2*ef_half_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-t1_time, initial_pulse=pi_ge_Q
    )

    pi_ge_I = Pulse(
        start=file_length
        - readout_dur
        - buffer
        -swap_time
        -2*ef_half_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-t1_time, initial_pulse=pi_ge_I
    )

    # Apply the first π/2 pulse at the start of the Ramsey sequence
    t2_ef_Q = Pulse(
        start=file_length
        - readout_dur
        - buffer
        -swap_time
        -ef_half_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-t1_time, initial_pulse=t2_ef_Q
    )

    t2_ef_I = Pulse(
        start=file_length
        - readout_dur
        - buffer
        -swap_time
        -ef_half_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-t1_time, initial_pulse=t2_ef_I
    )

    # Apply the second π/2 pulse, with an optional phase sweep for oscillations
    final_pi2_pulse_Q = Pulse(
        start=file_length - readout_dur - buffer-swap_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef,
    )
    ringupdown_seq.add_sweep(
            channel=4, sweep_name="none", initial_pulse=final_pi2_pulse_Q)
    final_pi2_pulse_I = Pulse(
        start=file_length - readout_dur - buffer-swap_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(
            channel=1, sweep_name="none", initial_pulse=final_pi2_pulse_I)
    swap = Pulse(
        start=file_length - readout_dur ,
        duration=-swap_time,
        amplitude=1.36,
        ssm_freq=swap_freq,
        phase=0,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # pi_ge_2_Q = Pulse(
    #     start=file_length
    #     - readout_dur,
    #     duration=-pi_ge_time,
    #     amplitude=ge_amp,
    #     ssm_freq=ssm_ge,
    #     phase=90+mixer_offset_ge,
    # )
    # ringupdown_seq.add_sweep(
    #     channel=4, sweep_name="start", start=0, stop=-t1_time, initial_pulse=pi_ge_2_Q
    # )

    # pi_ge_2_I = Pulse(
    #     start=file_length
    #     - readout_dur,
    #     duration=-pi_ge_time,
    #     amplitude=ge_amp,
    #     ssm_freq=ssm_ge,
    #     phase=0,
    # )
    # ringupdown_seq.add_sweep(
    #     channel=1, sweep_name="start", start=0, stop=-t1_time, initial_pulse=pi_ge_2_I
    # )
    # Readout pulse for the qubit
    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit1.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_1)

    main_pulse_2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_2)

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
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq

def ramsey_gf(
    qubit1: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    sweep_time=1000,
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
    ssm_ef = qubit1.ef_ssm
    ef_half_amp= qubit1.ef_half_amp
    ef_half_time= qubit1.ef_half_time
    ssm_ge = qubit1.ge_ssm
    readout_amp1 = qubit1.ro_amp
    readout_dur = qubit1.ro_dur
    ROIF1 = qubit1.ROIF
    mixer_offset_ge=qubit1.mixer_offset_ge
    mixer_offset_ef=qubit1.mixer_offset_ef
    buffer = 0 
    # Buffer time before readout
    pi_ge_Q_1 = Pulse(
        start=file_length
        - readout_dur-2*ef_half_time-3*pi_ge_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_Q_1
    )

    pi_ge_I_1 = Pulse(
        start=file_length
        - readout_dur-2*ef_half_time-3*pi_ge_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_I_1
    )
    # Apply the first π/2 pulse at the start of the Ramsey sequence
    t2_ge_Q = Pulse(
        start=file_length
        - readout_dur
        - buffer
        - 3*pi_ge_time
        -ef_half_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=t2_ge_Q
    )

    t2_ge_I = Pulse(
        start=file_length
        - readout_dur
        - buffer
        - 3*pi_ge_time
        -ef_half_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=t2_ge_I
    )
    pi_ge_Q_2 = Pulse(
        start=file_length
        - readout_dur-ef_half_time-2*pi_ge_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_Q_2
    )

    pi_ge_I_2 = Pulse(
        start=file_length
        - readout_dur-ef_half_time-2*pi_ge_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_I_2
    )
    pi_ge_Q_3 = Pulse(
        start=file_length
        - readout_dur-ef_half_time-pi_ge_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="none", initial_pulse=pi_ge_Q_3
    )

    pi_ge_I_3 = Pulse(
        start=file_length
        - readout_dur-ef_half_time-pi_ge_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="none", initial_pulse=pi_ge_I_3
    )
    # Apply the second π/2 pulse, with an optional phase sweep for oscillations
    final_pi2_pulse_Q = Pulse(
        start=file_length - readout_dur - buffer-pi_ge_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef,
    )

    ringupdown_seq.add_sweep(
            channel=4, sweep_name="none", initial_pulse=final_pi2_pulse_Q)
    
    final_pi2_pulse_I = Pulse(
        start=file_length - readout_dur - buffer-pi_ge_time,
        duration=-ef_half_time,
        amplitude=ef_half_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(
            channel=1, sweep_name="none", initial_pulse=final_pi2_pulse_I)
    pi_ge_Q_4 = Pulse(
        start=file_length
        - readout_dur,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="none", initial_pulse=pi_ge_Q_4
    )

    pi_ge_I_4 = Pulse(
        start=file_length
        - readout_dur,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="none", initial_pulse=pi_ge_I_4
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
        "10.225.208.207",
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
    Performs a Ramsey experiment with an additional readout pulse to measure
    quantum efficiency.

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
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq
