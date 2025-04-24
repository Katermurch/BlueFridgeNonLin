import sys

sys.path.append(
    r"C:\Users\quantum1\Documents\Python Scripts\Important Blue Fridge Python Files\New\nonlinear_QM"
)
from classes.generator import *
from hardware_control.wx_programs import *
import numpy as np
import matplotlib.pyplot as plt

wx_addr = get_wx_address()


def T1(qubit1: object, qubit2: object, gen_vals: dict, sweep_time=100000, verbose=True):
    """
    Measures the T1 relaxation time of a qubit by applying a π pulse and
    varying the wait time before readout.

    Args:
        qubit1 (object): The qubit under test.
        qubit2 (object): The other qubit used for readout.
        gen_vals (dict): General experiment parameters, including mixer offsets.
        sweep_time (int, optional): Maximum time delay for relaxation (in ns). Defaults to 100000.
        verbose (bool, optional): If True, generates visualizations of the pulse sequence. Defaults to True.

    Returns:
        Sequence: The generated pulse sequence for T1 measurement.
    """

    num_steps = 101  # Number of delay steps

    # Calculate total file length, ensuring it's a multiple of 10,000
    totlength = sweep_time + 4000  # Extra buffer time
    file_length = 10000 * (int(np.ceil(totlength / 10000)) + 1)

    # Extract qubit parameters
    readout_dur = qubit1.ro_dur
    ringupdown_seq = Sequence(file_length, num_steps)  # Create sequence object

    ge_amp = qubit1.ge_amp
    readout_amp1 = qubit1.ro_amp
    readout_amp2 = qubit2.ro_amp
    pi_ge_time = qubit1.ge_time
    ssm_ge = qubit1.ge_ssm
    ROIF1 = qubit1.ROIF
    ROIF2 = qubit2.ROIF
    phase_offset = gen_vals["mixer_offset"]
    mixer_offset_ge=qubit1.mixer_offset_ge
    mixer_offset_ef=qubit1.mixer_offset_ef

    # Apply π pulse followed by a variable wait time (T1 measurement)
    pi_ge_Q = Pulse(
        start=file_length - readout_dur,  # Start before readout
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90  + mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_Q
    )

    pi_ge_I = Pulse(
        start=file_length - readout_dur,  # Start before readout
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_I
    )
   
    # Readout pulses for both qubits
    main_pulse_q1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_q1)

    main_pulse_q2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp2,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_q2)

    # Trigger for Alazar data acquisition
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=50, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    # Plot pulse sequence if verbose mode is enabled
    if verbose:
        channel1_ch = ringupdown_seq.channel_list[0][0]  # Qubit control channel
        channel3_ch = ringupdown_seq.channel_list[2][0]  # Additional control channel
        marker1 = ringupdown_seq.channel_list[0][2]  # Marker signal

        channel = channel1_ch + channel3_ch + marker1

        plt.figure()
        plt.imshow(
            channel[:, file_length - 3000 - 300 : file_length - 3000 + 50],
            aspect="auto",
        )
        plt.show()

        plt.figure()
        plt.imshow(channel[:, :], aspect="auto")
        plt.colorbar()
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


def T1_ef(
    qubit1: object, qubit2: object, gen_vals: dict, sweep_time=15000, verbose=True
):
    """
    Measures the T1 relaxation time of the e → f transition by applying a π_ge
    and π_ef pulse, then varying the wait time before readout.

    Args:
        qubit1 (object): The qubit under test.
        qubit2 (object): The other qubit used for readout.
        gen_vals (dict): General experiment parameters, including mixer offsets.
        sweep_time (int, optional): Maximum time delay for relaxation (in ns). Defaults to 15000.
        verbose (bool, optional): If True, generates visualizations of the pulse sequence. Defaults to True.

    Returns:
        Sequence: The generated pulse sequence for T1_ef measurement.
    """

    num_steps = 101  # Number of delay steps

    # Calculate total file length, ensuring it's a multiple of 10,000
    totlength = sweep_time + 4000  # Extra buffer time
    file_length = 10000 * (int(np.ceil(totlength / 10000)) + 1)

    # Extract qubit parameters
    readout_dur = qubit1.ro_dur
    ringupdown_seq = Sequence(file_length, num_steps)

    ge_amp = qubit1.ge_amp
    ef_amp = qubit1.ef_amp
    pi_ge_time = qubit1.ge_time
    pi_ef_time = qubit1.ef_time
    ssm_ge = qubit1.ge_ssm
    ssm_ef = qubit1.ef_ssm
    readout_amp1 = qubit1.ro_amp
    readout_amp2 = qubit2.ro_amp
    ROIF1 = qubit1.ROIF
    ROIF2 = qubit2.ROIF
    mixer_offset_ge=qubit1.mixer_offset_ge
    mixer_offset_ef=qubit1.mixer_offset_ef  # Defaults to 0 if not provided

    # Apply π_ge and π_ef pulses followed by a variable wait time (T1_ef measurement)
    pi_ge_Q = Pulse(
        start=file_length - readout_dur- pi_ef_time,  # Start before readout
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90  + mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_Q
    )

    pi_ge_I = Pulse(
        start=file_length - readout_dur- pi_ef_time,  # Start before readout
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ge_I
    )

    pi_ef_Q = Pulse(
        start=file_length - readout_dur ,
        duration=-pi_ef_time,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=90+ mixer_offset_ef,
    )
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ef_Q
    )

    pi_ef_I = Pulse(
        start=file_length - readout_dur ,
        duration=-pi_ef_time,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="start", start=0, stop=-sweep_time, initial_pulse=pi_ef_I
    )

    # Readout pulses for both qubits
    main_pulse_q1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_q1)

    main_pulse_q2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp2,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_q2)

    # Trigger for Alazar data acquisition
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=50, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    # Plot pulse sequence if verbose mode is enabled
    if verbose:
        channel1_ch = ringupdown_seq.channel_list[0][0]  # Qubit control channel
        channel3_ch = ringupdown_seq.channel_list[2][0]  # Additional control channel
        marker1 = ringupdown_seq.channel_list[0][2]  # Marker signal

        channel = channel1_ch + channel3_ch + marker1

        plt.figure()
        plt.imshow(
            channel[:, file_length - 3000 - 300 : file_length - 3000 + 50],
            aspect="auto",
        )
        plt.show()

        plt.figure()
        plt.imshow(channel[:, :], aspect="auto")
        plt.colorbar()
        plt.show()

    # Save and load the sequence for execution (if ifload=True)
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
