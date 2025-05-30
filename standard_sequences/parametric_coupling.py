from classes.generator import *
from hardware_control.wx_programs import *


def parametric_coupling_time_domain(
    qubit1: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    ssm_para=0,
    spec_amp=0.5,
    sweep_time=0,
    phase=0,
    verbose=True,
):
    """
    Performs a time-domain parametric coupling experiment between qubits by
    applying a drive pulse and varying its duration.

    Args:
        qubit1 (object): The primary qubit being driven.
        qubit2 (object): The secondary qubit used for readout.
        gen_vals (dict): General experiment parameters, including mixer offsets.
        num_steps (int, optional): Number of steps in the parametric time sweep. Defaults to 101.
        ssm_para (float, optional): Sideband modulation frequency for parametric drive. Defaults to 0.
        spec_amp (float, optional): Amplitude of parametric drive. Defaults to 0.5.
        sweep_time (int, optional): Maximum duration for parametric drive pulse. Defaults to 0.
        verbose (bool, optional): If True, generates visualizations of the pulse sequence. Defaults to True.

    Returns: 
        Sequence: The generated pulse sequence for parametric coupling time-domain measurement.
    """

    totlength = sweep_time + 4000  # Extra buffer time
    file_length = 10000 * (int(np.ceil(totlength / 10000)) + 1)

    ringupdown_seq = Sequence(file_length, num_steps)  # Initialize sequence

    # Extract qubit parameters
    ge_amp = qubit1.ge_amp
    readout_amp1 = qubit1.ro_amp
    readout_amp2 = qubit2.ro_amp
    readout_dur = qubit1.ro_dur
    pi_ge = qubit1.ge_time
    ssm_ge = qubit1.ge_ssm
    ROIF1 = qubit1.ROIF
    ROIF2 = qubit2.ROIF
    mixer_offset_ge= qubit1.mixer_offset_ge
    buffer = 0

    # Apply π pulse on the selected qubit
    pi_ge_pulse_Q = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse_Q,
    )
    pi_ge_pulse_I = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse_I,
    )

    # Apply parametric drive with duration sweep
    parametric_drive = Pulse(
        start=file_length - readout_dur - buffer,
        duration=0,  # Initially zero, swept in time
        amplitude=spec_amp,
        ssm_freq=ssm_para,
        phase=phase,
        gaussian_bool= False
    )
    ringupdown_seq.add_sweep(
        channel=3,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=parametric_drive,
    )

    # Readout pulses for both qubits
    readout_pulse_q1 = Pulse(
        start=file_length - readout_dur - buffer,
        duration=readout_dur,
        amplitude=readout_amp1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(
        channel=2, sweep_name="none", initial_pulse=readout_pulse_q1
    )

    readout_pulse_q2 = Pulse(
        start=file_length - readout_dur - buffer,
        duration=readout_dur,
        amplitude=readout_amp2,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(
        channel=2, sweep_name="none", initial_pulse=readout_pulse_q2
    )

    # Trigger for Alazar data acquisition
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000 - buffer, duration=1000, amplitude=1
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
        plt.imshow(channel[:, 6000:8000], aspect="auto")
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
        '10.225.208.207',
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq
