from classes.generator import *
from hardware_control.wx_programs import *


def parametric_coupling(
    qubit1: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    ssm_ge=-0.2,
    pi_ge=20,
    ssm_start=-0.15,
    ssm_stop=-0.25,
    spec_amp=0.5,
    ROIF1=0,
    ROIF2=0,
    verbose=True,
):
    """
    Performs a parametric coupling experiment between qubits by applying a drive
    pulse and sweeping a sideband modulation frequency.

    Args:
        qubit1 (object): The primary qubit being driven.
        qubit2 (object): The secondary qubit used for readout.
        gen_vals (dict): General experiment parameters, including mixer offsets.
        num_steps (int, optional): Number of steps in the parametric sweep. Defaults to 101.
        ssm_ge (float, optional): Single sideband modulation frequency for π pulse. Defaults to -0.2.
        pi_ge (int, optional): Duration of the π pulse in ns. Defaults to 20.
        ssm_start (float, optional): Start frequency for sideband sweep. Defaults to -0.15.
        ssm_stop (float, optional): Stop frequency for sideband sweep. Defaults to -0.25.
        spec_amp (float, optional): Amplitude of parametric drive. Defaults to 0.5.
        ROIF1 (float, optional): Readout intermediate frequency for qubit 1. Defaults to 0.
        ROIF2 (float, optional): Readout intermediate frequency for qubit 2. Defaults to 0.
        verbose (bool, optional): If True, generates visualizations of the pulse sequence. Defaults to True.

    Returns:
        Sequence: The generated pulse sequence for parametric coupling measurement.
    """

    sweep_time = 1000  # Duration of parametric drive
    totlength = sweep_time + 4000  # Extra buffer time
    file_length = 10000 * (int(np.ceil(totlength / 10000)) + 1)

    ringupdown_seq = Sequence(file_length, num_steps)  # Initialize sequence

    # Extract qubit parameters
    ge_amp = qubit1.ge_amp
    readout_amp1 = qubit1.ro_amp
    readout_amp2 = qubit2.ro_amp
    readout_dur = qubit1.ro_dur
    phase_offset = gen_vals["mixer_offset"]

    # Apply π pulse on the selected qubit (q=0 or q=1)
    pi_ge_pulse = Pulse(
        start=file_length - readout_dur - sweep_time - 10,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    drive_channel = 4 if qubit1.qubit_id == "q1" else 3
    ringupdown_seq.add_sweep(
        channel=drive_channel, sweep_name="none", initial_pulse=pi_ge_pulse
    )

    # Apply parametric drive with frequency sweep
    parametric_drive = Pulse(
        start=file_length - readout_dur - 10,
        duration=-sweep_time,
        amplitude=spec_amp,
        ssm_freq=0,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=2,
        sweep_name="ssm_freq",
        start=ssm_start,
        stop=ssm_stop,
        initial_pulse=parametric_drive,
    )

    # Readout pulses for both qubits
    readout_pulse_q1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=1.3 * readout_amp1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="none", initial_pulse=readout_pulse_q1
    )

    readout_pulse_q2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=1.3 * readout_amp2,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(
        channel=1, sweep_name="none", initial_pulse=readout_pulse_q2
    )

    # Trigger for Alazar data acquisition
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    # Plot pulse sequence if verbose mode is enabled
    channel1_ch = ringupdown_seq.channel_list[0][0]  # Qubit control channel
    channel3_ch = ringupdown_seq.channel_list[2][0]  # Additional control channel
    marker1 = ringupdown_seq.channel_list[0][2]  # Marker signal

    channel = channel1_ch + channel3_ch + marker1
    if verbose:

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
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
    )

    return ringupdown_seq
