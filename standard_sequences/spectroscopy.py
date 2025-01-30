from generator import *
from wx_programs import *

wx_addr = get_wx_address()


def spectroscopy_ef(
    qubit: object,
    readout: object,
    gen_vals: dict,
    ssm_start=-0.15,
    ssm_stop=-0.25,
    spec_amp=0.5,
    ROIF1=0,
    ROIF2=0,
):
    """This function does the ef spectroscopy for the qubit by running over single sideband frequencies

    Args:
        qubit (object): qubit that we perform the spectroscopy on
        readout (object): readout parameters
        gen_vals (dict): other parameters that may be necessary
        ssm_start (float, optional): sweep start frequency Defaults to -0.15.
        ssm_stop (float, optional): sweep stop frequency Defaults to -0.25.
        spec_amp (float, optional): amplitude of pulse for spectroscopy Defaults to 0.5.
        ROIF1 (int, optional): actual sent frequency based on readout frequency and qubit readout for q1 Defaults to 0.
        ROIF2 (int, optional):  actual sent frequency based on readout frequency and qubit readout for q2. Defaults to 0.

    """

    sweep_time = 200
    file_length = 16000
    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    ## channels
    ge_amp = qubit.ge_amp
    ssm_ge = qubit.ge_ssm
    pi_ge = qubit.ge_time
    readout_amp = readout.readout_amp_1
    readout_dur = readout.ro_dur

    pi_ge_pulse = Pulse(
        start=file_length - readout_dur - sweep_time - 10,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)

    rabi_ge = Pulse(
        start=file_length - readout_dur - 10,
        duration=-sweep_time,
        amplitude=spec_amp,
        ssm_freq=0,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="ssm_freq",
        start=ssm_start,
        stop=ssm_stop,
        initial_pulse=rabi_ge,
    )

    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=50, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][
            0
        ]  # [channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]

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
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )
    return ringupdown_seq
