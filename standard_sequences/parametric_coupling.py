from generator import *
from wx_programs import *


def parametric_coupling_time_domain(
    qubit: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    ssm_para=0,
    spec_amp=0.5,
    sweep_time=0,
):
    
    pi_ge = qubit.ge_time
    ssm_ge = qubit.ge_ssm
    ge_amp = qubit.ge_amp
    readout_dur = qubit.ro_dur
    phase_offset = gen_vals['mixer_offset']
    phase_offset_ef = gen_vals['mixer_offset_ef']
    ROIF1 = qubit.ROIF
    ROIF2 = qubit2.ROIF

    totlength = sweep_time + 4000
    file_length = 10000 * (int(np.ceil(totlength / 10000)) + 1)
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class



    pi_ge_pulse = Pulse(
        start=file_length - readout_dur,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse,
    )

    rabi_ge = Pulse(
        start=file_length - readout_dur,
        duration=0,
        amplitude=spec_amp,
        ssm_freq=ssm_para,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=3,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ge,
    )
    # HET
    # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 )
    # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

    # Q1 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp_1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout_amp_2,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    #    main_pulse.phase = 90
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

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
