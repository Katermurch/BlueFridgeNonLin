from generator import *
from wx_programs import *


def rabi_ef_swap(
    qubit_rabi: object,
    qubit2: object, 
    gen_vals: dict,
    num_steps=51,
    sweep_time=200,
    swap_freq=0.1,
    swap_time=100,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class


    ge_amp = qubit_rabi.ge_amp
    phase_offset = gen_vals['mixer_offset']
    phase_offset_ef = gen_vals['mixer_offset_ef']
    readout_amp = qubit_rabi.ro_amp
    ROIF1 = qubit_rabi.ro_freq - qubit_rabi.RO_LO
    ROIF2 = qubit2.ro_freq - qubit2.RO_LO
    pi_ge = qubit_rabi.ge_time
    ssm_ge = qubit_rabi.ge_ssm
    ssm_ef = qubit_rabi.ef_ssm
    ef_amp = qubit_rabi.ef_amp
    readout_dur = qubit_rabi.ro_dur
    buffer = 0

    # first pi_ge pulse
    pi_ge_pulse = Pulse(
        start=file_length - readout_dur - pi_ge - buffer - swap_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse,
    )
    # drive rabi e-f
    rabi_ef = Pulse(
        start=file_length - readout_dur - pi_ge - buffer - swap_time,
        duration=0,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef,
    )
    
    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=1.7,
        ssm_freq=swap_freq,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)
        
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    # if q == 0:
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)
    #    main_pulse.phase = 90
    # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 )
    # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
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
