import sys

sys.path.append(
    r"C:\Users\quantum1\Documents\Python Scripts\Important Blue Fridge Python Files\New\nonlinear_QM"
)
from generator import *
from wx_programs import *


def full_tomo(
    q1: object,
    q2: object,
    gen_vals: dict,
    off=0,
    num_steps=51,
    sweep_time=200,
    tomo=1,
    amp_tomo=1,
    ph=90,
):

    file_length = 16000
    the_seq = Sequence(file_length, num_steps)

    ge_amp = q1.ge_amp

    readout_amp = q1.ro_amp  # 0.5# 1
    readout_dur = q1.ro_dur  # 8000 #13000 #1000
    ssm_ge = q1.ge_ssm
    pi_ge = q1.ge_time
    ROIF1 = q1.ROIF
    ROIF2 = q2.ROIF

    phase_offset = gen_vals["mixer_offset"]
    buffer = 5

    rabi_ge = Pulse(
        start=file_length - readout_dur - tomo * pi_ge / 2 - buffer * 2,
        duration=0,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(
        channel=1,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ge,
    )
    rabi_ge.phase = 90 + phase_offset
    the_seq.add_sweep(
        channel=2,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ge,
    )

    # pulse for x,y,z tomography
    pulse_tomo = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-tomo * pi_ge / 2,
        amplitude=amp_tomo,
        ssm_freq=ssm_ge,
        phase=0 + ph,
    )
    the_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pulse_tomo)
    pulse_tomo.phase = 90 + ph + phase_offset
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=pulse_tomo)

    # Readout
    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=q1.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_1)

    # Q2 Readout
    # if q == 0:
    main_pulse_2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=q2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_2)
    #    main_pulse.phase = 90
    # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 )
    # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    the_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][
            0
        ]  # [channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]

        channel = channel1_ch + channel3_ch + marker1

    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    the_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=0,
        write_binary=True,
    )
    the_seq.load_sequence_from_disk(
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )


def full_tomo_ef(
    q1: object,
    q2: object,
    gen_vals: dict,
    off=0,
    num_steps=51,
    sweep_time=200,
    tomo=1,
    amp_tomo=1,
    ph=90,
):
    file_length = 16000
    the_seq = Sequence(file_length, num_steps)

    ge_amp = q1.ge_amp
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    pi_ge = q1.ge_time
    pi_ef = q1.ef_time
    ROIF2 = q2.ROIF
    ROIF1 = q1.ROIF

    readout_amp = 1  # 0.5# 1
    readout_dur = q1.ro_dur  # 8000 #13000 #1000

    phase_offset = gen_vals["mixer_offset"]
    phase_offset_ef = gen_vals["mixer_offset_ef"]
    buffer = 5

    pi_pulse_ge = Pulse(
        start=file_length - readout_dur - tomo * pi_ef / 2 - buffer * 3,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(
        channel=3,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_pulse_ge,
    )
    pi_pulse_ge.phase = 90 + phase_offset
    the_seq.add_sweep(
        channel=4,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_pulse_ge,
    )

    # rabi_ef drive
    rabi_ef = Pulse(
        start=file_length - readout_dur - pi_ef / 2 * tomo - buffer * 2,
        duration=0,
        amplitude=ge_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    the_seq.add_sweep(
        channel=3,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef,
    )
    rabi_ef.phase = 90 + phase_offset_ef
    the_seq.add_sweep(
        channel=4,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef,
    )

    # pulse for x,y,z tomography
    pulse_tomo = Pulse(
        start=file_length - readout_dur - buffer,
        duration=-tomo * pi_ef / 2,
        amplitude=amp_tomo,
        ssm_freq=ssm_ef + 0.001,
        phase=0 + ph,
    )
    the_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=pulse_tomo)
    pulse_tomo.phase = 90 + ph + phase_offset_ef
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pulse_tomo)

    # Readout
    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=q1.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_1)

    # Q2 Readout
    # if q == 0:
    main_pulse_2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=q2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_2)
    #    main_pulse.phase = 90
    # main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 )
    # ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    the_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]

        channel = channel1_ch + channel3_ch + marker1

    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    the_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=0,
        write_binary=True,
    )
    the_seq.load_sequence_from_disk(
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )
