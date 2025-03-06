import sys

sys.path.append(
    r"C:\Users\quantum1\Documents\Python Scripts\Important Blue Fridge Python Files\New\nonlinear_QM"
)
from classes.generator import *
from hardware_control.wx_programs import *


def rabi_ef_swap_tomo(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=51,
    sweep_time=200,
    swap_freq=-0.21,
    swap_time=213.58765318403013,
    drive_amp_J = 1,
    tomo_comp = 'z'
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 3000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    qubit_rabi = qubit_rabi.deepcopy()
    qubit2 = qubit2.deepcopy()
    ef_amp = qubit_rabi.ef_amp
    ge_amp = qubit_rabi.ge_amp
    phase_offset = gen_vals["mixer_offset"]
    phase_offset_ef = gen_vals["mixer_offset_ef"]
    ROIF1 = qubit_rabi.ro_freq - qubit_rabi.RO_LO
    ROIF2 = qubit2.ro_freq - qubit2.RO_LO
    pi_ge = qubit_rabi.ge_time
    ssm_ge = qubit_rabi.ge_ssm
    ssm_ef = qubit_rabi.ef_ssm
    readout_dur = qubit_rabi.ro_dur
    buffer = 0
    tomo_time = qubit_rabi.ef_time * 1.5  #This adds a buffer for the tomography
    a_to_J = ((2 * np.pi) / (2 * (2 * qubit_rabi.ef_time * 10**-3)))/1.5
          # if J is in units of rad/micros
    J_to_a = 1 / a_to_J
    ef_amp = drive_amp_J * J_to_a
    qubit_rabi.ef_amp = ef_amp

    # first pi_ge pulse
    
    pi_ge_pulse = Pulse(
        start=file_length - readout_dur - buffer - swap_time - tomo_time,
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
        start=file_length - readout_dur - buffer - swap_time- tomo_time,
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
        start=file_length - readout_dur -buffer - tomo_time,
        duration=-swap_time,
        amplitude=1.23,
        ssm_freq=swap_freq,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)


    #tomoraphy pulses
    if tomo_comp == 'z':
        pass
    elif tomo_comp == 'x':
        tomo_pulse = Pulse(
            start=file_length - readout_dur - buffer,
            duration=-qubit_rabi.ef_time,
            amplitude=ef_amp/2,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse)
    elif tomo_comp == 'y':
        tomo_pulse = Pulse(
            start=file_length - readout_dur - buffer,
            duration=-qubit_rabi.ef_time,
            amplitude=ef_amp/2,
            ssm_freq=ssm_ef,
            phase=90,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse)


    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_1)

    # Q2 Readout
    # if q == 0:
    main_pulse_2 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_2)
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


def full_tomo_ef(
    q1: object,
    q2: object,
    gen_vals: dict,
    off=0,
    drive_amp = 1,
    num_steps=51,
    sweep_time=200,
    tomo=1,
    amp_tomo=1,
    ph=0,
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
        amplitude=readout_amp,
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
