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
    drive_amp=1,
    tomo_comp="z",
    state_comp="e",
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 50000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ef_amp = qubit_rabi.ef_amp
    ef_half_amp = qubit_rabi.ef_half_amp
    ef_half_time = qubit_rabi.ef_half_time
    ge_amp = qubit_rabi.ge_amp
    phase_offset = gen_vals["mixer_offset"]
    phase_offset_ef = gen_vals["mixer_offset_ef"]
    ROIF1 = qubit_rabi.ro_freq - qubit_rabi.RO_LO
    ROIF2 = qubit2.ro_freq - qubit2.RO_LO
    pi_ge = qubit_rabi.ge_time
    pi_ef = qubit_rabi.ef_time
    ssm_ge = qubit_rabi.ge_ssm
    ssm_ef = qubit_rabi.ef_ssm
    readout_dur = qubit_rabi.ro_dur
    mixer_offset_ge=qubit_rabi.mixer_offset_ge
    mixer_offset_ef=qubit_rabi.mixer_offset_ef
    # state prep can change 
    # state_prep_time =0# pi_ef
    buffer = 0
    if "z" == tomo_comp:
        tomo_time = 0
    elif "z-" == tomo_comp:
        tomo_time = pi_ef
    else:
        tomo_time = ef_half_time 
    
    if "e" == state_comp:
        state_prep_time = 0
    elif "f" == state_comp:
        state_prep_time = pi_ef
    else:
        state_prep_time = ef_half_time 
    
    ###########

    # comment this back in when we have found EP

    ##########
    # a_to_J = ((2 * np.pi) / (2 * (2 * qubit_rabi.ef_time * 10**-3))) / 1.5
    # # if J is in units of rad/micros
    # J_to_a = 1 / a_to_J
    # ef_amp = drive_amp_J * J_to_a
    # qubit_rabi.ef_amp = ef_amp

    # first pi_ge pulse

    pi_ge_pulse_Q = Pulse(
        start=file_length - readout_dur - swap_time - tomo_time- state_prep_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90  + mixer_offset_ge,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse_Q,
    )
    pi_ge_pulse_I = Pulse(
        start=file_length - readout_dur - swap_time - tomo_time- state_prep_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse_I,
    )
    if state_comp == "e":
        pass
    elif state_comp == "f":
        pi_ef_pulse_Q = Pulse(
            start=file_length - readout_dur - swap_time - tomo_time,
            duration=-state_prep_time,
            amplitude=ef_amp,
            ssm_freq=ssm_ef,
            phase=90+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(
            channel=4,
            sweep_name="start",
            start=0,
            stop=-sweep_time,
            initial_pulse=pi_ef_pulse_Q,
        )
        pi_ef_pulse_I = Pulse(
            start=file_length - readout_dur - swap_time - tomo_time,
            duration=-state_prep_time,
            amplitude=ef_amp,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(
            channel=1,
            sweep_name="start",
            start=0,
            stop=-sweep_time,
            initial_pulse=pi_ef_pulse_I,
        )
    elif state_comp == "x":
        half_ef_pulse_Q = Pulse(
            start=file_length - readout_dur - swap_time - tomo_time,
            duration=-state_prep_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=90+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(
            channel=4,
            sweep_name="start",
            start=0,
            stop=-sweep_time,
            initial_pulse=half_ef_pulse_Q,
        )
        half_ef_pulse_I = Pulse(
            start=file_length - readout_dur - swap_time - tomo_time,
            duration=-state_prep_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(
            channel=1,
            sweep_name="start",
            start=0,
            stop=-sweep_time,
            initial_pulse=half_ef_pulse_I,
        )
    else:
        raise ValueError("state_comp must be e, f, or x")
 
    rabi_ef_Q = Pulse(
        start=file_length - readout_dur - swap_time - tomo_time,
        duration=0,
        amplitude=drive_amp,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef_Q,
    )
    rabi_ef_I = Pulse(
        start=file_length - readout_dur - swap_time - tomo_time,
        duration=0,
        amplitude=drive_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef_I,
    )
    # tomoraphy pulses
    if tomo_comp == "z":
        pass
    elif tomo_comp == "z-":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-pi_ef,
            amplitude=ef_amp,
            ssm_freq=ssm_ef,
            phase=90+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-pi_ef,
            amplitude=ef_amp,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=tomo_pulse_I)

    elif tomo_comp == "x":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=90+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=tomo_pulse_I)
    elif tomo_comp == "x-":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=270+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=180,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=tomo_pulse_I)
    elif tomo_comp == "y":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=180+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=90,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=tomo_pulse_I)
    elif tomo_comp == "y-":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp,
            ssm_freq=ssm_ef,
            phase=270,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=tomo_pulse_I)
    else:
        raise ValueError("tomo_comp must be x,x-, y,y- or z,z-")

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=1.36,
        ssm_freq=swap_freq,
        phase=0,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)
    # swap.make()
    swap.show()
    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
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
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )

def sweep_pi_ef(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=51,
    pi_time=200,
    swap_freq=-0.21,
    swap_time=213.58765318403013,
    phase=0,
    amp_start=0.1,
    amp_stop=1.5,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ge_amp = qubit_rabi.ge_amp
    ROIF1 = qubit_rabi.ro_freq - qubit_rabi.RO_LO
    ROIF2 = qubit2.ro_freq - qubit2.RO_LO
    pi_ge = qubit_rabi.ge_time
    pi_ef = qubit_rabi.ef_time
    ssm_ge = qubit_rabi.ge_ssm
    ssm_ef = qubit_rabi.ef_ssm
    readout_dur = qubit_rabi.ro_dur
    mixer_offset_ge=qubit_rabi.mixer_offset_ge
    mixer_offset_ef=qubit_rabi.mixer_offset_ef
    
    pi_ge_pulse_Q = Pulse(
        start=file_length - readout_dur - swap_time-1*pi_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90  + mixer_offset_ge,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="none",
        initial_pulse=pi_ge_pulse_Q,
    )
    pi_ge_pulse_I = Pulse(
        start=file_length - readout_dur - swap_time-1*pi_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="none",
        initial_pulse=pi_ge_pulse_I,
    )
    # rabi_ef_Q_8 = Pulse(
    #     start=file_length - readout_dur - swap_time -7*pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_8,
    # )
    # rabi_ef_I_8 = Pulse(
    #     start=file_length - readout_dur - swap_time-7*pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_8,
    # )

    # rabi_ef_Q_7 = Pulse(
    #     start=file_length - readout_dur - swap_time-6*pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_7,
    # )
    # rabi_ef_I_7 = Pulse(
    #     start=file_length - readout_dur - swap_time-6*pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_7,
    # )

    # rabi_ef_Q_6 = Pulse(
    #     start=file_length - readout_dur - swap_time-5*pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_6,
    # )
    # rabi_ef_I_6 = Pulse(
    #     start=file_length - readout_dur - swap_time -5*pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_6,
    # )

    # rabi_ef_Q_5 = Pulse(
    #     start=file_length - readout_dur - swap_time -4*pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_5,
    # )
    # rabi_ef_I_5 = Pulse(
    #     start=file_length - readout_dur - swap_time-4*pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_5,
    # )





    # rabi_ef_Q_4 = Pulse(
    #     start=file_length - readout_dur - swap_time -3*pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_4,
    # )
    # rabi_ef_I_4 = Pulse(
    #     start=file_length - readout_dur - swap_time-3*pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_4,
    # )

    # rabi_ef_Q_3 = Pulse(
    #     start=file_length - readout_dur - swap_time-2*pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_3,
    # )
    # rabi_ef_I_3 = Pulse(
    #     start=file_length - readout_dur - swap_time-2*pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_3,
    # )

    # rabi_ef_Q_2 = Pulse(
    #     start=file_length - readout_dur - swap_time-pi_time ,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_2,
    # )
    # rabi_ef_I_2 = Pulse(
    #     start=file_length - readout_dur - swap_time -pi_time,
    #     duration=-pi_time,
    #     amplitude=0,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="amplitude",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_2,
    # )

    rabi_ef_Q_1 = Pulse(
        start=file_length - readout_dur - swap_time ,
        duration=-pi_time,
        amplitude=0,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef+phase,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="amplitude",
        start=amp_start,
        stop=amp_stop,
        initial_pulse=rabi_ef_Q_1,
    )
    rabi_ef_I_1 = Pulse(
        start=file_length - readout_dur - swap_time ,
        duration=-pi_time,
        amplitude=0,
        ssm_freq=ssm_ef,
        phase=phase,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="amplitude",
        start=amp_start,
        stop=amp_stop,
        initial_pulse=rabi_ef_I_1,
    )

    # rabi_ef_Q_1 = Pulse(
    #     start=file_length - readout_dur - swap_time ,
    #     duration=-pi_time,
    #     amplitude=1.29,
    #     ssm_freq=ssm_ef,
    #     phase=90+mixer_offset_ef+phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=4,
    #     sweep_name="none",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_Q_1,
    # )
    # rabi_ef_I_1 = Pulse(
    #     start=file_length - readout_dur - swap_time ,
    #     duration=-pi_time,
    #     amplitude=1.29,
    #     ssm_freq=ssm_ef,
    #     phase=phase,
    # )  # pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(
    #     channel=1,
    #     sweep_name="none",
    #     start=amp_start,
    #     stop=amp_stop,
    #     initial_pulse=rabi_ef_I_1,
    # )

    swap = Pulse(
        start=file_length - readout_dur ,
        duration=-swap_time,
        amplitude=1.36,
        ssm_freq=swap_freq,
        phase=0,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)
    # swap.make()
    # swap.show()
    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
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
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )

def rabi_tomo_phase(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=51,
    sweep_time=200,
    swap_freq=-0.21,
    swap_time=213.58765318403013,
    drive_amp_J=1,
    tomo_comp="z",
    phase_start=0,
    phase_stop=-10,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ef_amp = qubit_rabi.ef_amp
    ef_half_amp_x = qubit_rabi.ef_half_amp_x
    ef_half_amp_y = qubit_rabi.ef_half_amp_y
    ef_half_time = qubit_rabi.ef_half_time
    ge_amp = qubit_rabi.ge_amp
    phase_offset = gen_vals["mixer_offset"]
    phase_offset_ef = gen_vals["mixer_offset_ef"]
    ROIF1 = qubit_rabi.ro_freq - qubit_rabi.RO_LO
    ROIF2 = qubit2.ro_freq - qubit2.RO_LO
    pi_ge = qubit_rabi.ge_time
    pi_ef = qubit_rabi.ef_time
    ssm_ge = qubit_rabi.ge_ssm
    ssm_ef = qubit_rabi.ef_ssm
    readout_dur = qubit_rabi.ro_dur
    mixer_offset_ge=qubit_rabi.mixer_offset_ge
    mixer_offset_ef=qubit_rabi.mixer_offset_ef
    # state prep can change 
    state_prep_time =0# pi_ef
    buffer = 0
    if "z" == tomo_comp:
        tomo_time = 0
    else:
        tomo_time = ef_half_time  # This adds a buffer for the tomography
    
    ###########

    # comment this back in when we have found EP

    ##########
    # a_to_J = ((2 * np.pi) / (2 * (2 * qubit_rabi.ef_time * 10**-3))) / 1.5
    # # if J is in units of rad/micros
    # J_to_a = 1 / a_to_J
    # ef_amp = drive_amp_J * J_to_a
    # qubit_rabi.ef_amp = ef_amp

    # first pi_ge pulse

    pi_ge_pulse_Q = Pulse(
        start=file_length - readout_dur - swap_time - tomo_time- state_prep_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90  + mixer_offset_ge,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="none",
        initial_pulse=pi_ge_pulse_Q,
    )
    pi_ge_pulse_I = Pulse(
        start=file_length - readout_dur - swap_time - tomo_time- state_prep_time,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="none",
        initial_pulse=pi_ge_pulse_I,
    )
    
    # tomoraphy pulses
    if tomo_comp == "z":
        pass
    elif tomo_comp == "x":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp_x,
            ssm_freq=ssm_ef,
            phase=90+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="phase",start=phase_start,stop=phase_stop, initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp_x,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="phase",start=phase_start,stop=phase_stop, initial_pulse=tomo_pulse_I)
    elif tomo_comp == "y":
        tomo_pulse_Q = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp_y,
            ssm_freq=ssm_ef,
            phase=180+mixer_offset_ef,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse_Q)
        tomo_pulse_I = Pulse(
            start=file_length - readout_dur  - swap_time,
            duration=-ef_half_time,
            amplitude=ef_half_amp_y,
            ssm_freq=ssm_ef,
            phase=90,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=tomo_pulse_I)
    else:
        raise ValueError("tomo_comp must be x, y, or z")

    swap = Pulse(
        start=file_length - readout_dur - tomo_time,
        duration=-swap_time,
        amplitude=1.36,
        ssm_freq=swap_freq,
        phase=0,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)
    # swap.make()
    swap.show()
    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
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
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )