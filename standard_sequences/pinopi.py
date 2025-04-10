from classes.generator import *
from hardware_control.wx_programs import *

wx_addr = get_wx_address()


def pi_nopi_ge(
    coef: float, offset: float, qubit: object, qubit2: object, gen_vals: dict
):
    """
    This function should run the pi nopi sequence for a qubit on the ge
    manifold, using the qubit's properties.

    Args:
        coef (float): coefficient for the pi pulse
        offeset (float): offset for the pi pulse
        qubit (class): takes a qubit class instance and runs on this qubit
        qubit2 (class): takes a second qubit class instance this is just for readout
        gen_vals (dict): takes a general_vals class instance and runs on these values
        save_dir (str): directory where the sequence will be saved
    """
    readout_dur = qubit.ro_dur
    pi_ge = qubit.ge_time
    ssm_ge = qubit.ge_ssm
    ROIF1 = qubit.ROIF
    ROIF2 = qubit2.ROIF
    ge_amp = qubit.ge_amp
    file_length = 16000
    num_steps = 3
    the_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    # channel 4 is for qubit control
    rabi_ge_x = Pulse(
        start=file_length - readout_dur,
        duration=pi_ge * coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90,
    )
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge_x)
    rabi_ge_y = Pulse(
        start=file_length - readout_dur,
        duration=pi_ge * coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=rabi_ge_y)
    # Readout
    # HET readout
    # Q1 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # channel 3 marker 1 is for pc trigger
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    the_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    the_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=offset,
        write_binary=True,
    )
    the_seq.load_sequence_from_disk(
        wx_addr, base_name="foo", file_path=write_dir, num_offset=0, ch_amp=[1, 1, 1, 1]
    )


def pi_nopi_ef(
    coef: float, offset: float, qubit: object, qubit2: object, gen_vals: dict
):
    """
    This function should run the pi nopi sequence for a qubit on the ef
    manifold, using the qubit's properties.

    Args:
        coef (float): coefficient for the pi pulse
        offeset (float): offset for the pi pulse
        qubit (class): takes a qubit class instance and runs on this qubit
        qubit2 (class): takes a second qubit class instance this is just for readout
        gen_vals (dict): takes a general_vals class instance and runs on these values
        save_dir (str): directory where the sequence will be saved
    """
    readout_dur = qubit.ro_dur
    pi_ef = qubit.ef_time
    pi_ge = qubit.ge_time
    ssm_ge = qubit.ge_ssm
    ssm_ef = qubit.ef_ssm
    ROIF1 = qubit.ROIF
    ROIF2 = qubit2.ROIF
    ge_amp = qubit.ge_amp
    ef_amp = qubit.ef_amp
    file_length = 16000
    num_steps = 3
    the_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    # we need to write it seperately here since the ge T1 time is short, so we want to make sure the readout is close to the end of the pi pulse for the e state
    if 0 == coef:
        # pulse to e state
        rabi_ge = Pulse(
            start=file_length - readout_dur - (pi_ge),
            duration=pi_ge * coef,
            amplitude=ge_amp,
            ssm_freq=ssm_ge,
            phase=0,
        )
        the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge)

        # HET readout
        # Q1 Readout
        main_pulse = Pulse(
            start=file_length - readout_dur,
            duration=readout_dur,
            amplitude=qubit.ro_amp,
            ssm_freq=ROIF1,
            phase=-file_length * ROIF1 * 360,
        )
        the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

        # Q2 Readout
        main_pulse = Pulse(
            start=file_length - readout_dur,
            duration=readout_dur,
            amplitude=qubit2.ro_amp,
            ssm_freq=ROIF2,
            phase=-file_length * ROIF2 * 360,
        )
        the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

        # channel 3 marker 1 is for pc trigger
        alazar_trigger = Pulse(
            start=file_length - readout_dur - 1000, duration=1000, amplitude=1
        )
        the_seq.add_sweep(
            channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
        )

    else:

        # pulse to e state
        rabi_ge = Pulse(
            start=file_length - readout_dur - (pi_ge + coef * pi_ef),
            duration=pi_ge * coef,
            amplitude=ge_amp,
            ssm_freq=ssm_ge,
            phase=0,
        )
        the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge)

        # pulse to f state
        rabi_ef = Pulse(
            start=file_length - readout_dur - coef * pi_ef,
            duration=pi_ef * coef,
            amplitude=ef_amp,
            ssm_freq=ssm_ef,
            phase=0,
        )
        the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ef)
        # Readout
        # HET readout
        # Q1 Readout
        main_pulse = Pulse(
            start=file_length - readout_dur,
            duration=readout_dur,
            amplitude=qubit.ro_amp,
            ssm_freq=ROIF1,
            phase=-file_length * ROIF1 * 360,
        )
        the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

        # Q2 Readout
        main_pulse = Pulse(
            start=file_length - readout_dur,
            duration=readout_dur,
            amplitude=qubit2.ro_amp,
            ssm_freq=ROIF2,
            phase=-file_length * ROIF2 * 360,
        )
        the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

        # channel 3 marker 1 is for pc trigger
        alazar_trigger = Pulse(
            start=file_length - readout_dur - 1000, duration=1000, amplitude=1
        )
        the_seq.add_sweep(
            channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
        )

    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    the_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=offset,
        write_binary=True,
    )
    the_seq.load_sequence_from_disk(
        wx_addr, base_name="foo", file_path=write_dir, num_offset=0, ch_amp=[1, 1, 1, 1]
    )


def pipi_pi_nopi(
    ge_coef: float,
    ef_coef: float,
    offset: float,
    qubit: object,
    qubit2: object,
    gen_vals: dict,
):
    """
    This function should run the pipi pi nopi sequence (three state) for a
    qubit using the qubit's properties.

    Args:
        coef (float): coefficient for the pi pulse
        offeset (float): offset for the pi pulse
        qubit (class): takes a qubit class instance and runs on this qubit
        qubit2 (class): takes a second qubit class instance this is just for readout
        gen_vals (dict): takes a general_vals class instance and runs on these values
        save_dir (str): directory where the sequence will be saved
    """
    readout_dur = qubit.ro_dur
    pi_ef = qubit.ef_time
    pi_ge = qubit.ge_time
    ssm_ge = qubit.ge_ssm
    ssm_ef = qubit.ef_ssm
    ROIF1 = qubit.ROIF
    ROIF2 = qubit2.ROIF
    ge_amp = qubit.ge_amp
    ef_amp = qubit.ef_amp
    file_length = 16000
    num_steps = 3
    the_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    # pulse to e state
    rabi_ge_x = Pulse(
        start=file_length - readout_dur - (pi_ge + ef_coef * pi_ef),
        duration=pi_ge * ge_coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90,
    )
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge_x)
    rabi_ge_y = Pulse(
        start=file_length - readout_dur - (pi_ge + ef_coef * pi_ef),
        duration=pi_ge * ge_coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=rabi_ge_y)
    # pulse to f state
    rabi_ef_x = Pulse(
        start=file_length - readout_dur - ef_coef * pi_ef,
        duration=pi_ef * ef_coef,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=90,
    )
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ef_x)
    rabi_ef_y = Pulse(
        start=file_length - readout_dur - ef_coef * pi_ef,
        duration=pi_ef * ef_coef,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    the_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=rabi_ef_y)
    # Readout
    # HET readout
    # Q1 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # channel 3 marker 1 is for pc trigger
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    the_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    the_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=offset,
        write_binary=True,
    )
    the_seq.load_sequence_from_disk(
        wx_addr, base_name="foo", file_path=write_dir, num_offset=0, ch_amp=[1, 1, 1, 1]
    )


def pi_nopi_swap(
    q1: object,
    q2: object,
    gen_vals_dict: dict,
    copief,
    coswap,
    copige,
    num_steps=3,
    swap_freq=0,
    swap_time=60,
    phase=0,
    swap_amp=0,
    off=0
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    num_steps = num_steps
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    off = 0
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 0
    readout_dur = q1.ro_dur
    pi_ge_pulse_x = Pulse(
        start=file_length
        - readout_dur
        - copief * pi_ef
        - coswap * swap_time
        - 2 * buffer,
        duration=-pi_ge * copige,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse_x)

    pi_ge_pulse_y = Pulse(
        start=file_length
        - readout_dur
        - copief * pi_ef
        - coswap * swap_time
        - 2 * buffer,
        duration=-pi_ge * copige,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pi_ge_pulse_y)
    pi_ef_pulse_x = Pulse(
        start=file_length - readout_dur - coswap * swap_time - buffer,
        duration=-pi_ef * copief,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=90,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse_x)
    pi_ef_pulse_y = Pulse(
        start=file_length - readout_dur - coswap * swap_time - buffer,
        duration=-pi_ef * copief,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pi_ef_pulse_y)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time * coswap,
        amplitude=swap_amp,
        ssm_freq=swap_freq,
        phase=phase,
        gaussian_bool= False
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # Readout

    # HET
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
    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    #
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=off,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        '10.225.208.207',
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )

def pi_nopi_z(
    q1: object,
    q2: object,
    gen_vals_dict: dict,
    copief,
    num_steps=3,
    swap_freq=0,
    swap_time=60,
    phase=0,
    swap_amp=0,
    pi_time=10,
    pi_amp=1,
    pi_phase=0,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    num_steps = num_steps
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    off = 0
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 0
    readout_dur = q1.ro_dur
    pi_ge_pulse = Pulse(
        start=file_length
        - readout_dur
        - copief * pi_ef
        -  swap_time
        - 2 * buffer,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)
    pi_ef_pulse = Pulse(
        start=file_length - readout_dur -  swap_time - buffer,
        duration=-pi_ef * copief,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=swap_amp,
        ssm_freq=swap_freq,
        phase=phase,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # Readout

    # HET
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
    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    #
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=off,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )


def pi_nopi_x(
    q1: object,
    q2: object,
    gen_vals_dict: dict,
    copief,
    num_steps=3,
    swap_freq=0,
    swap_time=60,
    phase=0,
    swap_amp=0,
    pi_time=10,
    pi_amp=1,
    pi_phase=0,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    num_steps = num_steps
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    off = 0
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 0
    
    readout_dur = q1.ro_dur
    pi_ge_pulse = Pulse(
        start=file_length
        - readout_dur
        - pi_time
        -  swap_time
        - 2 * buffer,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)
    pi_ef_pulse = Pulse(
        start=file_length - readout_dur -  swap_time - buffer,
        duration=-pi_time,
        amplitude=pi_amp,
        ssm_freq=ssm_ef,
        phase=copief*180 + pi_phase,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=swap_amp,
        ssm_freq=swap_freq,
        phase=phase,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # Readout

    # HET
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
    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    #
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=off,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )

def pi_nopi_y(
    q1: object,
    q2: object,
    gen_vals_dict: dict,
    copief,
    num_steps=3,
    swap_freq=0,
    swap_time=60,
    phase=0,
    swap_amp=0,
    pi_time=10,
    pi_amp=1,
    pi_phase=0,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    num_steps = num_steps
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    off = 0
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 0
    readout_dur = q1.ro_dur
    pi_ge_pulse = Pulse(
        start=file_length
        - readout_dur
        -  pi_time
        -  swap_time
        - 2 * buffer,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)
    pi_ef_pulse = Pulse(
        start=file_length - readout_dur -  swap_time - buffer,
        duration=-pi_time,
        amplitude=pi_amp,
        ssm_freq=ssm_ef,
        phase=copief*180+90+ pi_phase,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=swap_amp,
        ssm_freq=swap_freq,
        phase=phase,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # Readout

    # HET
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
    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    #
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=off,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        "128.252.134.31",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )

def pi_nopi_sweep_phase(
    q1: object,
    q2: object,
    gen_vals_dict: dict,
    copief,
    num_steps=3,
    swap_freq=0,
    swap_time=60,
    phase=0,
    swap_amp=0,
    pi_time=10,
    pi_amp=1,
    pi_phase=0,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    num_steps = num_steps
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    off = 0
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 0
    
    readout_dur = q1.ro_dur
    pi_ge_pulse = Pulse(
        start=file_length
        - readout_dur
        - 2*pi_time
        -  swap_time
        - 2 * buffer,
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)

    pi_ef_pulse = Pulse(
        start=file_length - readout_dur -pi_time -  swap_time - buffer,
        duration=-pi_time,
        amplitude=1.46,
        ssm_freq=ssm_ef,
        phase=153,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse)

    pi_ef_pulse_2 = Pulse(
        start=file_length - readout_dur -  swap_time - buffer,
        duration=-pi_time,
        amplitude=1.49,
        ssm_freq=ssm_ef,
        phase=153 + pi_phase,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse_2)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=swap_amp,
        ssm_freq=swap_freq,
        phase=phase,
        gaussian_bool=False,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # Readout

    # HET
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
    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    #
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        base_name="foo",
        file_path=write_dir,
        use_range_01=False,
        num_offset=off,
        write_binary=True,
    )
    ringupdown_seq.load_sequence_from_disk(
        "10.225.208.207",
        base_name="foo",
        file_path=write_dir,
        num_offset=0,
        ch_amp=[1, 1, 1, 1],
    )
def resonator_sweep_q2(q1,q2,num_steps=3,ro_amplitude=1,co_pi=0,ssm_ge = -0.2,pi_ge =20,ROIF1=0,ROIF2=0):
    file_length = 16000 #ns
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 0
    
    readout_dur = 5000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class

    #2000 #13000 #1000

    ef_amp=1.5
    # pi_ge_pulse = Pulse(start=file_length-readout_dur-10, duration=-pi_ge*co_pi, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
    # pi_ef_pulse = Pulse(start=file_length-readout_dur-10, duration=-pi_ef*co_pief, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    # ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ef_pulse)

    
  
    
#READOUT:::HET    

#    Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=ro_amplitude,ssm_freq=ROIF2, phase=-file_length*ROIF2*360)
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)



    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk(
        '10.225.208.207', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1, 1, 1, 1])
    return ringupdown_seq