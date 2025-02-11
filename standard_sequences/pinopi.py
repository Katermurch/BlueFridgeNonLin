from generator import *
from wx_programs import *

wx_addr = get_wx_address()


def pi_nopi_ge(
    coef: float, offset: float, qubit: object, qubit2: object, gen_vals: dict
):
    """
    This function should run the pi nopi sequence for a qubit on the ge manifold, using the qubit's properties.

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
    rabi_ge = Pulse(
        start=file_length - readout_dur - coef * pi_ge,
        duration=pi_ge * coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge)
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
    This function should run the pi nopi sequence for a qubit on the ef manifold, using the qubit's properties.

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
    This function should run the pipi pi nopi sequence (three state) for a qubit using the qubit's properties.

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
    ### This is to get to e state###
    ge_pulse = Pulse(start=file_length - readout_dur)

    ### This is to get to f state###
    # pulse to e state
    rabi_ge = Pulse(
        start=file_length - readout_dur - (pi_ge + ef_coef * pi_ef),
        duration=pi_ge * ge_coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge)

    # pulse to f state
    rabi_ef = Pulse(
        start=file_length - readout_dur - ef_coef * pi_ef,
        duration=pi_ef * ef_coef,
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
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
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

    readout_dur = q1.ro_dur
    pi_ge_pulse = Pulse(
        start=file_length - readout_dur - copief * pi_ef - coswap * swap_time,
        duration=-pi_ge * copige,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)
    pi_ef_pulse = Pulse(
        start=file_length - readout_dur - coswap * swap_time,
        duration=-pi_ef * copief,
        amplitude=ge_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time * coswap,
        amplitude=1.7,
        ssm_freq=swap_freq,
        phase=0,
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
