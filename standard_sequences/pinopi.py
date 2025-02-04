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
