from generator import *
from wx_programs import *

wx_addr = get_wx_address()


def pi_nopi_ge(
    coef: float, offset: float, qubit: object, qubit2: object, readout: object, gen_vals: dict
):
    """
    This function should run the pipi pi nopi sequence for a qubit, using the qubit's properties.

    Args:
        coef (float): coefficient for the pi pulse
        offeset (float): offset for the pi pulse
        qubit (class): takes a qubit class instance and runs on this qubit
        readout (class): takes a readout class instance and runs on this readout
        gen_vals (dict): takes a general_vals class instance and runs on these values
        save_dir (str): directory where the sequence will be saved
    """
    readout_dur = readout.ro_dur
    pi_ge = qubit.ge_time
    ssm_ge = qubit.ge_ssm
    ROIF1 = qubit.ro_freq - readout.RO_LO
    ROIF2 = qubit2.ro_freq - readout.RO_LO
    ge_amp = qubit.ge_amp
    file_length = 16000
    num_steps = 3
    the_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    # channel 4 is for qubit control
    rabi_ge = Pulse(
        start=file_length - readout_dur -coef * pi_ge,
        duration=pi_ge * coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    the_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=rabi_ge)
     #Readout
    # HET readout
    # Q1 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout.readout_amp_1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    the_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout.readout_amp_2,
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
