import sys

sys.path.append(
    r"C:\Users\quantum1\Documents\Python Scripts\Important Blue Fridge Python Files\New\nonlinear_QM"
)
from classes.generator import *
from hardware_control.wx_programs import *


def gef_tomo(
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
    tomo_comp="z",
    buffer=0,
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
    ef_amp = q1.ef_amp
    tomo_time = q1.ge_time * 1.5  # This adds a buffer for the tomography
    readout_dur = q1.ro_dur
    pi_ge_pulse = Pulse(
        start=file_length
        - readout_dur
        - copief * pi_ef
        - coswap * swap_time
        - tomo_time,
        duration=-pi_ge * copige,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)
    pi_ef_pulse = Pulse(
        start=file_length - readout_dur - coswap * swap_time,
        duration=-pi_ef * copief,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_pulse)
    # tomoraphy pulses
    if tomo_comp == "z":
        pass
    elif tomo_comp == "x":
        tomo_pulse = Pulse(
            start=file_length - readout_dur - buffer - swap_time,
            duration=-q1.ef_time,
            amplitude=ef_amp / 2,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse)
    elif tomo_comp == "y":
        tomo_pulse = Pulse(
            start=file_length - readout_dur - buffer - swap_time,
            duration=-q1.ef_time,
            amplitude=ef_amp / 2,
            ssm_freq=ssm_ef,
            phase=90,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse)
    else:
        raise ValueError("tomo_comp must be x, y, or z")

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time * coswap,
        amplitude=swap_amp,
        ssm_freq=swap_freq,
        phase=phase,
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
