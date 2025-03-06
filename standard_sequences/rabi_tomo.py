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
    drive_amp_J=1,
    tomo_comp="z",
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    # qubit_rabi = qubit_rabi.copy()
    # qubit2 = qubit2.copy()
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
    tomo_time = qubit_rabi.ef_time * 1.5  # This adds a buffer for the tomography
    a_to_J = ((2 * np.pi) / (2 * (2 * qubit_rabi.ef_time * 10**-3))) / 1.5
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
        start=file_length - readout_dur - buffer - swap_time - tomo_time,
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
    # tomoraphy pulses
    if tomo_comp == "z":
        pass
    elif tomo_comp == "x":
        tomo_pulse = Pulse(
            start=file_length - readout_dur - buffer,
            duration=-qubit_rabi.ef_time,
            amplitude=ef_amp / 2,
            ssm_freq=ssm_ef,
            phase=0,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse)
    elif tomo_comp == "y":
        tomo_pulse = Pulse(
            start=file_length - readout_dur - buffer,
            duration=-qubit_rabi.ef_time,
            amplitude=ef_amp / 2,
            ssm_freq=ssm_ef,
            phase=90,
        )
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=tomo_pulse)
    else:
        raise ValueError("tomo_comp must be x, y, or z")

    swap = Pulse(
        start=file_length - readout_dur - buffer - tomo_time,
        duration=-swap_time,
        amplitude=1.23,
        ssm_freq=swap_freq,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    main_pulse_1 = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse_1)
