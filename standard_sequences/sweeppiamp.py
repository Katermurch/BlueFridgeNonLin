import sys

sys.path.append(
    r"C:\Users\quantum1\Documents\Python Scripts\Important Blue Fridge Python Files\New\nonlinear_QM"
)
from classes.generator import *
from hardware_control.wx_programs import *


def pi_ge_amp(
    q1: object,
    q2: object,
    gen_vals: dict,
    num_steps=3,
    amp=0.8,
    pi_ge=24,
    swap_freq=-0.21,
    swap_time=213.58765318403013,
):
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge_time = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 3
    readout_dur = q1.ro_dur

    pi_ge_1 = Pulse(
        start=file_length - readout_dur - 4 * buffer - 3 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_1)

    pi_ge_2 = Pulse(
        start=file_length - readout_dur - 3 * buffer - 2 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_2)

    pi_ge_3 = Pulse(
        start=file_length - readout_dur - 2 * buffer - 1 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_3)

    pi_ge_4 = Pulse(
        start=file_length - readout_dur - 1 * buffer - 0 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_4)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=1.36,
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


def pi_ef_amp(
    q1: object,
    q2: object,
    gen_vals: dict,
    num_steps=3,
    amp=0.8,
    pi_ge=24,
    swap_freq=-0.21,
    swap_time=213.58765318403013,
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    # (coef=0,coefpief=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,pi_ef=20,ssm_ge = -0.04575,ssm_ef=-0.04575-0.15)
    # (num_steps = 3,coef=0,coefpief=0,coefpifh=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,pi_ef=20,pi_fh=20,ssm_ge = -0.04,ssm_ef=-0.04-0.15,ssm_fh=-0.04-0.3,phase_offset=0,phase_offset_ef=0)
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge_time = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 3
    readout_dur = q1.ro_dur

    pi_ge_pulse = Pulse(
        start=file_length - 8 * buffer - 8 * pi_ge - readout_dur - swap_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_pulse)

    pi_ef_1 = Pulse(
        start=file_length - readout_dur - 8 * buffer - 7 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_1)

    pi_ef_2 = Pulse(
        start=file_length - readout_dur - 7 * buffer - 6 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_2)

    pi_ef_3 = Pulse(
        start=file_length - readout_dur - 6 * buffer - 5 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_3)

    pi_ef_4 = Pulse(
        start=file_length - readout_dur - 5 * buffer - 4 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_4)

    pi_ef_5 = Pulse(
        start=file_length - readout_dur - 4 * buffer - 3 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_5)

    pi_ef_6 = Pulse(
        start=file_length - readout_dur - 3 * buffer - 2 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_6)

    pi_ef_7 = Pulse(
        start=file_length - readout_dur - 2 * buffer - 1 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_7)

    pi_ef_8 = Pulse(
        start=file_length - readout_dur - 1 * buffer - 0 * pi_ge - swap_time,
        duration=-pi_ge,
        amplitude=amp,
        ssm_freq=ssm_ef,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_8)

    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=1.36,
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


def pi_ge_sweep(
    q1: object,
    q2: object,
    gen_vals: dict,
    num_steps=3,
):
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge_time = q1.ge_time
    ef_amp = q1.ef_amp
    buffer = 3
    readout_dur = q1.ro_dur
    mixer_offset_ge = q1.mixer_offset_ge
    mixer_offset_ef = q1.mixer_offset_ef

    pi_ge_Q = Pulse(
        start=file_length - readout_dur,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90 + mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_Q)

    pi_ge_I = Pulse(
        start=file_length - readout_dur,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pi_ge_I)

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


def pi_ef_sweep(
    q1: object,
    q2: object,
    gen_vals: dict,
    num_steps=3,
):
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class
    ssm_ge = q1.ge_ssm
    ssm_ef = q1.ef_ssm
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    ge_amp = q1.ge_amp
    pi_ef = q1.ef_time
    pi_ge_time = q1.ge_time
    pi_ef_time = q1.ef_time
    ef_amp = q1.ef_amp
    buffer = 3
    readout_dur = q1.ro_dur
    mixer_offset_ge = q1.mixer_offset_ge
    mixer_offset_ef = q1.mixer_offset_ef

    pi_ge_Q = Pulse(
        start=file_length - readout_dur - pi_ef_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90 + mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_Q)

    pi_ge_I = Pulse(
        start=file_length - readout_dur - pi_ef_time,
        duration=-pi_ge_time,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pi_ge_I)
    pi_ef_Q = Pulse(
        start=file_length - readout_dur,
        duration=-pi_ef_time,
        amplitude=ef_amp,
        ssm_freq=ssm_ge,
        phase=90 + mixer_offset_ef,
    )
    ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ef_Q)

    pi_ef_I = Pulse(
        start=file_length - readout_dur,
        duration=-pi_ef_time,
        amplitude=ef_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pi_ef_I)
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


def pi_ge_gaussian_tune_iq(
    q1,
    q2,
    gen_vals,
    num_steps=3,
    amp=0.8,
    pi_ge_time=24,
    swap_freq=-0.21,
    swap_time=213.58765318403013,
    file_length=16000,
    verbose=True,
):
    ringupdown_seq = Sequence(file_length, num_steps)

    ssm_ge = q1.ge_ssm
    mixer_offset_ge = q1.mixer_offset_ge
    readout_amp_1 = q1.ro_amp
    ROIF1 = q1.ROIF
    readout_amp_2 = q2.ro_amp
    ROIF2 = q2.ROIF
    readout_dur = q1.ro_dur
    buffer = 3

    # Loop to add 4 pi pulses (I and Q channels)
    for i in range(4):
        start_time = (
            file_length
            - readout_dur
            - (4 - i) * buffer
            - (3 - i) * pi_ge_time
            - swap_time
        )
        pi_ge_I = Pulse(
            start=start_time,
            duration=-pi_ge_time,
            amplitude=amp,
            ssm_freq=ssm_ge,
            phase=0,
            gaussian_bool=True,
        )
        pi_ge_Q = Pulse(
            start=start_time,
            duration=-pi_ge_time,
            amplitude=amp,
            ssm_freq=ssm_ge,
            phase=90 + mixer_offset_ge,
            gaussian_bool=True,
        )
        ringupdown_seq.add_sweep(channel=1, sweep_name="none", initial_pulse=pi_ge_I)
        ringupdown_seq.add_sweep(channel=4, sweep_name="none", initial_pulse=pi_ge_Q)

    # Swap pulse
    swap = Pulse(
        start=file_length - readout_dur,
        duration=-swap_time,
        amplitude=1.36,
        ssm_freq=swap_freq,
        phase=0,
    )
    ringupdown_seq.add_sweep(channel=3, sweep_name="none", initial_pulse=swap)

    # Readout pulses for both qubits
    for amp, ROIF in [(readout_amp_1, ROIF1), (readout_amp_2, ROIF2)]:
        readout = Pulse(
            start=file_length - readout_dur,
            duration=readout_dur,
            amplitude=amp,
            ssm_freq=ROIF,
            phase=-file_length * ROIF * 360,
        )
        ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=readout)

    # Alazar trigger
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    # Plot if verbose
    if verbose:
        ch1, ch3, m1 = (
            ringupdown_seq.channel_list[0][0],
            ringupdown_seq.channel_list[2][0],
            ringupdown_seq.channel_list[0][2],
        )
        full = ch1 + ch3 + m1
        plt.figure()
        plt.imshow(
            full[
                :, file_length - readout_dur - 1000 - 4000 : file_length - readout_dur
            ],
            aspect="auto",
        )
        plt.show()

    # Save and load the sequence
    write_dir = (
        r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    )
    ringupdown_seq.write_sequence_to_disk(
        "foo", write_dir, use_range_01=False, num_offset=0, write_binary=True
    )
    ringupdown_seq.load_sequence_from_disk(
        "128.252.134.31", "foo", write_dir, num_offset=0, ch_amp=[1, 1, 1, 1]
    )

    return ringupdown_seq
