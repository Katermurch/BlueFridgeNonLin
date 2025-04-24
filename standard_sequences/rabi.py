from classes.generator import *
from hardware_control.wx_programs import *

ge_amp_setting = 1.2  # 1.2#1.9 for 1in driving
ssm_ge_setting = (
    -0.15
)  # -.162#-0.17 #0.0875#-0.17 #0.0875 #0.170#0.1675 # 0.150 #0.1779 + 0.004#0.1751
pi_ge_time_setting = 25  # 33#24#76#17 #176 #102   #50#18 #64
pi2_ge_time_setting = 16  # 12 #88 #51

ssm_ef_setting = -0.0675  # -0.20555
pi_ef_time_setting = 48
pi2_ef_time_setting = 24

ro_pulse_dur = 4000  # 5000
readout_amp_1 = 0.6  # .5#1.2#.7#1
readout_amp_2 = 0.4  # .55#.7#.4#1.0#0.3077#0.769#0.134*0.6

mixer_offset = 0  # 7.7
mixer_offset_ef = 20


def rabi_ge(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    sweep_time=200,
):
    """
    This function should run the rabi ge sequence for a qubit, using the
    qubit's properties.

    Args:
        qubit_rabi (_type_): this is the qubit you are performing the rabi on
        qubit2 (_type_): this qubit exists for homodyne readout
        gen_vals (dict): a dictionary of general values for readout and hardware control
        num_steps (int, optional): _description_. Defaults to 101.
        sweep_time (int, optional): _description_. Defaults to 200.
    """

    # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    ## channels
    pi_ge = qubit_rabi.ge_time
    ge_amp = qubit_rabi.ge_amp
    ROIF1 = qubit_rabi.ROIF
    ROIF2 = qubit2.ROIF
    readout_dur = qubit_rabi.ro_dur
    mixer_offset_ge=qubit_rabi.mixer_offset_ge
    mixer_offset_ef=qubit_rabi.mixer_offset_ef
    phase_offset = gen_vals["mixer_offset"]
    buffer = 0
    rabi_ge_x = Pulse(
        start=file_length - readout_dur  ,
        duration=0,
        amplitude=ge_amp,
        ssm_freq=qubit_rabi.ge_ssm,
        phase=90+mixer_offset_ge,
    )
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ge_x,
    )
    rabi_ge_y = Pulse(
        start=file_length - readout_dur ,
        duration=0,
        amplitude=ge_amp,
        ssm_freq=qubit_rabi.ge_ssm,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ge_y,
    )
    # HET readout
    # Q1 Readout
    main_pulse = Pulse(
        start=file_length  - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length  - readout_dur - 1000, duration=50, amplitude=1
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
        base_name="foo",  # rabi_oop_test
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
        ch_amp=gen_vals["wx_amps"],
    )


def rabi_ef(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=51,
    sweep_time=200,
):
    """
    This function should run the rabi ef sequence for a qubit, using the
    qubit's properties.

    Args:
        qubit_rabi (object): this is the qubit you are performing the rabi on
        qubit2 (object): this qubit exists for homodyne readout
        readout (object): this object chooses ther readout parameters
        gen_vals (dict): a dictionary of general values for readout and hardware control
        num_steps (int, optional): _description_. Defaults to 51.
        sweep_time (int, optional): _description_. Defaults to 200.
    """

    # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
    #    num_steps = 101
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    #    sweep_time = 200#20000# #300 #1000 #3000
    ## channels

    ## channels
    pi_ge = qubit_rabi.ge_time
    ef_amp = qubit_rabi.ef_amp
    ge_amp = qubit_rabi.ge_amp
    ssm_ge = qubit_rabi.ge_ssm
    ssm_ef = qubit_rabi.ef_ssm
    ROIF1 = qubit_rabi.ROIF
    ROIF2 = qubit2.ROIF
    mixer_offset_ge=qubit_rabi.mixer_offset_ge
    mixer_offset_ef=qubit_rabi.mixer_offset_ef
    readout_dur = qubit_rabi.ro_dur
    phase_offset = gen_vals["mixer_offset"]
    buffer = 0

    # Pi ge amp
    pi_ge_pulse_x = Pulse(
        start=file_length
        - readout_dur,  # the 100 is the match the buffer below
        duration=-pi_ge,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=90+mixer_offset_ge,
    
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="start",
        start=0,
        stop=-sweep_time,
        initial_pulse=pi_ge_pulse_x,
    )
    pi_ge_pulse_y = Pulse(
        start=file_length
        - readout_dur,  # the 100 is the match the buffer below
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
        initial_pulse=pi_ge_pulse_y,
    )
    # drive rabi e-f
    rabi_ef_x = Pulse(
        start=file_length
        - readout_dur
        ,  # buffer here to make sure rabi doesnt bleed into readout
        duration=0,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=90+mixer_offset_ef,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef_x,
    )
    rabi_ef_y = Pulse(
        start=file_length
        - readout_dur
        ,  # buffer here to make sure rabi doesnt bleed into readout
        duration=0,
        amplitude=ef_amp,
        ssm_freq=ssm_ef,
        phase=0,
    )  # pulse is also a class p is an instance
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ef_y,
    )
    # Rabi Qubit Readout

    main_pulse = Pulse(
        start=file_length - readout_dur ,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Other qubit Readout

    main_pulse = Pulse(
        start=file_length - readout_dur ,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur  - 1000, duration=1000, amplitude=1
    )
    ringupdown_seq.add_sweep(
        channel=3, marker=1, sweep_name="none", initial_pulse=alazar_trigger
    )

    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0]
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


def amp_rabi_ge(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    amp_range=1,
):
    """
    This function should run the rabi ge sequence for a qubit, using the
    qubit's properties.

    Args:
        qubit_rabi (_type_): this is the qubit you are performing the rabi on
        qubit2 (_type_): this qubit exists for homodyne readout
        gen_vals (dict): a dictionary of general values for readout and hardware control
        num_steps (int, optional): _description_. Defaults to 101.
        sweep_time (int, optional): _description_. Defaults to 200.
    """

    # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    ## channels
    pi_ge = qubit_rabi.ge_time
    ge_amp = qubit_rabi.ge_amp
    ROIF1 = qubit_rabi.ROIF
    ROIF2 = qubit2.ROIF
    readout_dur = qubit_rabi.ro_dur
    phase_offset = gen_vals["mixer_offset"]
    buffer = 500
    rabi_ge = Pulse(
        start=file_length - readout_dur - buffer,
        duration=0,
        amplitude=ge_amp,
        ssm_freq=qubit_rabi.ge_ssm,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="amplitude",
        start=0,
        stop=amp_range,
        initial_pulse=rabi_ge,
    )
    # HET readout
    # Q1 Readout
    main_pulse = Pulse(
        start=file_length - buffer - readout_dur,
        duration=readout_dur,
        amplitude=qubit_rabi.ro_amp,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    main_pulse = Pulse(
        start=file_length - buffer - readout_dur,
        duration=readout_dur,
        amplitude=qubit2.ro_amp,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - buffer - readout_dur - 1000, duration=50, amplitude=1
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
        base_name="foo",  # rabi_oop_test
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
        ch_amp=gen_vals["wx_amps"],
    )

def test_rabi_ge(
    qubit_rabi: object,
    qubit2: object,
    gen_vals: dict,
    num_steps=101,
    phase_offset=0,
):
    """
    This function should run the rabi ge sequence for a qubit, using the
    qubit's properties.

    Args:
        qubit_rabi (_type_): this is the qubit you are performing the rabi on
        qubit2 (_type_): this qubit exists for homodyne readout
        gen_vals (dict): a dictionary of general values for readout and hardware control
        num_steps (int, optional): _description_. Defaults to 101.
        sweep_time (int, optional): _description_. Defaults to 200.
    """

    # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    ## channels
    ef_amp = qubit_rabi.ef_amp
    pi_ge = qubit_rabi.ge_time
    ge_amp = qubit_rabi.ge_amp
    ROIF1 = qubit_rabi.ROIF
    ROIF2 = qubit2.ROIF
    readout_dur = qubit_rabi.ro_dur
    
    buffer =0
    rabi_ge_Q = Pulse(
        start=file_length - 16000  ,
        duration=16000,
        amplitude=ef_amp,
        ssm_freq=qubit_rabi.ef_ssm,
        phase=90+phase_offset,
    )
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="none",
        
        initial_pulse=rabi_ge_Q,
    )
    rabi_ge_I = Pulse(
        start=file_length - 16000,
        duration=16000,
        amplitude=ef_amp,
        ssm_freq=qubit_rabi.ef_ssm,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=1,
        sweep_name="none",
        
        initial_pulse=rabi_ge_I,
    )
    

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