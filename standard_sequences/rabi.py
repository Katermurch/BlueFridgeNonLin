from generator import *
from wx_programs import *


def rabi_ge(
     qubit_rabi, qubit2, readout, num_steps=51, sweep_time=200,
    mixer_offset =0
):  # this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    #    num_steps = 51
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    #    sweep_time = 200#6000 #300 #1000 #3000
    ## channels
    pi_ge = qubit_rabi.ge_time

    ge_amp = qubit_rabi.ge_amp
    #    ssm_ge = ssm_ge_setting
    ROIF1 = qubit_rabi.ro_freq - readout.RO_LO
    ROIF2 = qubit2.ro_freq - readout.RO_LO
    readout_dur = readout.ro_dur  # 8000 #13000 #1000

    phase_offset = mixer_offset

   
    rabi_ge = Pulse(
        start=file_length - readout_dur - 100,
        duration=0,
        amplitude=ge_amp,
        ssm_freq=qubit_rabi.ge_ssm,
        phase=0,
    )
    ringupdown_seq.add_sweep(
        channel=4,
        sweep_name="width",
        start=0,
        stop=-sweep_time,
        initial_pulse=rabi_ge,
    )
    # HET
    # Q1 Readout
    # if q == 1:
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout.readout_amp_1,
        ssm_freq=ROIF1,
        phase=-file_length * ROIF1 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    # Q2 Readout
    # if q == 0:
    main_pulse = Pulse(
        start=file_length - readout_dur,
        duration=readout_dur,
        amplitude=readout.readout_amp_2,
        ssm_freq=ROIF2,
        phase=-file_length * ROIF2 * 360,
    )
    ringupdown_seq.add_sweep(channel=2, sweep_name="none", initial_pulse=main_pulse)

    ## markers
    alazar_trigger = Pulse(
        start=file_length - readout_dur - 1000, duration=50, amplitude=1
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
        # plt.figure()
        # plt.imshow(
        #     channel[:, file_length - 3000 - 300 : file_length - 3000 + 50],
        #     aspect="auto",
        # )
        # plt.show()

        # plt.figure()
        # plt.imshow(channel[:, 6000:8000], aspect="auto")
        # plt.show()

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
