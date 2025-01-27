from generator import *


def pipi_pi_nopi(
    coef=0, off=0, ro_dur=8000, readout_amp=1, pi_ge=24, ssm_ge=-0.04575, ge_amp, readout_dur, sweep_time = 200
): 
    file_length = 16000
    num_steps = 3
    ringupdown_seq = Sequence(
        file_length, num_steps
    )  # this creates something called rabi_seq that is an instance of a sequence class

    rabi_ge = Pulse(
        start=file_length - readout_dur,
        duration=-pi_ge * coef,
        amplitude=ge_amp,
        ssm_freq=ssm_ge,
        phase=0,
    )  # pulse is also a class p is an instance

    ringupdown_seq.add_sweep(
        channel=3, sweep_name="none", start=0, stop=-sweep_time, initial_pulse=rabi_ge
    )
    rabi_ge.phase = 90
    ringupdown_seq.add_sweep(
        channel=4, sweep_name="none", start=0, stop=-sweep_time, initial_pulse=rabi_ge
    )

    # trigger pulse to open switch gate
    #    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
    #    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )

    # Readout
    main_pulse = Pulse(
        start=file_length - readout_dur, duration=readout_dur, amplitude=readout_amp
    )  # original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(
        channel=1, marker=2, sweep_name="none", initial_pulse=main_pulse
    )  # , marker=2

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
    #        plt.figure()
    #        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
    #        plt.show()
    #
    #        plt.figure()
    #        plt.imshow(channel[:,:], aspect='auto')
    #        plt.show()
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
def pipi_pi_nopi(coef=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,ssm_ge = -0.04575): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    num_steps = 3
    the_seq = gen.Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    ge_amp = ge_amp_setting
   
    readout_amp = ro_amp#0.5# 1
    readout_dur = ro_pulse_dur#ro_dur #8000#13000 #1000
    
 
    
#    rabi_ge = gen.Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
#    the_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    the_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    rabi_ge = gen.Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
    the_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #Readout
#    main_pulse = gen.Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    the_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = gen.Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    the_seq.load_sequence_from_disk(wx_addr, base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom    