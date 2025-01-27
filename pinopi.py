from generator import *
from wx_programs import *
wx_addr =  get_wx_address()


def pipi_pi_nopi(coef=0,off=0,pi_ge=24, ge_amp = 1, ssm_ge = -0.04575,readout_dur = 1):
    file_length = 16000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200
    
 
    rabi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
    the_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)

    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] 
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1

    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    the_seq.load_sequence_from_disk(wx_addr, base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom    