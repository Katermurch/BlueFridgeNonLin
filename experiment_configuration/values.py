readout_dict = {
    "RO_LO": 6.6247,
    "RO_LO_pwr": 16,
    "ro_dur": 5000,
}
swap_gate = {
    "swap_freq": 0,
    "swap_amp": 0,
    "swap_time": 0,
}
q1_dict = {
    "qubit_id": "q1",
    "ro_freq": 6.72739,
    "ge_ssm": -0.1144,
    "ef_ssm": -0.25684,#-0.2568,#-0.2581,#
    "ge_time": 77,
    "ef_time": 44,
    "ef_half_time": 25,
    "ef_half_amp": 1.3,
    "ge_amp": 1.01,#.5
    "ef_amp":1.5,
    "IQ_angle":60,
    "ro_amp": .25,#0.2#.35,
    "qubit_thr": [-10000, -600],
    "mixer_offset_ge":2.5,
    "mixer_offset_ef":6,
}
q2_dict = {
    "qubit_id": "q2",
    "ro_freq": 6.65554,
    "ge_ssm": -0.154,
    "ef_ssm": -0.2962,
    "ge_time": 74,
    "ge_amp": 0.4,
    "ef_amp": 1,
    "IQ_angle": 25,
    "ro_amp": .45,#5,
    "qubit_thr": [-10000, 1900],
}
q3_dict = {
    "qubit_id": "q3",
    "ro_freq": 6.65839,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
    "RO_LO": 6.4804,
}
q4_dict = {
    "qubit_id": "q4",
    "ro_freq": 6.5113,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
}
q5_dict = {
    "qubit_id": "q5",
    "ro_freq": 6.44436,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
}

general_vals_dict = {
    "mixer_offset": 0,
    "mixer_offset_ef": 20,
    "wx_amps": [1,1,1.7,1.015],  # maximum 1.9
    "coupler_off_value": 0.7,
    "wx_offs": [-0.036, 0, -0.08, -0.098],
    "qubit_bnc": 4.6,
    "TWPA_freq": 4.5,#4.45,
    "TWPA_pwr": -5.6#-5.5,
}
bnc_address = {
    "target_bnc_black": "GPIB0::19::INSTR",
    "big_agilent": "GPIB0::11::INSTR",
    "agilent_function_generator": "GPIB0::30::INSTR",
    "target_bnc_6": "USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR",
    "wx_address": "10.225.208.204",
}
