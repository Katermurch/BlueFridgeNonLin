import sys

main_directory = r"C:\Users\quantum1\Documents"

hardware_paths = {
    "instrument_path": r"C:\Users\quantum1\Documents\Python\instr\analyzer",
    "analyzer_path": r"C:\Users\quantum1\Documents\Python\instr\python_interface\python_without_WX2184C",
}


def add_paths_to_sys():
    for path in hardware_paths.values():
        if path not in sys.path:
            sys.path.append(path)


bnc_address = {
    "target_bnc_black": "GPIB0::19::INSTR",
    "big_agilent": "GPIB0::30::INSTR",
    "agilent_function_generator": "GPIB0::10::INSTR",
    "target_bnc_6": "USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR",
}

add_paths_to_sys()
