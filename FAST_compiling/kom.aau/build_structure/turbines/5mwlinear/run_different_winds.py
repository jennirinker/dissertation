import os
import numpy as np
import multiprocessing as mp

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num - 1] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

def run_fast_simulation(wind_speed):
    source_dir = os.getcwd()
    wind_str = str(round(wind_speed, 2)).replace('.', '-') + "ms" 

    wind_file = "../wind/" + wind_str + ".wnd"
    wind_file_str = '"' + wind_file + '"'
    irpt_file = "NRELOffshrBsline5MW_AeroDyn.ipt"
    fst_file = "5mw_lin"
    lin_file = "NRELOffshrBsline5MW_Linear.dat"

    run_dir = "./tmp" + wind_str
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    os.system("cp FAST " + run_dir)
    os.system("cp " + irpt_file + " " + run_dir)
    os.system("cp " + fst_file + ".fst " + run_dir)
    os.system("cp *.dat " + run_dir)
    os.system("cp -r AeroData/ " + run_dir)
    os.system("cp -r PlatformDesigns/ " + run_dir)
    os.system("cp DISCON.so " + run_dir)
    os.chdir(run_dir)

    rate_wind_speed = 12.5
    replace_line(irpt_file, 10, wind_file_str + '  WindFile    - Name of file containing wind data (quoted string)\n')
    if wind_speed >= rate_wind_speed:
        replace_line(lin_file, 6, "   3        TrimCase    - Trim case {1: find nacelle yaw, 2: find generator torque, 3: find collective blade pitch} (switch) [used only when CalcStdy=True and GenDOF=True]\n")
    else:
        replace_line(lin_file, 6, "   2        TrimCase    - Trim case {1: find nacelle yaw, 2: find generator torque, 3: find collective blade pitch} (switch) [used only when CalcStdy=True and GenDOF=True]\n")

    os.system("./FAST " + fst_file + ".fst")
    os.system("mv " + fst_file + ".lin ../outfiles/" + wind_str + ".out" )
    os.chdir(source_dir)
    os.system("rm -rf " + run_dir)
    return

def main():
    cores = mp.cpu_count()
    print cores
    p = mp.Pool(cores)
    min_wind = 15
    max_wind = 16
    wind_step = 1
    wind_speeds = np.arange(min_wind, max_wind + wind_step, wind_step)
    p.map(run_fast_simulation, wind_speeds)
    #for wind_speed in wind_speeds:
    #    run_fast_simulation(wind_speed)
    return

if __name__=='__main__':
    main()

