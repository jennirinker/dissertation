
import numpy as np


def main():
    min_wind = 4.8
    max_wind = 25
    wind_step = 0.1
    wind_speeds = np.arange(min_wind, max_wind + wind_step, wind_step)
    wind_path = "./wind/"
    header = []
    header.append("!\tSteady wind file created 2013-03-21 by anb\n")
    header.append("!\tTime\tWind\tWind\tVert.\tHoriz.\tVert.\tLinV\tGust\n")
    header.append("!\t\t\tSpeed\tDir\t\tSpeed\tShear\tShear\tShear\tSpeed\n")
    times = [0.0, 999]
    wind_dir = 0
    vert_speed = 0
    horiz_shear = 0
    vert_shear = 0.2
    linv_shear = 0
    gust_speed = 0
    for wind_speed in wind_speeds:
        fil_str = str(round(wind_speed,2)).replace('.', '-') + "ms" + ".wnd"
        fil_str = wind_path + fil_str
        f = open(fil_str, 'w')
        f.writelines(header)
        for time in times:
            line = "\t" + str(time) + "\t\t" + str(wind_speed) + "\t\t" \
                    + str(wind_dir) + "\t\t" + str(vert_speed) + "\t\t" \
                    + str(horiz_shear) + "\t\t" + str(vert_shear) + "\t\t" \
                    + str(linv_shear) + "\t\t" + str(gust_speed) + "\n"
            f.write(line)


if __name__=='__main__':
    main()
