import numpy as np
import matplotlib.pyplot as plt 
import sys
from numpy.core.numeric import cross

list_file_names = ["i1", "i2", "i3", "i4",
"j1", "k1"]
spinning_speeds = [13, 13, 13, 13,
13, 13]
analysis_file_name = "analysis.txt"

analysis_handler = open(analysis_file_name, "w")
analysis_handler.write(f"files: {list_file_names}\n")
analysis_handler.write(f"speeds: {spinning_speeds}\n")

list_files = []

def generate_list_of_files():
    for i in range(len(list_file_names)):
        curr_file_name = list_file_names[i] + ".list"
        curr_file_obj = open(curr_file_name, "r")
        list_files.append(curr_file_obj)

generate_list_of_files()

for j in range(len(list_files)):
    cross_peaks = []
    peak_source_lst = []
    final_peak_list = []

    tmp_str = list_file_names[j].split(".")
    content = list_files[j].read()
    content_lines = content.splitlines()
    del content_lines[0]
    for i in range(len(content_lines)):
        if(i < 2):
            continue
        split_line = content_lines[i].split()
        w1 = split_line[1]
        w2 = split_line[2]
        cross_peaks.append([w1,w2])
        peak_source_lst.append([tmp_str[0]])

    #filtration by spinning sidebands and diagonal peaks
    x_list,y_list = map(list, zip(*cross_peaks))

    x_np_tmp = np.array(x_list)
    y_np_tmp = np.array(y_list)

    x_np = np.asarray(x_np_tmp, dtype = float)
    y_np = np.asarray(y_np_tmp, dtype = float)

    color_lst = ['blue' for i in cross_peaks]

    threhsold = 1.5

    speed = float(spinning_speeds[j])
    speed_ppm = (speed*1000)/150
    for k in range(len(cross_peaks)):
        curr_peak = cross_peaks[k]
        a = float(curr_peak[0])
        b = float(curr_peak[1])
        cond = [(((a - speed_ppm - threhsold <= x <= a - speed_ppm + threhsold) and (b - threhsold <= y <= b + threhsold)) or ((b - speed_ppm - threhsold <= x <= b - speed_ppm + threhsold) and (a - threhsold <= y <= a + threhsold))) and (abs(x-y) <= 0.3) for x,y in zip(x_np, y_np)]

        if(any(cond)):
            color_lst[k] = 'r'

        if(abs(a-b) <= 0.3):
            color_lst[k] = 'g'

    red_count = color_lst.count('r')
    green_count = color_lst.count('g')
    blue_count = color_lst.count('blue')
    print(f"total: {len(cross_peaks)}, diagonal: {green_count}, spinning sidebands: {red_count}, regular: {blue_count}")
    analysis_handler.write(f"file: {list_file_names[j]}, total: {len(cross_peaks)}, diagonal: {green_count}, spinning sidebands: {red_count}, regular: {blue_count}\n")

    for t in range(len(cross_peaks)):
        if(color_lst[t] == 'blue'):
            final_peak_list.append(cross_peaks[t])

    output_file_name = list_file_names[j] + "rss.list"
    output_handler = open(output_file_name,"w")

    output_handler.write("CC\n")
    output_handler.write("      Assignment         w1         w2\n")
    output_handler.write("\n")
    for s in range(len(final_peak_list)):
        curr_str = "              ?-?    " + str(final_peak_list[s][0]) + "     " + str(final_peak_list[s][1]) + "\n"
        output_handler.write(curr_str)

    plt.clf()
    plt.scatter(x_np, y_np, color = color_lst ,s = 0.3)
    plt.axis([0, 200, 0, 200])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.xlabel('13C [ppm]')
    plt.ylabel('13C [ppm]')
    plt.title(list_file_names[j])
    plt.savefig(list_file_names[j])