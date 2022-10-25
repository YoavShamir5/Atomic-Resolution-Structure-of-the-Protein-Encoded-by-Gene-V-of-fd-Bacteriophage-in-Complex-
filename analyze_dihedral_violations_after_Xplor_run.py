import csv 
import sys 

dataName = "runMaxAmbigFive_8.75A"
talosName = "runMaxAmbigFive_8.5A"
stats_file_name = "MaxAmbigFive_8.75A.stats"
orig_restraints_csv_file_name = "runMaxAmbigFive_9A_no_diViolations.csv"
per_threshold = 50
amount_threshold = 10
####

new_talos_file_name = "talos_" + talosName + ".tbl"
no_violations_file_name = dataName + "_no_diViolations.csv"
marked_file_name = dataName + "_dViolations_marked.csv"
viol_data_file_name = dataName + "_diViol_data.csv"
param_files_name = dataName + "_nodiViol_params.txt"
viol_data = []
no_violations = []
viol_counter = 0

def generate_new_talos_file():
    new_talos_file = open(new_talos_file_name, 'w')
    for i in range(len(no_violations)):
        curr_str_a = "assign (resid    " + no_violations[i][0] + " and name " + no_violations[i][1] + "    ) (resid    " + no_violations[i][2] + " and name " + no_violations[i][3] + "    )\n"
        new_talos_file.write(curr_str_a)
        curr_str_b = "       (resid    " + no_violations[i][4] + " and name " + no_violations[i][5] + "   ) (resid    " + no_violations[i][6] + " and name " + no_violations[i][7] + "    )    " + no_violations[i][8] + "  " + no_violations[i][9] + "   " + no_violations[i][10] + " " + no_violations[i][11] + "\n"
        new_talos_file.write(curr_str_b)

def generate_viol_data():
    stats_file = open(stats_file_name, "r")
    content = stats_file.read()
    start_ind = content.find("violations in XplorPot term CDIH")
    end_ind = content.find("violations in XplorPot term IMPR")
    sub_content = content[start_ind:end_ind-2]
    sub_content_lines = sub_content.splitlines()

    for j in range(6):
        sub_content_lines.pop(0)

    for i in range(len(sub_content_lines)):
        split_line_tmp = sub_content_lines[i].split()
        split_line = [float(split_line_tmp[0]), float(split_line_tmp[1]), int(split_line_tmp[2]), int(split_line_tmp[3]), split_line_tmp[4], split_line_tmp[5], int(split_line_tmp[6]), split_line_tmp[7], split_line_tmp[8], int(split_line_tmp[9]), split_line_tmp[10], split_line_tmp[11], int(split_line_tmp[12]), split_line_tmp[13], split_line_tmp[14]]
        viol_data.append(split_line)

def write_viol_data():
    with open(viol_data_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(viol_data)

def rule_out_viol_from_orig():
    global viol_counter
    with open(orig_restraints_csv_file_name, newline = '') as f:
        reader = csv. reader(f)
        orig_restraints_data = list(reader)

    for i in range(len(viol_data)):
        if((viol_data[i][0] >= per_threshold) and (viol_data[i][1] >= amount_threshold)):
            orig_restraints_data[viol_data[i][2]].append("viol")
            viol_counter = viol_counter + 1

    for l in range(len(orig_restraints_data)):
        if(len(orig_restraints_data[l]) == 13):
                if (orig_restraints_data[l][12] == "viol"):
                    continue
        no_violations.append(orig_restraints_data[l])

    #write marked file
    with open(marked_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(orig_restraints_data)
    #write no violations file
    with open(no_violations_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(no_violations)

    print(f"viol counter: {viol_counter}")

def write_params_file():
    param_file = open(param_files_name, "w")
    param_file.write(f"data name: {dataName} \n")
    param_file.write(f"per threhsold: {per_threshold} \n")
    param_file.write(f"amount threhsold: {amount_threshold} \n")
    param_file.write(f"stats file name: {stats_file_name} \n")
    param_file.write(f"orig_restraints_csv_file_name: {orig_restraints_csv_file_name} \n")
    param_file.write(f"viol counter: {viol_counter} \n")

#RUN PROCESS
generate_viol_data()
write_viol_data()
rule_out_viol_from_orig()
generate_new_talos_file()
write_params_file()
##