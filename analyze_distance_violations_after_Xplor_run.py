import csv
import sys 
import re
import json 

##PARAMETETRS
data_name = "runMaxAmbigFive_8.75A"
xplor_data_name = "runMaxAmbigFive_8.5A"
violations_file_name = "MaxAmbigFive_8.75A.stats"
restraints_file_name = "runMaxAmbigFive_8.75A_with_cutoff8.75A_xplor.csv"
#[%]
viol_threshold = 50
#[A]
amount_threshold = 2
###

xplor_file_name = xplor_data_name + "_noDistViol.tbl"
no_violations_csv_file_name = data_name + "_no_DistViol.csv"
marked_file_name = data_name + "_Distviol_marked.csv"
all_violations_file_name = data_name + "_all_DistViol.csv"
parameters_file_name = data_name + "_distViol_params.txt"
len_restraints_data = 0

single_bond_contacts = [["Ca", "Cb"],["C", "Ca"], ["Cb", "Cg"],["Cg", "Cd"],["Cg", "Cd2"],["Cb","Cg1"],["Cb","Cg2"],["Cg1","Cd1"],["Cd","Ce"],["Cg","Cd1"],["Cd1","Ce1"],["Cd2","Ce2"],["Ce1","Cz"],["Ce2","Cz"], ["Cd2","Ce3"],["Ce3","Cz3"], ["Cz3", "Ch2"], ["Cz2", "Ch2"], ["Ce2", "Cz2"]]
low_restraint = 2
high_restraint = 7.5

viol_file = open(violations_file_name, "r")
no_violations = []
all_violations_output = []
viol_counter = 0

def read_viol_and_rule_out():
    global viol_counter
    global len_restraints_data
    content = viol_file.read()
    start_ind = content.find("violations in NOEPot term")
    end_ind = content.find("violations in RepelPot term")
    sub_content = content[start_ind:end_ind-2]
    sub_content_lines = sub_content.splitlines()
    for j in range(6):
        sub_content_lines.pop(0)

    with open(restraints_file_name, newline = '') as f:
        reader = csv. reader(f)
        restraints_data = list(reader)
        len_restraints_data = len(restraints_data)

    for i in range(len(sub_content_lines)):
        curr_line_split = sub_content_lines[i].split()
        curr_per = float(curr_line_split[0])
        curr_amount = float(curr_line_split[1])
        curr_index = int(curr_line_split[2])
        ind_a = int(curr_line_split[6])
        name_a = curr_line_split[9]
        ind_b = int(curr_line_split[15])
        name_b = curr_line_split[18]
        curr_lst = [curr_per, curr_amount, curr_index, ind_a, name_a, ind_b, name_b]
        all_violations_output.append(curr_lst)
        if((curr_per >= viol_threshold) and (curr_amount >= amount_threshold)):
            restraints_data[curr_index].append("viol")
            viol_counter = viol_counter + 1

    with open(all_violations_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(all_violations_output)

    ##############################
    with open(marked_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(restraints_data)

    for l in range(len(restraints_data)):
        if(len(restraints_data[l]) == 6):
                if (restraints_data[l][5] == "viol"):
                    continue
        no_violations.append(restraints_data[l])

    with open(no_violations_csv_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(no_violations)

def generate_data_for_xplor(aggregated_data):
    aggregated_data_for_xplor = []

    for i in range(len(aggregated_data)):
        print(aggregated_data[i])
        peaks = aggregated_data[i][0]
        possible_assignments = []
        for k in range(len(aggregated_data[i][1])):
            curr_possibility = []
            num_str_a = aggregated_data[i][1][k][0][0]
            name_str_a = aggregated_data[i][1][k][0][1]
            shift_a = aggregated_data[i][1][k][0][2]
            total_name_a = num_str_a + name_str_a
            num_str_b = aggregated_data[i][1][k][1][0]
            name_str_b = aggregated_data[i][1][k][1][1]
            shift_b = aggregated_data[i][1][k][1][2]
            total_name_b = num_str_b + name_str_b
            curr_possibility = [total_name_a, total_name_b]
            possible_assignments.append(curr_possibility)
        aggregated_data_for_xplor.append([peaks, possible_assignments])

    return aggregated_data_for_xplor

def is_single_bond(current):
    global single_bond_contacts
    reversed_current = current[::-1]
    if current in single_bond_contacts or reversed_current in single_bond_contacts:  
        return True
    else:
        return False

def generate_xplor_restraints(xplor_data):
    xplor_file = open(xplor_file_name, "w")
    xplor_file.write("set echo=false end\nset wrnlev=0 end\n")
    for i in range(len(xplor_data)):
        lines_to_write = []
        single_bond_helper_counter = 0
        for j in range(len(xplor_data[i][1])):
            atom_a = str(xplor_data[i][1][j][0])
            atom_b = str(xplor_data[i][1][j][1])
            split_a_tmp = re.split('(\d+)',atom_a)
            split_a = [x for x in split_a_tmp if x] 
            split_b_tmp = re.split('(\d+)',atom_b)
            split_b = [x for x in split_b_tmp if x] 
            #if atom name does not include any digits
            if(len(split_a) == 2):
                resid_a = split_a[0]
                name_a = split_a[1].capitalize()
            #if atom name does include digits
            else:
                resid_a = split_a[0]
                name_a = split_a[1].capitalize() + split_a[2]
            if(len(split_b) == 2):
                resid_b = split_b[0]
                name_b = split_b[1].capitalize()
            else:
                resid_b = split_b[0]
                name_b = split_b[1].capitalize() + split_b[2]

            #distance parameters provided by the user
            a = (int(low_restraint) + float(high_restraint))/2
            b = a - int(low_restraint)
            c = float(high_restraint) - a
            #for cases of single-bond contacts
            d = a - 1
            current = [name_a, name_b]
            if(j == 0):
                #delete entire restraint if one of the options is an X-X restraint
                if((resid_a == resid_b) and (name_a == name_b)):
                    lines_to_write = []
                    break
                #handle single-bond cntact - the unambig. case 
                if((len(xplor_data[i][1]) == 1) and (abs(int(resid_a) - int(resid_b)) == 0) and (is_single_bond(current))):
                    lines_to_write = []
                    break
                else:
                    write_str = "assign (( resid " + resid_a + " and name " + name_a + " )) (( resid " + resid_b + " and name " + name_b + " )) "
                    lines_to_write.append(write_str)
            else:
                #delete entire restraint if one of the options is an X-X restraint
                if((resid_a == resid_b) and (name_a == name_b)):
                    lines_to_write = []
                    break
                #handle single-bond cntact - the ambig. case 
                if((single_bond_helper_counter == 0) and (abs(int(resid_a) - int(resid_b)) == 0) and (is_single_bond(current))):
                    single_bond_helper_counter = single_bond_helper_counter + 1
                    tmp_str = lines_to_write[0]
                    lines_to_write[0] =  tmp_str + str(a) + " " + str(d) + " " + str(c) + "\n"
                write_str = "or (( resid " + resid_a + " and name " + name_a + " )) (( resid " + resid_b  + " and name " + name_b + " ))\n"
                lines_to_write.append(write_str)
        if((single_bond_helper_counter == 0) and (len(lines_to_write) > 0)):
            tmp_str = lines_to_write[0]
            lines_to_write[0] =  tmp_str + str(a) + " " + str(b) + " " + str(c) + "\n"
        if(len(lines_to_write) > 0):
            for s in range (len(lines_to_write)):
                xplor_file.write(lines_to_write[s])

    xplor_file.write("set echo=true end\nset wrnlev=5 end\n")

def reformat_no_violations_data(unformatted):
    reformatted = []
    for i in range(len(unformatted)):
        print(unformatted[i])
        print(len(unformatted[i]))
        peaks_list_replaced = unformatted[i][0].replace("'", '"')
        peaks_lst = json.loads(peaks_list_replaced)
        assign_list_replaced = unformatted[i][1].replace("'", '"')
        assign_lst = json.loads(assign_list_replaced)
        ambig_replaced = unformatted[i][2].replace("'", '"')
        ambig = json.loads(ambig_replaced)
        ambig_type_replaced = unformatted[i][3].replace("'", '"')
        ambig_type = json.loads(ambig_type_replaced)
        sources_replaced = unformatted[i][4].replace("'", '"')
        sources = json.loads(sources_replaced)
        curr_entry = [peaks_lst, assign_lst, ambig, ambig_type, sources]
        reformatted.append(curr_entry)
    return reformatted

def write_params_file():
    param_files = open(parameters_file_name, "w")
    param_files.write(f"violations input: {violations_file_name}\n")
    param_files.write(f"total num. of viol.: {len(all_violations_output)}\n")
    param_files.write(f"restraints input: {restraints_file_name}\n")
    param_files.write(f"violation threshold: {viol_threshold}\n")
    param_files.write(f"amount threshold: {amount_threshold}\n")
    param_files.write(f"original num. of restraints: {len_restraints_data}\n")
    param_files.write(f"num. of violations by threshold: {viol_counter}\n")
    param_files.write(f"num. of restraints after filtration: {len(no_violations)}\n")

#RUN SEQUENCE

read_viol_and_rule_out()
no_violations_formatted = reformat_no_violations_data(no_violations)
aggregated_data_for_xplor = generate_data_for_xplor(no_violations_formatted)
generate_xplor_restraints(aggregated_data_for_xplor)
write_params_file()