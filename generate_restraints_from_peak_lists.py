import sys
import numpy as np
import re
import csv
from collections import defaultdict
import json

curr_scheme = "fully"

error = 0.3
threshold = 0.4
max_ambig = 5
low_restraint = 2
high_restraint = 7.5
file_name = "runE_fully"

list_file_names = ["a1rss.list", "a2rss.list", "a3rss.list", "a4rss.list",
"b1rss.list", "b2rss.list", "b3rss.list", "b4rss.list"]

cross_peaks = []
cs_lst = []
possible_assign_per_peak_lst = []
possible_assign_per_peak_lst_rand_del = []
cross_peaks_with_neighbors = []
csv_xplor_data = []
peak_source_lst = []

analysis_file_path = "analysis_" + file_name + ".csv"
analysis_aggregated_file_path = "analysis_aggregated_" + file_name + ".csv"

single_bond_contacts = [["Ca", "Cb"],["C", "Ca"], ["Cb", "Cg"],["Cg", "Cd"],["Cg", "Cd2"],["Cb","Cg1"],["Cb","Cg2"],["Cg1","Cd1"],["Cd","Ce"],["Cg","Cd1"],["Cd1","Ce1"],["Cd2","Ce2"],["Ce1","Cz"],["Ce2","Cz"], ["Cd2","Ce3"],["Ce3","Cz3"], ["Cz3", "Ch2"], ["Cz2", "Ch2"], ["Ce2", "Cz2"]]

onethreeglycerol_dict = {"G": {"C": 1, "CA": 0}, 
                "A":{"C": 1, "CA": 0, "CB": 1},
                "S":{"C": 1, "CA": 0, "CB": 1}, 
                "C":{"C": 1, "CA": 0, "CB": 1},
                "P":{"C": 2/5,"CA": 4/5, "CB": 3/5, "CG": 1,"CD": 0}, 
                "R":{"C": 2/5, "CA": 4/5, "CB": 3/5, "CG": 1,"CD": 0, "CZ": 0}, #CZ prob. data not available, chosen to be prob=0
                "D":{"C": 3/6, "CA": 4/6, "CB": 5/6,"CG": 2/6}, 
                "N":{"C": 3/6, "CA": 4/6, "CB": 5/6,"CG": 2/6}, 
                "M":{"C": 3/6,"CA": 4/6,"CB": 5/6,"CG": 2/6,"CE": 1}, 
                "T":{"C": 3/6, "CA": 4/6, "CB": 5/6, "CG2": 2/6}, 
                "I":{"C": 3/6,"CA": 4/6,"CB": 0, "CG1": 5/6, "CG2": 1, "CD1":2/6}, #good,check gamma numbering
                "K":{"C": 9/12,"CA": 5/12,"CB": 11/12, "CG": 4/12,"CD": 10/12, "CE": 5/12}, #good, check D
                "F":{"C": 1, "CA": 0, "CB": 1, "CG": 0, "CD1": 1, "CD2": 1, "CE1": 0, "CE2": 1, "CZ":1}, #good, check CE1/CE2 difference
                "Y":{"C": 1, "CA": 0, "CB": 1, "CG": 0, "CD1": 1, "CD2": 1, "CE1": 0, "CE2": 1, "CZ":1},
                #"W":{},
                "H":{"C": 1,"CA": 0,"CB": 1, "CG": 1/2, "CD2": 1/2, "CE1": 0}, #check nomenclature
                "E":{"C": 2/5,"CA": 4/5,"CB": 3/5,"CG": 1, "CD": 0},
                "Q":{"C": 2/5,"CA": 4/5,"CB": 3/5,"CG": 1, "CD": 0}, 
                "V":{"C": 1, "CA": 0, "CB": 0, "CG1": 1, "CG2": 1},
                "L":{"C": 0, "CA": 1, "CB": 0, "CG": 0, "CD1": 1, "CD2": 1} 
                }

list_files = []

for i in range(len(list_file_names)):
    curr_file_obj = open(list_file_names[i], "r")
    list_files.append(curr_file_obj)

shifts_file = open("shifts_data_new.shifts", "r")
sequence_file = open("seq_file.txt", "r")
xplor_file_name = "xplor_" + file_name + ".tbl"
xplor_file = open(xplor_file_name, "w")
sequence = sequence_file.read()

def is_single_bond(current):
    global single_bond_contacts
    reversed_current = current[::-1]
    if current in single_bond_contacts or reversed_current in single_bond_contacts:  
        return True
    else:
        return False

def update_distance_counters(cross_peak_name):
    dist_a = int(cross_peak_name[0][0])
    dist_b = int(cross_peak_name[1][0])
    if(dist_a == dist_b):
        #intra
        return 1
    elif(abs(dist_a - dist_b) == 1):
        #sequential
        return 2
    elif(2 <= abs(dist_a - dist_b) <= 4):
        #medium
        return 3
    else:
        #long
        return 4

#generate cs_lst
def generate_cs_lst():
    content = shifts_file.read()
    content_lines = content.splitlines()
    del content_lines[0]
    for i in range(len(content_lines)):
        split_line = content_lines[i].split()
        index = split_line[0][1:]
        name = split_line[1]
        shift = split_line[2]
        if((name != 'N') and (name != 'NE2')):
            cs_lst.append([index, name, shift])

#generate cross_peaks list
def generate_cross_peaks_lst():
    for j in range(len(list_files)):
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

def analyze_ambiguity():
    cs_shifts = np.array([obj[2] for obj in cs_lst])
    cs_shifts = cs_shifts.astype(np.float)
    for i in range(len(cross_peaks)):
        f2_val = float(cross_peaks[i][0])
        f1_val = float(cross_peaks[i][1])
        f2_ambig_cs_lst = np.where((cs_shifts >= f2_val - error) & (cs_shifts <= f2_val + error))
        f1_ambig_cs_lst = np.where((cs_shifts >= f1_val - error) & (cs_shifts <= f1_val + error))
        possible_assignments_for_peak = []
        cs_arr = np.array(cs_lst)
        for k in range(len(cs_shifts[f2_ambig_cs_lst])):
            for l in range(len(cs_shifts[f1_ambig_cs_lst])): 
                if (curr_scheme == 'fully'):
                    possible_assignments_for_peak.append([cs_arr[f2_ambig_cs_lst][k].tolist(),cs_arr[f1_ambig_cs_lst][l].tolist()])
                else:
                    atom_f2_lst = cs_arr[f2_ambig_cs_lst][k].tolist()
                    atom_f1_lst = cs_arr[f1_ambig_cs_lst][l].tolist()

                    atom_f2_name =  atom_f2_lst[1]
                    atom_f1_name =  atom_f1_lst[1]
                    if(curr_scheme == '1,3-glycerol'):
                        atom_f2_prob = float(onethreeglycerol_dict[sequence[int(atom_f2_lst[0])-1]][atom_f2_name])
                        atom_f1_prob = float(onethreeglycerol_dict[sequence[int(atom_f1_lst[0])-1]][atom_f1_name])
                    elif(curr_scheme == '2-glycerol'):
                        atom_f2_prob = 1 - float(onethreeglycerol_dict[sequence[int(atom_f2_lst[0])-1]][atom_f2_name])
                        atom_f1_prob = 1 - float(onethreeglycerol_dict[sequence[int(atom_f1_lst[0])-1]][atom_f1_name])
                    multiplied = atom_f1_prob*atom_f2_prob
                    if((multiplied < threshold) or (multiplied == 0)):
                        continue
                    else:
                        possible_assignments_for_peak.append([cs_arr[f2_ambig_cs_lst][k].tolist(),cs_arr[f1_ambig_cs_lst][l].tolist()])
        dist = 0
        ambig = len(possible_assignments_for_peak)
        if(ambig == 1):
            dist = update_distance_counters(possible_assignments_for_peak[0])
        
        if(ambig <= int(max_ambig)):
            test_lst = [cross_peaks[i], possible_assignments_for_peak, ambig, dist]
            possible_assign_per_peak_lst.append([cross_peaks[i], possible_assignments_for_peak, ambig, dist, peak_source_lst[i]])
            possible_correlations = []
            for l in range(len(possible_assignments_for_peak)):
                atom_a_name = str(possible_assignments_for_peak[l][0][0]) + str(possible_assignments_for_peak[l][0][1])
                atom_b_name = str(possible_assignments_for_peak[l][1][0]) + str(possible_assignments_for_peak[l][1][1])
                current_correlation = [atom_a_name, atom_b_name] 
                possible_correlations.append(current_correlation)
            cross_peaks_with_neighbors.append([cross_peaks[i], possible_correlations])

def generate_xplor_restraints(xplor_data):
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

def write_analysis_file():
    with open(analysis_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(possible_assign_per_peak_lst)

def aggregate_by_dup_perm_rev(seq, ambig_dict, contact_type_dict):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        sorted_lst = sorted(item[1])
        sorted_lst_rev = [elem[::-1] for elem in sorted_lst]
        rev_and_not_rev = [sorted_lst, sorted_lst_rev]
        lst_for_key = sorted(rev_and_not_rev)[0]
        tally[json.dumps(lst_for_key)].append(i)
        ambig_dict[json.dumps(lst_for_key)].append(item[2])
        contact_type_dict[json.dumps(lst_for_key)].append(item[3])
    return ((key,locs) for key,locs in tally.items())

def remove_dup_perm_rev(source):
    ambig_dict = defaultdict(list)
    contact_type_dict = defaultdict(list)
    final_dup_data = []
    for dup in sorted(aggregate_by_dup_perm_rev(source, ambig_dict, contact_type_dict)):
        assignments_lst = json.loads(dup[0])
        cross_peaks_lst = [source[i][0] for i in range(len(source)) if i in dup[1]]
        cross_peaks_sources_lst = [source[i][4] for i in range(len(source)) if i in dup[1]]
        ambig = ambig_dict[dup[0]][0]
        contact_type = contact_type_dict[dup[0]][0]
        entry = [cross_peaks_lst, assignments_lst, ambig, contact_type, cross_peaks_sources_lst]
        final_dup_data.append(entry)
    return final_dup_data

def generate_aggregated_data_for_xplor(aggregated_data):
    aggregated_data_for_xplor = []

    for i in range(len(aggregated_data)):
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

def write_params_file():
    params_file_name = file_name + "_params.txt"
    file_handler = open(params_file_name, "w")
    file_handler.write("run name: " + file_name + "\n")
    file_handler.write("labeling scheme: " + curr_scheme + "\n")
    file_handler.write("low restraint: " + str(low_restraint) + "\n")
    file_handler.write("high restraint: " + str(high_restraint) + "\n")
    file_handler.write("error: " + str(error) + "\n")
    file_handler.write("threshold: " + str(threshold) + "\n")
    file_handler.write("list of files: " + str(list_file_names) + "\n")
    file_handler.write("max ambig.: " + str(max_ambig) + "\n")
    file_handler.close()

def generate_csv_of_xplor_restraints(agg_data):
    for i in range(len(agg_data)):
        ambig_level = agg_data[i][2]
        if(ambig_level == 1):
            first_ind = agg_data[i][1][0][0][0]
            first_name = agg_data[i][1][0][0][1]
            second_ind = agg_data[i][1][0][1][0]
            second_name = agg_data[i][1][0][1][1]
            if((first_ind == second_ind) and (first_name == second_name)):
                continue     
            #if unambiguous single-bond restraints  delete restraint
            current_names = [first_name.capitalize(), second_name.capitalize()]
            if((abs(int(first_ind) - int(second_ind)) == 0)) and (is_single_bond(current_names)):
                continue
        #AMBIGUOUS CASES
        else:
            single_bond_count = 0
            includes_atom_atom_count = 0
            for k in range(ambig_level):
                index_a = agg_data[i][1][k][0][0]
                name_a = agg_data[i][1][k][0][1]
                index_b = agg_data[i][1][k][1][0]
                name_b = agg_data[i][1][k][1][1]
                #if at least one of the possibilities is X-X - delete the entire restraint
                if((index_a == index_b) and (name_a == name_b)):
                    includes_atom_atom_count = includes_atom_atom_count + 1
                current_names = [name_a.capitalize(), name_b.capitalize()]
                #count the number of possibilites that are single-bond. if it is larger than 1 - delete the entire restraint.
                if((abs(int(index_a) - int(index_b)) == 0)) and (is_single_bond(current_names)):
                    single_bond_count = single_bond_count + 1
            if((includes_atom_atom_count != 0)):
                continue     
        csv_xplor_data.append(aggregated_data[i])

generate_cs_lst()
generate_cross_peaks_lst()
analyze_ambiguity()
write_analysis_file()

aggregated_data = remove_dup_perm_rev(possible_assign_per_peak_lst)
aggregated_input_for_xplor_restraints = generate_aggregated_data_for_xplor(aggregated_data)
generate_xplor_restraints(aggregated_input_for_xplor_restraints)
generate_csv_of_xplor_restraints(aggregated_data)

for k in range(len(aggregated_data)):
    if(aggregated_data[k][2] == 0):
        aggregated_data.pop(k)

with open(analysis_aggregated_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(aggregated_data)

for i in range(len(csv_xplor_data)):
    if(csv_xplor_data[i][2] == 0):
        csv_xplor_data.pop(i)

csv_xplor_file_name = file_name + "_xplor.csv"

with open(csv_xplor_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(csv_xplor_data)

print(f"Length aggregated data: {len(aggregated_data)}, length xplor data: {len(csv_xplor_data)}")
for i in range(len(list_files)):
    list_files[i].close()
shifts_file.close()
sequence_file.close()
xplor_file.close()
sequence_file.close()

write_params_file()