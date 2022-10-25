
import json
import csv
import sys
import numpy as np
import copy
from collections import defaultdict
import re

file_name = "runE_all"
fully_file_name = "analysis_aggregated_runE_fully.csv"
sparse_file_name = "analysis_aggregated_runE_sparse.csv"
low_restraint = 2
high_restraint = 7.5
delta = 0.3

records_file = file_name + "_with_records.csv"
xplor_file_name = file_name + "_xplor.tbl"
xplor_file = open(xplor_file_name, "w")

single_bond_contacts = [["Ca", "Cb"],["C", "Ca"], ["Cb", "Cg"],["Cg", "Cd"],["Cg", "Cd2"],["Cb","Cg1"],["Cb","Cg2"],["Cg1","Cd1"],["Cd","Ce"],["Cg","Cd1"],["Cd1","Ce1"],["Cd2","Ce2"],["Ce1","Cz"],["Ce2","Cz"], ["Cd2","Ce3"],["Ce3","Cz3"], ["Cz3", "Ch2"], ["Cz2", "Ch2"], ["Ce2", "Cz2"]]

cross_peaks = []
sparse_data = []
fully_data = []
csv_xplor_data = []
#

def populate_data_list(file_name, empty_lst):
    with open(file_name, newline = '') as f:
        reader = csv. reader(f)
        data = list(reader)
        for i in range(len(data)):
            curr_entry_cross_peaks = data[i][0]
            curr_entry_assignments = data[i][1]
            curr_entry_ambig = data[i][2]
            curr_entry_contact_type = data[i][3]
            curr_entry_cross_peaks_replaced = curr_entry_cross_peaks.replace("'", '"')
            curr_entry_assignments_replaced = curr_entry_assignments.replace("'", '"')
            curr_entry_ambig_replaced = curr_entry_ambig.replace("'", '"')
            curr_entry_contact_type_replaced = curr_entry_contact_type.replace("'", '"')
            lst_currr_entry_cross_peaks = json.loads(curr_entry_cross_peaks_replaced)
            lst_currr_entry_assignments = json.loads(curr_entry_assignments_replaced)
            lst_currr_entry_ambig = json.loads(curr_entry_ambig_replaced)
            lst_currr_entry_contact_type = json.loads(curr_entry_contact_type_replaced)
            curr_entry_sources = data[i][4]
            curr_entry_sources_replaced = curr_entry_sources.replace("'", '"')
            curr_entry_sources_lst = json.loads(curr_entry_sources_replaced)
            entry = [lst_currr_entry_cross_peaks, lst_currr_entry_assignments, lst_currr_entry_ambig, lst_currr_entry_contact_type, curr_entry_sources_lst]
            empty_lst.append(entry)

#create a list with one representative peak for each entry in the data_list from csv
def representative_peaks(data):
    rep_lst = []
    #iterate over all cross_peaks in the .csv file
    for i in range(len(data)):
        rep_cross_peak = data[i][0][0]
        rep_lst.append(rep_cross_peak)
    return rep_lst

def get_indexes_of_fully_to_override(sparse, fully, fully_copy, delta, fully_entire):

    fully_first_shifts = np.array([float(elem[0]) for elem in fully])
    fully_second_shifts = np.array([float(elem[1]) for elem in fully])
    ret_list = []

    for i in range(len(sparse)):
        first_shift = float(sparse[i][0])
        second_shift = float(sparse[i][1])

        first_second = np.where((fully_first_shifts >= first_shift - delta) & (fully_first_shifts <= first_shift + delta) & (fully_second_shifts >= second_shift - delta) & (fully_second_shifts <= second_shift + delta))
        second_first = np.where((fully_first_shifts >= second_shift - delta) & (fully_first_shifts <= second_shift + delta) & (fully_second_shifts >= first_shift - delta) & (fully_second_shifts <= first_shift + delta))

        joined_lists = list(set(list(first_second[0]) + list(second_first[0])))
        #keep record - which cross-peak from the sparse file caused this fully cross-peak to be overriden
        for s in range(len(joined_lists)):
            fully_copy[joined_lists[s]].append([first_shift, second_shift])
            fully_entire[joined_lists[s]][0] = ""
        ret_list.append([i, joined_lists])
    return ret_list

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

def append_non_overriden_fully_to_sparse(fully ,sparse):
    for i in range(len(fully)):
        if(fully[i][0] != ""):
            sparse.append(fully[i])

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
    sources_dict = defaultdict(list)
    final_dup_data = []
    for dup in sorted(aggregate_by_dup_perm_rev(source, ambig_dict, contact_type_dict)):
        assignments_lst = json.loads(dup[0])
        cross_peaks_lst = [source[i][0] for i in range(len(source)) if i in dup[1]]
        cross_peaks_sources_lst = [source[i][4] for i in range(len(source)) if i in dup[1]]
        flat_cross_peak_list = [item for sublist in cross_peaks_lst for item in sublist]
        flat_sources_list = [item for sublist in cross_peaks_sources_lst for item in sublist]
        ambig = ambig_dict[dup[0]][0]
        contact_type = contact_type_dict[dup[0]][0]
        entry = [flat_cross_peak_list, assignments_lst, ambig, contact_type, flat_sources_list]
        if(len(flat_cross_peak_list) != len(flat_sources_list)):
            print('yes')
        final_dup_data.append(entry)
    return final_dup_data

def is_single_bond(current):
    global single_bond_contacts
    reversed_current = current[::-1]
    if current in single_bond_contacts or reversed_current in single_bond_contacts:  
        return True
    else:
        return False

#XPLOR RESTRAINTS
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


def generate_csv_of_xplor_restraints(agg_data):
    for i in range(len(agg_data)):
        ambig_level = agg_data[i][2]
        #UNAMBIGUOUS CASES
        if(ambig_level == 1):
            first_ind = agg_data[i][1][0][0][0]
            first_name = agg_data[i][1][0][0][1]
            second_ind = agg_data[i][1][0][1][0]
            second_name = agg_data[i][1][0][1][1]
            if((first_ind == second_ind) and (first_name == second_name)):
                continue     
            #if unambiguous single-bond restraints, delete restraint
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
        csv_xplor_data.append(sparse_data_reduced[i])

populate_data_list(sparse_file_name, sparse_data)
populate_data_list(fully_file_name, fully_data)

fully_data_copy_with_records = copy.deepcopy(fully_data)

sparse_rep = representative_peaks(sparse_data)
fully_rep = representative_peaks(fully_data)

joined_lists = get_indexes_of_fully_to_override(sparse_rep, fully_rep, fully_data_copy_with_records, delta, fully_data)

append_non_overriden_fully_to_sparse(fully_data, sparse_data)

print("length sparse: ", len(sparse_data))
sparse_data_reduced = remove_dup_perm_rev(sparse_data)
print("length sparse reduced: ", len(sparse_data_reduced))

sparse_data_reduced_for_xplor = generate_aggregated_data_for_xplor(sparse_data_reduced)

generate_xplor_restraints(sparse_data_reduced_for_xplor)
generate_csv_of_xplor_restraints(sparse_data_reduced)

with open(records_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(fully_data_copy_with_records)

override_file_name = file_name + "_override.csv"
final_file_name = file_name + "_final.csv"

with open(override_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(fully_data)

with open(final_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(sparse_data_reduced)

csv_xplor_file_name = file_name + "_xplor.csv"

with open(csv_xplor_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(csv_xplor_data)

def write_params_file():
    params_file_name = file_name + "_params.txt"
    file_handler = open(params_file_name, "w")
    file_handler.write("run name: " + file_name + "\n")
    file_handler.write("fully file: " + fully_file_name + "\n")
    file_handler.write("sparse file: " + sparse_file_name + "\n")
    file_handler.write("low restraint: " + str(low_restraint) + "\n")
    file_handler.write("high restraint: " + str(high_restraint) + "\n")
    file_handler.write("delta: " + str(delta) + "\n")
    file_handler.close()

write_params_file()