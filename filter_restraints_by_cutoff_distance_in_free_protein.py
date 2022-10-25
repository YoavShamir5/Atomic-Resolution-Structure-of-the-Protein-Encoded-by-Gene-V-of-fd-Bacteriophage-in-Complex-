import csv
import sys
from Bio.Data.IUPACData import protein_letters_3to1
import json 
import Bio
from Bio.PDB import PDBParser
import Bio.PDB
from Bio.PDB.Structure import Structure
from Bio.Data.IUPACData import protein_letters_1to3
import copy
import re  
from collections import defaultdict

longer_test_counter = 0
possibilities_test_counter = 0

max_dist = 16
all_possibilites = []
possibilites_csv_file_name = "runE_all_final.csv"
pdb_file_name = "1VQB.pdb1"
data_name = "runE"
tmp_updated_data_file_name = data_name + "_filtered_by_" + str(max_dist) + "A"
updated_data_file_name = tmp_updated_data_file_name + ".csv"
updated_agregated_data_file_name = "aggregated_" + data_name + "_filtered_by_" + str(max_dist) + "A.csv"
xplor_file_name = "runMaxAmbigFive_initial.tbl"

updated_data = []
csv_xplor_data = []
single_bond_contacts = [["Ca", "Cb"],["C", "Ca"], ["Cb", "Cg"],["Cg", "Cd"],["Cg", "Cd2"],["Cb","Cg1"],["Cb","Cg2"],["Cg1","Cd1"],["Cd","Ce"],["Cg","Cd1"],["Cd1","Ce1"],["Cd2","Ce2"],["Ce1","Cz"],["Ce2","Cz"], ["Cd2","Ce3"],["Ce3","Cz3"], ["Cz3", "Ch2"], ["Cz2", "Ch2"], ["Ce2", "Cz2"]]
low_restraint = 2
high_restraint = 7.5

parser = Bio.PDB.PDBParser(QUIET=True)
structures = parser.get_structure('protein', pdb_file_name)
structure_a = structures[0]
chain_a = structure_a['A']

def find_dist(index_a, name_a, index_b, name_b):
    global longer_test_counter
    arginine_21_unmodelled_atoms = ['21CG', '21CD', '21CZ']
    #some carbons are not modelled in 1VQB.PDB1 - handeling those specifically - 21CG, 21CD, 21CZ
    if(str(index_a) == '21' and name_a == 'CG' and str(index_b) != '21'):
        dist = chain_a[21]['CB'] - chain_a[index_b][name_b] + 1.5
    elif((str(index_a) == '21' and name_a == 'CD') and str(index_b) != '21'):
        dist = chain_a[21]['CB'] - chain_a[index_b][name_b] + 3
    elif((str(index_a) == '21' and name_a == 'CZ')  and str(index_b) != '21'):
        dist = chain_a[21]['CB'] - chain_a[index_b][name_b] + 6
    ##
    elif(str(index_b) == '21' and name_b == 'CG' and str(index_a) != '21'):
        dist = chain_a[index_a][name_a] - chain_a[21]['CB'] + 1.5
    elif((str(index_b) == '21' and name_b == 'CD') and str(index_a) != '21'):
        dist = chain_a[index_a][name_a] - chain_a[21]['CB'] + 3
    elif((str(index_b) == '21' and name_b == 'CZ')  and str(index_a) != '21'):
        dist = chain_a[index_a][name_a] - chain_a[21]['CB'] + 6
    elif(str(index_a) == '21' and str(index_b) == '21'):
        dist = 5
    else:
        dist = chain_a[index_a][name_a] - chain_a[index_b][name_b]
    if(dist > max_dist):
        longer_test_counter = longer_test_counter + 1
    return dist

def calc_type_of_unambig(index_first, index_second):
    #if intra
    if(int(index_first) == int (index_second)):
        return 1  
    #if seq  
    elif(abs(int(index_first) - int(index_second)) == 1):
        return 2
    #if medium
    elif(abs(int(index_first) - int(index_second)) <= 4):
        return 3
    #else - long
    else:
        return 4

def write_updated_data_to_csv():
    with open(updated_data_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(updated_data)

def import_all_possibilites_from_csv():
    global possibilities_test_counter
    with open(possibilites_csv_file_name, newline = '') as f:
        reader = csv. reader(f)
        data = list(reader)
        for i in range(len(data)):
            curr_entry_cross_peaks = data[i][0]
            curr_entry_assignments = data[i][1]
            curr_entry_ambig = data[i][2]
            curr_entry_contact_type = data[i][3]
            curr_entry_sources_lst = data[i][4]
            curr_entry_cross_peaks_replaced = curr_entry_cross_peaks.replace("'", '"')
            curr_entry_assignments_replaced = curr_entry_assignments.replace("'", '"')
            curr_entry_ambig_replaced = curr_entry_ambig.replace("'", '"')
            curr_entry_contact_type_replaced = curr_entry_contact_type.replace("'", '"')
            curr_entry_sources_replaced = curr_entry_sources_lst.replace("'", '"')
            lst_currr_entry_cross_peaks = json.loads(curr_entry_cross_peaks_replaced)
            lst_currr_entry_assignments = json.loads(curr_entry_assignments_replaced)
            lst_currr_entry_ambig = json.loads(curr_entry_ambig_replaced)
            lst_curr_entry_sources = json.loads(curr_entry_sources_replaced)
            lst_currr_entry_contact_type = json.loads(curr_entry_contact_type_replaced)
            curr_entry = [lst_currr_entry_cross_peaks, lst_currr_entry_assignments, lst_currr_entry_ambig, lst_currr_entry_contact_type, lst_curr_entry_sources]

            updated_possibilites = []
            for j in range(len(curr_entry[1])):
                possibilities_test_counter = possibilities_test_counter + 1
                first_index = int(curr_entry[1][j][0][0])
                first_name = curr_entry[1][j][0][1]
                second_index = int(curr_entry[1][j][1][0])
                second_name = curr_entry[1][j][1][1]
                dist = find_dist(first_index, first_name, second_index, second_name)
                if(dist <= max_dist):
                    updated_possibilites.append(curr_entry[1][j])
            updated_ambig = len(updated_possibilites)
            if(updated_ambig == 1):
                updated_ambig_type = calc_type_of_unambig(int(updated_possibilites[0][0][0]), int(updated_possibilites[0][1][0]))
            else:
                updated_ambig_type = 0

            updated_entry = [lst_currr_entry_cross_peaks, updated_possibilites, updated_ambig, updated_ambig_type, lst_curr_entry_sources]
            if(len(lst_curr_entry_sources) != len(lst_currr_entry_cross_peaks)):
                print("YES")
            updated_data.append(updated_entry)

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

import_all_possibilites_from_csv()
write_updated_data_to_csv()
print("length before: ", len(updated_data))
aggregated_data = remove_dup_perm_rev(updated_data)
print("length after: ", len(aggregated_data))
aggregated_data_for_xplor = generate_aggregated_data_for_xplor(aggregated_data)
generate_xplor_restraints(aggregated_data_for_xplor)
generate_csv_of_xplor_restraints(aggregated_data)

for i in range(len(csv_xplor_data)):
    if(csv_xplor_data[i][2] == 0):
        csv_xplor_data.pop(i)

csv_xplor_file_name = tmp_updated_data_file_name + "_xplor.csv"

with open(csv_xplor_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(csv_xplor_data)

for i in range(len(aggregated_data)):
    if(aggregated_data[i][2] == 0):
        aggregated_data.pop(i)

with open(updated_agregated_data_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(aggregated_data)
