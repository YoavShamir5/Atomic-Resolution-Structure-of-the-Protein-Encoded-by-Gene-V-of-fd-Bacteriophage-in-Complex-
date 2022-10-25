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

from numpy.lib.arraysetops import _in1d_dispatcher

input_csv_xplor_file_name = "runE_filtered_by_16A_xplor_no_sb"
updated_csv_file_name = input_csv_xplor_file_name + "_no_dimer_contacts.csv"
ruled_out_csv_file_name = input_csv_xplor_file_name + "_ruled_out_restraints.csv"
updated_tbl_file_name = "runMaxAmbigFive_initial_no_DR.tbl"
low_restraint = 2
high_restraint = 8
dimer_dist = 7

single_bond_contacts = [["Ca", "Cb"],["C", "Ca"], ["Cb", "Cg"],["Cg", "Cd"],["Cg", "Cd2"],["Cb","Cg1"],["Cb","Cg2"],["Cg1","Cd1"],["Cd","Ce"],["Cg","Cd1"],["Cd1","Ce1"],["Cd2","Ce2"],["Ce1","Cz"],["Ce2","Cz"], ["Cd2","Ce3"],["Ce3","Cz3"], ["Cz3", "Ch2"], ["Cz2", "Ch2"], ["Ce2", "Cz2"]]
imported_csv_data = []
updated_csv_data_no_sb = []
ruled_out_contacts = []
structure_file = "1vqb.pdb1"
parser = Bio.PDB.PDBParser(QUIET=True)
structures = parser.get_structure('protein', structure_file)

monomer_a_structure = structures[0]
monomer_b_structure = structures[1]
chain_a = monomer_a_structure["A"]
chain_b = monomer_b_structure["A"]

def import_data_from_csv_and_update_restraints():
    global possibilities_test_counter
    file_name = input_csv_xplor_file_name + ".csv"
    with open(file_name, newline = '') as f:
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
            imported_csv_data.append(curr_entry)

def is_single_bond(current):
    global single_bond_contacts
    reversed_current = current[::-1]
    if current in single_bond_contacts or reversed_current in single_bond_contacts:  
        return True
    else:
        return False

def find_inter_monomer_dist(corr):

    index_a = int(corr[0][0])
    name_a = corr[0][1]
    index_b = int(corr[1][0])
    name_b = corr[1][1]

    dist = chain_a[index_a][name_a] - chain_b[index_b][name_b]

    if(dist < dimer_dist):
        print(f"{index_a}/{name_a} - {index_b}/{name_b}, {dist}")
    return dist


#iterate over all restraints - if a restraint includes even one single-bond possibility, rule out the entire restraint
def filter_csv_data_by_single_bond_correlations():
    for k in range(len(imported_csv_data)):        
        try:
            if(any((find_inter_monomer_dist(i) < dimer_dist) for i in imported_csv_data[k][1])):
                ruled_out_contacts.append(imported_csv_data[k])
                continue
        except KeyError as e:
            updated_csv_data_no_sb.append(imported_csv_data[k])
        else:
            updated_csv_data_no_sb.append(imported_csv_data[k])

def generate_updated_data_for_xplor(aggregated_data):
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

def generate_xplor_restraints(xplor_data):
    xplor_file = open(updated_tbl_file_name, "w")
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

import_data_from_csv_and_update_restraints()
filter_csv_data_by_single_bond_correlations()
updated_data_for_tbl = generate_updated_data_for_xplor(updated_csv_data_no_sb)
generate_xplor_restraints(updated_data_for_tbl)

with open(updated_csv_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(updated_csv_data_no_sb)


with open(ruled_out_csv_file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(ruled_out_contacts)

def write_params_file():
    params_file_name = input_csv_xplor_file_name + "_params.txt"
    file_handler = open(params_file_name, "w")

    file_handler.write("dimer dist: " + str(dimer_dist) + "\n")
    file_handler.write(f"input file: {input_csv_xplor_file_name}")
    file_handler.close()

write_params_file()