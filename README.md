Atomic-Resolution Structure of the Protein Encoded by Gene V of fd Bacteriophage in Complex with Viral ssDNA Determined by Magic-Angle Spinning Solid-State NMR
Python Scripts used in the paper Atomic-Resolution Structure of the Protein Encoded by Gene V of fd Bacteriophage in Complex with Viral ssDNA Determined by Magic-Angle Spinning Solid-State NMR



1. remove_spinning_sidebands.py - remove peaks attributed to the first spinning sideband from Sparky peak lists provided as input.
2. generate_restraints_from_peak_lists.py - generate Ambiguous Distance Restraints based on input peak lists an user-provided parameters.
3. aggregate_fully_and_sparse_restraints.py - provided with two sets of generated restraints, one from fully labeled data and the other 
   from sparsely labeled data, generates an aggreagted set of distance restraints.
4. filter_restraints_by_cutoff_distance_in_free_protein.py - uses the input of distance restraints and rules out possible assignments
   corresponding to distance restraints larger than a user-defined cutoff in the free protein structure.
5. rule_out_potential_single_bond_restraints.py - provided with an input of distance restraints, rules out any restraint which includes
   one or more single-bond contacts in its list of possible assignments. 
6. rule_out_potential_dimer_restraints.py - given an input of distance restraints, measures distance in the free protein homodimer, and
   rules out restraints which include one or more possible assignments corresponding to a short (by a cutoff) inter-monomer distance in the
   free protein homodimer.
7. analyze_dihedral_violations_after_Xplor_run.py - given a set of torsion angle restraints that was provided to an Xplor run as well as the Xplor
   statistics output, rules out restraints that are violated beyond a user-defined cutoff value.
8. analyze_distance_violations_after_Xplor_run.py - given a set of distance restraints that was provided to an Xplor run as well as the Xplor
   statistics output, rules out restraints that are violated beyond a user-defined cutoff value.
9. analyze_by_cutoff_from_average_structure_after_Xplor_run.py - given a set of distance restraints that was provided to an Xplor run as well as
the average of the k lowest energy structures of the Xplor run, rules out possible assignments further apart than some cutoff distance in the 
average structure. 
