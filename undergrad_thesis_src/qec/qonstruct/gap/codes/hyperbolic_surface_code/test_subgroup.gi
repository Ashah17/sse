#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Precondition: run init_surface.gi
#
# Load:
#   use_group: index of subgroup.
#   verbose: print debug
#

H := normal_low_index_subgroups[use_group];
tfn_subgroup := Grp(H);

G_rs_mod_H := G_rs / tfn_subgroup;

iso_bc := IsomorphicSubgroups(G_rs_mod_H, CyclicGroup(r));
iso_ca := IsomorphicSubgroups(G_rs_mod_H, CyclicGroup(s));
iso_ab := IsomorphicSubgroups(G_rs_mod_H, CyclicGroup(2));

if verbose then
    Print("\t# of isomorphisms for BC, CA, AB: ",
                                Length(iso_bc), ", ",
                                Length(iso_ca), ", ",
                                Length(iso_ab), "\n");
fi;
