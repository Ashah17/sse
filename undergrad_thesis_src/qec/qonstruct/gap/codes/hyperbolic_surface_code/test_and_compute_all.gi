#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Precondition: run init_surface.gi
#

verbose := false;

# Create the output folder.
LoadPackage("IO");
output_folder := Concatenation("codes/hy_sc/", String(r), "_", String(s));
Print("Make output folder ", output_folder , ": ", IO_mkdir(output_folder, 448), "\n");

while use_group <= Length(normal_low_index_subgroups) do
    Print("Subgroup ", use_group, ", Index = ", Index(normal_low_index_subgroups[use_group]), "\n");

    Read("gap/codes/hyperbolic_surface_code/test_subgroup.gi");
    Read("gap/codes/hyperbolic_surface_code/compute_code.gi");

    use_group := use_group+1;
od;
