#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Load:
#   r and s -- the Schlaffi symbols for the surface.
#   max_index -- max index of normal subgroup
#

# Get functions.
Read("gap/utils/groups.gi");
Read("gap/utils/qec.gi");
Read("gap/utils/operations.gi");

LoadPackage("LINS");

# Compute normal subgroups.

G_base := FreeGroup(3);
G_base := G_base / [ G_base.1^2, G_base.2^2, G_base.3^2, 
                    (G_base.1*G_base.2)^2, 
                    (G_base.2*G_base.3)^r,
                    (G_base.3*G_base.1)^s ];
G_rs := Subgroup(G_base, [G_base.2*G_base.3, G_base.3*G_base.1]);

gr := LowIndexNormalSubgroupsSearch(G_rs, max_index);
normal_low_index_subgroups := ComputedNormalSubgroups(gr);

Print("# of normal subgroups = ", Length(normal_low_index_subgroups), "\n");
Print("Indices:");
for H in normal_low_index_subgroups do
    Print(" ", Index(H));
od;
Print("\n");

# A quirk of LINS: the root of the graph is different from G_rs,
# so we need to update G_rs.

gr_root := LinsRoot(gr);
_G_rs := G_rs;
G_rs := Grp(gr_root);

use_group := 1;
