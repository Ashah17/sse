#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Precondition: run test_subgroup.gi.
#
# Load: verbose (true if print out debugging messages).

f_ComputeIndicatorVectors :=
    function(check_cosets, qubit_cosets)
        local G, H, common, v_check_list, v_check;

        v_check_list := [];
        for G in check_cosets do
            v_check := [];
            for H in qubit_cosets do
                common := Intersection(AsList(H), AsList(G));
                if Length(common) > 0 then
                    Add(v_check, One(GF(2)));
                else
                    Add(v_check, Zero(GF(2)));
                fi;
            od;
            Add(v_check_list, v_check);
        od;
        return v_check_list;
    end;

f_CodeFilename :=
    function(n, k)
        local filename;

        filename := Concatenation(output_folder, "/", 
                                    String(n), "_", String(k), 
                                    "_g", String(use_group),
                                    ".txt");
        return filename;
    end;

done := false;

ii := 1;
while ii <= Length(iso_bc) do
    jj := 1;
    while jj <= Length(iso_ca) do
        kk := 1;
        while kk <= Length(iso_ab) do
            H_bc := Image(iso_bc[ii]);
            H_ca := Image(iso_ca[jj]);
            H_ab := Image(iso_ab[kk]);

            faces := RightCosets(G_rs_mod_H, H_bc);
            vertices := RightCosets(G_rs_mod_H, H_ca);
            edges := RightCosets(G_rs_mod_H, H_ab);

            genus := 2 - (Length(vertices) - Length(edges) + Length(faces));

            n_data := Length(edges);
            # Compute checks. A data qubit is in the support of an X (Z) check if
            # its edge coset shares an element with the check's face (vertex) coset.
            x_checks := f_ComputeIndicatorVectors(vertices, edges);
            z_checks := f_ComputeIndicatorVectors(faces, edges);
            # Check that the check weights are correct.
            if not All(List(x_checks, v -> f_OperatorWeight(v) = s)) 
            or not All(List(z_checks, v -> f_OperatorWeight(v) = r)) 
            then 
                kk := kk+1;
                continue;
            fi;
            # Compute stabilizer generators:
            x_checks := f_ComputeStabilizerGenerators(x_checks);
            z_checks := f_ComputeStabilizerGenerators(z_checks);
            # Compute logical operators:
            x_operators := f_CssLXZOperators(x_checks, z_checks);
            z_operators := f_CssLXZOperators(z_checks, x_checks);

            n_x_checks := Length(x_checks);
            n_z_checks := Length(z_checks);
            n_checks := n_x_checks + n_z_checks;
            n_x_ops := Length(x_operators);
            n_z_ops := Length(z_operators);

            genus := n_data - n_checks;
            if n_data - n_checks = genus and n_x_ops = genus and n_z_ops = genus then
                Print("\tData qubits: ", n_data, "\n");
                Print("\tLogical qubits: ", genus, "\n");
                Print("\tChecks: ", n_checks, "\n");
                Print("\t\tX: ", n_x_checks, ", Z: ", n_z_checks, "\n");
                Print("\tX/Z Operators: ", Length(x_operators), ", ", Length(z_operators), "\n");
                Print("\t\tX min weight: ", f_MinOperatorWeight(x_operators), "\n");
                Print("\t\tZ min weight: ", f_MinOperatorWeight(z_operators), "\n");

                f_WriteCodeToFile(f_CodeFilename(n_data, genus),    
                                    x_checks,
                                    z_checks,
                                    x_operators,
                                    z_operators);
            fi;
            kk := kk+1;
        od;
        if done then break; fi;
        jj := jj+1;
    od;
    if done then break; fi;
    ii := ii+1;
od;
