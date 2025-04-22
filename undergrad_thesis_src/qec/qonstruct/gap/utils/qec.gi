#
# Author:   Suhas Vittal
# date:     19 January 2024
#
# Utility code for QEC.
#

LoadPackage("IO");

f_OperatorWeight :=
    function(operator)
        local w, x;
        w := 0;
        for x in operator do
            if x = One(GF(2)) then
                w := w+1;
            fi;
        od;
        return w;
    end;

f_ComputeStabilizerGenerators :=
    function(check_list)
        local M;
        M := SemiEchelonMat(check_list);
        return check_list{ Filtered(M.heads, x -> x > 0) };
    end;

f_AreCommuting :=
    function(v, w)
        return (v * w) = Zero(GF(2));
    end;

# Computes the logical X/Z operators. To compute logical X operators,
# the first argument should be the X checks and the second should be
# the Z checks. For logical Z, vice versa.
f_CssLXZOperators := 
    function(same_type_checks, opposing_checks)
        local im_of_same, im_of_opp, ker, op_list;

        im_of_same := VectorSpace(GF(2), same_type_checks); 
        
        op_list := [];
        Append(op_list, same_type_checks);
        Append(op_list, NullspaceMat(TransposedMat(opposing_checks)));

        ker := VectorSpace(GF(2), op_list);
        return Filtered(BasisVectors(Basis(ker)), v -> not (v in im_of_same));
    end;

f_MinOperatorWeight :=
    function(operators) 
        local min_w, op, w;

        min_w := f_OperatorWeight(operators[1]);
        for op in operators do
            w := f_OperatorWeight(op);
            if w < min_w then   min_w := w; fi;
        od;
        return min_w;
    end;

# Prints generators for tanner graph.
f_PrintGenerators :=
    function(stabilizers, type, tabs)
        local i, j, s, x;

        i := 0;
        for s in stabilizers do
            Print(tabs, type, i);
            j := 0;
            for x in s do
                if x = One(GF(2)) then
                    Print(",", j);
                fi;
                j := j+1;
            od;
            Print("\n");
            i := i+1;
        od;
    end;

f_WriteChecksAndOperatorsToFile :=
    function(fout, checks, operators, type)
        local line, i, j, s, x;

        i := 0;
        for s in checks do
            line := Concatenation(type, String(i));
            j := 0;
            for x in s do
                if x = One(GF(2)) then
                    line := Concatenation(line, ",", String(j));
                fi;
                j := j+1;
            od;
            i := i+1;

            IO_WriteLine(fout, line);
        od;
        i := 0;
        for s in operators do
            line := Concatenation("O", type, String(i));
            j := 0;
            for x in s do
                if x = One(GF(2)) then
                    line := Concatenation(line, ",", String(j));
                fi;
                j := j+1;
            od;
            i := i+1;

            IO_WriteLine(fout, line);
        od;
    end;

f_WriteCodeToFile :=
    function(file, x_checks, z_checks, x_operators, z_operators)
        local fout;

        fout := IO_File(file, "w");
        Print("Creating file ", file, ": ", IsFile(fout), "\n");
        f_WriteChecksAndOperatorsToFile(fout, x_checks, x_operators, "X");
        f_WriteChecksAndOperatorsToFile(fout, z_checks, z_operators, "Z");
    end;
