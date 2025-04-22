#
# Author:   Suhas Vittal
# date:     19 January 2024
#
# Utility code for Groups.
#

LoadPackage("Gauss");

f_ListOrders :=
    function(G)
        local orders, x;
        orders := List(Elements(G), x -> Order(x));

        Print("Orders:");
        for x in orders do
            Print(" ", x);
        od;
        Print("\n");
    end;

# Returns true if H is suspected to be Torsion-Free.
f_IsTorsionFree := 
    function(H)
        local g, ord, result;
        for g in GeneratorsOfGroup(H) do
            result := IO_CallWithTimeoutList(rec(seconds := 20), Order, [g]);
            if result[1] then
                ord := result[2];
                if not (ord = infinity) then
                    return false;
                fi;
            fi;
        od;
        return true;
    end;

# Returns a record of X and Z operators (.x and .z respectively).
f_ComputeOperators :=
    function(x_check_list, z_check_list)
        local im_x,
                im_z,
                ker_x,
                ker_z,
                z_op_list,
                x_op_list;

        im_x := VectorSpace(GF(2), x_check_list);
        im_z := VectorSpace(GF(2), z_check_list);

        ker_x := VectorSpace(GF(2), NullspaceMat(TransposedMat(x_check_list)));
        ker_z := VectorSpace(GF(2), NullspaceMat(TransposedMat(z_check_list)));
        
        x_op_list := Filtered(BasisVectors(Basis(ker_x)), v -> not (v in im_x));
        z_op_list := Filtered(BasisVectors(Basis(ker_z)), v -> not (v in im_z));

        return rec(x := x_op_list, z := z_op_list);
    end;

