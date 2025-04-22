#
# Author:   Suhas Vittal
# date:     19 January 2024
#
# Useful operators.
#

Any :=
    function(li)
        local x;
        for x in li do
            if x then return true; fi;
        od;
        return false;
    end;

All :=
    function(li)
        local x;
        for x in li do
            if not x then return false; fi;
        od;
        return true;
    end;
