#
# Basic operations 
#

a := 2+5;
b := (a=7);
c := (a=7);     # IsIdenticalObj(b, c) = false
d := b;         # IsIdenticalObj(b, d) = true

#
# Lists
#

primes := [2, 3, 5, 7, 9];  # Note: these are 1-indexed.
Add(primes, 11);
Append(primes, [13, 17]);

primes_immut := Immutable(primes);      # Cannot be modified.
position_of_3 := Position(primes, 3);   # Get's the index of the element
n_primes := Length(primes);             # Note: GAP handles resizing.

#
# Sets
#

z_mod_4_set := Set([0, 1, 2, 3]);   # Lists and Sets are functionally
                                    # identicable, except Sets are sorted.
z_mod_4_set_is_a_set := IsSSortedList(z_mod_4_set); 
z_mod_5_set := Set([0, 1, 2, 3]);
AddSet(z_mod_5_set, 4);

#
# Ranges
#

zero_to_hundred := [0..100]; # Also functionally identical to Lists.

#
# Loops
# 

sum_of_first_four_ints := 0;
for x in [1..4] do
    sum_of_first_four_ints := sum_of_first_four_ints + x;
od;

# While loops are similar: while (cond) do ... od;
# Vector operations:
#   (1) Product(List)
#   (2) Sum(List)
#   (3) List(List, function)
#   (4) Filtered(List, function)
#   (5) List{ [indices] } --> gets sublist

#
# Functions
# 

multiply_acc := 
    function(a, b, c)
        return a*b + c;
    end;
