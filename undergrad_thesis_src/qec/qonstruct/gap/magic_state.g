# Author: Suhas Vittal
# date:    30 January 2024


PrintMatrix :=
    function(M)
        local i, j;
        
        i := 1;
        while i <= Length(M) do
            j := 1;
            while j <= Length(M[i]) do
                if M[i][j] = One(GF(2)) then
                    Print("1 ");
                else
                    Print("0 ");
                fi;
                j := j+1;
            od;
            i := i+1;
            Print("\n");
        od;
    end;

M := [[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
      [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1],
      [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
      [1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
      [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]];
M := M * One(GF(2));
ech := SemiEchelonMat(M);
A := ech.vectors;
Print(ech.heads, "\n");

PrintMatrix(A);
