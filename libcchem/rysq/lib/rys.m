(*gx[-2, 0] = 0
gx[-1, 0] = 0
gx[0, j_:0] = 1
gx[i_, j_:0] := (i-1) (B t[a]) gx[i-2, j] + (rx - Rx t[a]) gx[i-1, j]

gy[-2, 0] = 0
gy[-1, 0] = 0
gy[0, j_:0] = 1
gy[i_, j_:0] := (i-1) (B t[a]) gy[i-2, j] + (ry - Ry t[a]) gy[i-1, j]

gz[-2, 0] = 0
gz[-1, 0] = 0
gz[0, j_:0] = W[a]
gz[i_, j_:0] := (i-1) (B t[a]) gz[i-2, j] + (rz - Rz t[a]) gz[i-1, j]
*)

gx[-2, i_] := 0
gx[-1, i_] := 0
gx[i_, -1] := 0
gx[i_, -2] := 0
gx[0, 0] := W[a]

B00 := 0.5 t[a]
B10 := 1/(2 A) (1 - B t[a])
B01 := 1/(2 B) (1 - A t[a])

Cx10 := xAi - B xAB t[a]
Cx01 := xBk + A xAB t[a]  

gx[i_, 0] := (i-1) B10 gx[i-2,0] + 0 B00 gx[i-1,0-1] + Cx10 gx[i-1,0]
gx[i_, j_:1] := (j-1) B01 gx[i,j-2] + i B00 gx[i-1,j-1] + Cx01 gx[i,j-1]


gy[-2, i_] := 0
gy[-1, i_] := 0
gy[i_, -1] := 0
gy[i_, -2] := 0
gy[0, 0] := 1

Cy10 := yAi - B yAB t[a]
Cy01 := yBk + A yAB t[a]  

gy[i_, 0] := (i-1) B10 gy[i-2,0] + 0 B00 gy[i-1,0-1] + Cy10 gy[i-1,0]
gy[i_, j_:1] := (j-1) B01 gy[i,j-2] + i B00 gy[i-1,j-1] + Cy01 gy[i,j-1]

gz[-2, i_] := 0
gz[-1, i_] := 0
gz[i_, -1] := 0
gz[i_, -2] := 0
gz[0, 0] := 1

Cz10 := zAi - B zAB t[a]
Cz01 := zBk + A zAB t[a]  

gz[i_, 0] := (i-1) B10 gz[i-2,0] + 0 B00 gz[i-1,0-1] + Cz10 gz[i-1,0]
gz[i_, j_:1] := (j-1) B01 gz[i,j-2] + i B00 gz[i-1,j-1] + Cz01 gz[i,j-1]

    ix [i_, 0] := gx [i,0]
    ix [i_,j_:0]:= dx ix [i +1,j -1] + ix [i,j -1]


    iy [i_, 0] := gy [i,0]
    iy [i_,j_:0]:= dy iy [i +1,j -1] + iy [i,j -1]

    iz [i_, 0] := gz [i,0]
    iz [i_,j_:0]:= dz iz [i +1,j -1] + iz [i,j -1]
