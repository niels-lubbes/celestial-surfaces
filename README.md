# Celestial surfaces


## Introduction

We present Mathematica code for the automatic verification
of some of the proofs in the article [Self-intersections of surfaces that contain two circles through each point](https://arxiv.org/abs/?).
Below we refer to the environments in this article.

For running the code copy paste the code presented below in [Mathematica](https://www.wolfram.com/mathematica/trial/).
The same code can be found in the Mathematica files
[dP-surface-lattice.nb](https://github.com/niels-lubbes/celestial-surfaces/blob/master/dP-surface-lattice.nb)
and
[parametric-type.nb](https://github.com/niels-lubbes/celestial-surfaces/blob/master/parametric-type.nb).


## Parametric type

The following code verifies that the map associated to a parametric type is
the map defined at [Remark 9](https://arxiv.org/abs/?).

```Mathematica
Remove["Global`*"]

(* projectivized Hamiltonian product *)
hp[a_, b_] := {
   a[[1]]*b[[1]],
   a[[2]]*b[[2]] - a[[3]]*b[[3]] - a[[4]]*b[[4]] - a[[5]]*b[[5]],
   a[[2]]*b[[3]] + a[[3]]*b[[2]] + a[[4]]*b[[5]] - a[[5]]*b[[4]],
   a[[2]]*b[[4]] + a[[4]]*b[[2]] + a[[5]]*b[[3]] - a[[3]]*b[[5]],
   a[[2]]*b[[5]] + a[[5]]*b[[2]] + a[[3]]*b[[4]] - a[[4]]*b[[3]] };

(* inverse stereographic projection *)
isp[p_] := {
   p[[2]]^2 + p[[3]]^2 + p[[4]]^2 + p[[1]]^2,
   2*p[[1]]*p[[2]],
   2*p[[1]]*p[[3]],
   2*p[[1]]*p[[4]],
   p[[2]]^2 + p[[3]]^2 + p[[4]]^2 - p[[1]]^2 };

(* stereographic projection *)
sp[q_] := {q[[1]] - q[[5]], q[[2]], q[[3]], q[[4]]};

(* dehomogenization map from P3 into R3 *)
dh[p_] := {p[[2]], p[[3]], p[[4]]}/p[[1]];
```

We construct the maps `psi1` and `psi2` in two different ways.

```Mathematica

(* Hamiltonian product of circles *)
A = isp@{1, t0*Cos[a] + t1, t0*Sin[a] + t2, t3};
B = isp@{1, u0*Cos[b] + u1, u0*Sin[b] + u2, u3};
psi1 = dh@sp@hp[A, B];

(* encoding of parametric type *)
paramType2Map[t0_, t1_, t2_, t3_, u0_, u1_, u2_, u3_] := Module[
    {aa, bb, cc, dd, ee, ff, gg, hh, qq},
    aa = t1 + t0*Cos[a]; bb = t2 + t0*Sin[a]; cc = t3;
    dd = (aa^2 + bb^2 + cc^2 - 1)/2;
    ee = u1 + u0*Cos[b]; ff = u2 + u0*Sin[b]; gg = u3;
    hh = (ee^2 + ff^2 + gg^2 - 1)/2;
    qq = -aa*hh - bb*gg + cc*ff - dd*ee + (dd + 1)*(hh + 1);
    Return[{aa*ee - bb*ff - cc*gg - dd*hh,
            aa*ff + bb*ee + cc*hh - dd*gg,
            aa*gg - bb*hh + cc*ee + dd*ff}/qq ];
];
psi2 = paramType2Map[t0, t1, t2, t3, u0, u1, u2, u3];
```

We verify that the numerator and denominators
of the entries of `psi1` and `psi2` are the equal up to constant factor 4,
and thus the rational maps they represent are equal.

```Mathematica
Table[ Expand@Numerator@psi1[[i]] == Expand[4*Numerator@psi2[[i]]], {i, 3}]
Table[ Expand@Denominator@psi1[[i]] == Expand[4*Denominator@psi2[[i]]], {i, 1, 3}]
```

Output:

    {True, True, True}
    {True, True, True}


## Experiment with Cliffordian surfaces

We initialize the parametric types in [Table 2](https://arxiv.org/abs/?).

```Mathematica
D3a = {   1,  -1, 0, 0}~Join~{  1,   0, 3/2, 0};
D3b = {1/10, 1/2, 0, 0}~Join~{1/2,   0, 1/2, 0};
D3c = {3/10, 1/2, 0, 0}~Join~{1/2,   0, 1/2, 0};
D3d = {1/10,   0, 0,-4}~Join~{  1,   0,   0,-1};
D3e = { 1/5, 1/2, 0, 0}~Join~{1/2,   0, 1/2, 0};

D4a = {   1,   1, 0, 0}~Join~{  1,   0,   1, 0};
D4b = { 1/2,   0, 1, 0}~Join~{1/2,   1,   0, 0};
D4c = { 1/2, 1/2, 0, 0}~Join~{1/2,   0, 1/2, 0};
D4d = {   1,-3/2, 0, 0}~Join~{  1,   0,-3/2, 0};

D5a = {   1,   0, 0, 0}~Join~{  1, 3/2,   0, 0};
D5b = {   1,   0, 0, 0}~Join~{  1,   2,   0, 0};
D5c = {   1,   0, 0, 0}~Join~{  1, 5/2,   0, 0};
```

Plot surface with parametric type `D5a`

```Mathematica
{S1, T1X, T1Y, T1Z, S2, T2X, T2Y, T2Z} = D5a

ParametricPlot3D[
    paramType2Map[S1, T1X, T1Y, T1Z, S2, T2X, T2Y, T2Z], {a, 0, 2 Pi}, {b, 0, 2 Pi},
    PlotRange -> All,
    Mesh -> False,
    PlotPoints -> 50,
    PlotStyle -> Directive[Opacity[0.8]],
    Boxed -> False,
    Axes -> False
]
```

Output:

    {1, 0, 0, 0, 1, 3/2, 0, 0}

![output image](https://raw.githubusercontent.com/niels-lubbes/celestial-surfaces/master/D5a.png?token=GHSAT0AAAAAACCLVZOYCSWJRK5RF7HMFOZMZCY6PDA "Cliffordian surface D5a")

See the Mathematica file
[parametric-type.nb](https://raw.githubusercontent.com/niels-lubbes/cyclides/master/cyclides.nb?token=AF74RI5WHEIOYDZWOJAIARTBUUO6G)
for additional code that enables to construct Cliffordian surfaces by adjusting sliders.


## Initialization of classes  in the Neron-Severi lattice

We follow the notations and definitions at [Section 2](https://arxiv.org/abs/?).
The code below encodes elements in the Neron-Severi lattice N(X) of a dP surface
in terms of a list such that the class `a0*l0+a1*l1+c1*e1+...+c4*e4`
is represented as the list `{a0,a1,c1,c2,c3,c4}`.

```Mathematica
Remove["Global`*"]

(* The intersection product of classes u and v is computed as u.M.v. *)
M = { {0, 1,  0,  0,  0,  0},
      {1, 0,  0,  0,  0,  0},
      {0, 0, -1,  0,  0,  0},
      {0, 0,  0, -1,  0,  0},
      {0, 0,  0,  0, -1,  0},
      {0, 0,  0,  0,  0, -1} };

ak = {2, 2, -1, -1, -1, -1}; (* anticanonical class *)
l0 = {1, 0,  0,  0,  0,  0};
l1 = {0, 1,  0,  0,  0,  0};
```

We construct all possible elements for B(X), E(X) and G(X).

```
perm = Permutations[Array[Mod[#, 4] - 1 &, 4*6], {6}];
blist = Select[perm, #.M.# == -2 && #.M.ak == 0 && #.M.l0 >= 0 && #.M.l1 >= 0 &]
glist = Select[perm, #.M.# ==  0 && #.M.ak == 2 &]
elist = Select[perm, #.M.# == -1 && #.M.ak == 1 &]
```

Output

    { {0, 1,-1, 0,-1, 0}, {0, 1,-1, 0, 0,-1}, {0, 1,-1,-1, 0, 0}, {0, 1, 0,-1,-1, 0},
      {0, 1, 0,-1, 0,-1}, {0, 1, 0, 0,-1,-1}, {0, 0, 1,-1, 0, 0}, {0, 0, 1, 0,-1, 0},
      {0, 0, 1, 0, 0,-1}, {0, 0,-1, 1, 0, 0}, {0, 0,-1, 0, 1, 0}, {0, 0,-1, 0, 0, 1},
      {0, 0, 0, 1,-1, 0}, {0, 0, 0, 1, 0,-1}, {0, 0, 0,-1, 1, 0}, {0, 0, 0,-1, 0, 1},
      {0, 0, 0, 0, 1,-1}, {0, 0, 0, 0,-1, 1}, {1, 0,-1, 0,-1, 0}, {1, 0,-1, 0, 0,-1},
      {1, 0,-1,-1, 0, 0}, {1, 0, 0,-1,-1, 0}, {1, 0, 0,-1, 0,-1}, {1, 0, 0, 0,-1,-1},
      {1, 1,-1,-1,-1,-1} }

    { {0, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0}, {1, 2,-1,-1,-1,-1}, {1, 1, 0,-1, 0,-1},
      {1, 1, 0,-1,-1, 0}, {1, 1, 0, 0,-1,-1}, {1, 1,-1, 0, 0,-1}, {1, 1,-1, 0,-1, 0},
      {1, 1,-1,-1, 0, 0}, {2, 1,-1,-1,-1,-1} }

    { {0, 1,-1, 0, 0, 0}, {0, 1, 0,-1, 0, 0}, {0, 1, 0, 0,-1, 0}, {0, 1, 0, 0, 0,-1},
      {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1},
      {1, 0,-1, 0, 0, 0}, {1, 0, 0,-1, 0, 0}, {1, 0, 0, 0,-1, 0}, {1, 0, 0, 0, 0,-1},
      {1, 1, 0,-1,-1,-1}, {1, 1,-1, 0,-1,-1}, {1, 1,-1,-1, 0,-1}, {1, 1,-1,-1,-1, 0} }


We define shorthand notation for the classes in the sets B(X), E(X) and G(X).

```Mathematica
(* declare classes in E(X) *)
e1={0,0,1,0,0,0};e2={0,0,0,1,0,0};e3={0,0,0,0,1,0};e4={0,0,0,0,0,1};
e01={1,0,-1,0,0,0};e02={1,0,0,-1,0,0};e03={1,0,0,0,-1,0};e04={1,0,0,0,0,-1};
e11={0,1,-1,0,0,0};e12={0,1,0,-1,0,0};e13={0,1,0,0,-1,0};e14={0,1,0,0,0,-1};
ep1={1,1,0,-1,-1,-1};ep2={1,1,-1,0,-1,-1};ep3={1,1,-1,-1,0,-1};ep4={1,1,-1,-1,-1,0};
(* declare classes in G(X) *)
g0={1,0,0,0,0,0};g1={0,1,0,0,0,0};g2={2,1,-1,-1,-1,-1};g3={1,2,-1,-1,-1,-1};
g12={1,1,-1,-1,0,0};g13={1,1,-1,0,-1,0};g14={1,1,-1,0,0,-1};
g23={1,1,0,-1,-1,0};g24={1,1,0,-1,0,-1};g34={1,1,0,0,-1,-1};
(* declare classes in B(X) *)
b12={1,0,-1,-1,0,0};b13={1,0,-1,0,-1,0};b14={1,0,-1,0,0,-1};b23={1,0,0,-1,-1,0};
b24={1,0,0,-1,0,-1};b34={1,0,0,0,-1,-1};
bp12={0,1,-1,-1,0,0};bp13={0,1,-1,0,-1,0};bp14={0,1,-1,0,0,-1};
bp23={0,1,0,-1,-1,0};bp24={0,1,0,-1,0,-1};bp34={0,1,0,0,-1,-1};
b0={1,1,-1,-1,-1,-1};
bt12={0,0,1,-1,0,0};bt13={0,0,1,0,-1,0};bt14={0,0,1,0,0,-1};bt23={0,0,0,1,-1,0};
bt24={0,0,0,1,0,-1};bt34={0,0,0,0,1,-1};
```

The following methods converts a class or real structure to a string.

```Mathematica
str[q_]:=Module[{},
(* E(X) *)
If[q==e1,Return["e1"]];If[q==e2,Return["e2"]];If[q==e3,Return["e3"]];If[q==e4,Return["e4"]];
If[q==e01,Return["e01"]];If[q==e02,Return["e02"]];If[q==e03,Return["e03"]];
If[q==e04,Return["e04"]];If[q==e11,Return["e11"]];If[q==e12,Return["e12"]];
If[q==e13,Return["e13"]];If[q==e14,Return["e14"]];If[q==ep1,Return["ep1"]];
If[q==ep2,Return["ep2"]];If[q==ep3,Return["ep3"]];If[q==ep4,Return["ep4"]];
(* G(X) *)
If[q==g0,Return["g0"]];If[q==g1,Return["g1"]];If[q==g2,Return["g2"]];If[q==g3,Return["g3"]];
If[q==g12,Return["g12"]];If[q==g13,Return["g13"]];If[q==g14,Return["g14"]];
If[q==g23,Return["g23"]];If[q==g24,Return["g24"]];If[q==g34,Return["g34"]];
(* B(X) *)
If[q==b12,Return["b12"]];If[q==b13,Return["b13"]];If[q==b14,Return["b14"]];
If[q==b23,Return["b23"]];If[q==b24,Return["b24"]];If[q==b34,Return["b34"]];
If[q==bp12,Return["bp12"]];If[q==bp13,Return["bp13"]];If[q==bp14,Return["bp14"]];
If[q==bp23,Return["bp23"]];If[q==bp24,Return["bp24"]];If[q==bp34,Return["bp34"]];
If[q==b0,Return["b0"]];If[q==bt12,Return["bt12"]];
If[q==bt13,Return["bt13"]];If[q==bt14,Return["bt14"]];If[q==bt23,Return["bt23"]];
If[q==bt24,Return["bt24"]];If[q==bt34,Return["bt34"]];
(* apply recursively str to each element in the list *)
If[ q!=Flatten[q], Return@ToString[str/@q]];
Return[ToString[q]];
];
```

## Automatic verification for Lemma 15

The following code is for the proof of [Lemma 15](https://arxiv.org/abs/?).

```Mathematica
(*
Construct all possible quartets {a,b,u,v} in B(X) such that:
l0.a>0, l0.b>0,
l1.u>0, l1.v>0,
a.b=0, u.v=0,
a.u!=-1, a.v!=-1.
*)
checkQuartet[q_] := Return[
   l0.M.q[[1]] == l0.M.q[[2]] == 1 &&
   l1.M.q[[3]] == l1.M.q[[4]] == 1 &&
   q[[1]].M.q[[2]] == q[[3]].M.q[[4]] == 0 &&
   q[[1]].M.q[[3]] != -1 && q[[1]].M.q[[4]] != -1
];
quartetList = Select[Permutations[blist, {4}], checkQuartet[#] &];

(*
Identify elements in quartetList that are equivalent up to permuting {e1,e2,e3,e4}.
*)
doPerm[c_,pm_]:=Return[{c[[1]],c[[2]],c[[pm[[1]]]],c[[pm[[2]]]],c[[pm[[3]]]],c[[pm[[4]]]]}];

checkQuartetEq[a_, b_] := Module[{pmList, i},
    If[{a[[2]], a[[1]], a[[3]], a[[4]]} == b, Return[True]];
    pmList = Permutations[{3, 4, 5, 6}, {4}];
    For[i = 1, i < Length[pmList], i++,
        If[doPerm[#, pmList[[i]]] & /@ a == b, Return[True]];
    ];
    Return[False];
];

reducedQuartetList = DeleteDuplicates[quartetList, checkQuartetEq[#1, #2] &]
```

Output

    {{{0, 1, -1, 0, -1, 0}, {0, 1, 0, -1, 0, -1}, {1, 0, -1, 0, 0, -1}, {1, 0, 0, -1, -1, 0}}}

From the quartet in `reducedQuartetList`, we construct B(X), G(X) and E(X).

```Mathematica
{q1, q2, q3, q4} = reducedQuartetList[[1]];
BX4 = Select[ blist, #.M.q1 != -1 && #.M.q2 != -1 && #.M.q3 != -1 && #.M.q4 != -1 &]
GX4 = Select[ glist, #.M.q1 >=  0 && #.M.q2 >=  0 && #.M.q3 >=  0 && #.M.q4 >=  0 &]
EX4 = Select[ elist, #.M.q1 >=  0 && #.M.q2 >=  0 && #.M.q3 >=  0 && #.M.q4 >=  0 &]
```

Output

    {{0, 1, -1, 0, -1, 0}, {0, 1, 0, -1, 0, -1}, {1, 0, -1, 0, 0, -1}, {1, 0, 0, -1, -1, 0}}
    {{0, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, -1, -1}, {1,  1, -1, -1, 0, 0}}
    {{0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}}


Conversion of output to shorthand notation.

```Mathematica
str@BX4
str@GX4
str@EX4
```

Output

    {bp13, bp24, b14, b23}
    {g1, g0, g34, g12}
    {e1, e2, e3, e4}


## Intersection numbers for Lemma 21

We display the intersection products of classes on dP surface of degree 5
as considered in [Lemma 21](https://arxiv.org/abs/?).

```Mathematica
allBX5 = {b12, b13, b23, bp12, bp13, bp23, bt12, bt13, bt23};
allGX5 = {g0, g1, g12, g13, g23};
allEX5 = {e1, e2, e3, e01, e02, e03, e11, e12, e13};

(*display table with column and row headers*)
BGE5 = allBX5 ~Join~ allGX5 ~Join~ allEX5;
tab = Table[ BGE5[[i]].M.BGE5[[j]], {i, 1, Length[BGE5]}, {j, 1, Length[BGE5]}];
TableForm@({{x} ~Join~ (str/@BGE5)} ~Join~ Transpose[{str/@BGE5} ~Join~ Transpose@tab])
```

Output:

    x    b12 b13 b23 bp12 bp13 bp23 bt12 bt13 bt23 g0 g1 g12 g13 g23 e1 e2 e3 e01 e02 e03 e11 e12 e13
    b12  -2  -1  -1  -1    0    0    0    1    1   0  1  -1   0   0   1  1  0 -1  -1   0   0   0   1
    b13  -1  -2  -1   0   -1    0    1    0   -1   0  1   0  -1   0   1  0  1 -1   0  -1   0   1   0
    b23  -1  -1  -2   0    0   -1   -1   -1    0   0  1   0   0  -1   0  1  1  0  -1  -1   1   0   0
    bp12 -1   0   0  -2   -1   -1    0    1    1   1  0  -1   0   0   1  1  0  0   0   1  -1  -1   0
    bp13  0  -1   0  -1   -2   -1    1    0   -1   1  0   0  -1   0   1  0  1  0   1   0  -1   0  -1
    bp23  0   0  -1  -1   -1   -2   -1   -1    0   1  0   0   0  -1   0  1  1  1   0   0   0  -1  -1
    bt12  0   1  -1   0    1   -1   -2   -1    1   0  0   0   1  -1  -1  1  0  1  -1   0   1  -1   0
    bt13  1   0  -1   1    0   -1   -1   -2   -1   0  0   1   0  -1  -1  0  1  1   0  -1   1   0  -1
    bt23  1  -1   0   1   -1    0    1   -1   -2   0  0   1  -1   0   0 -1  1  0   1  -1   0   1  -1
    g0    0   0   0   1    1    1    0    0    0   0  1   1   1   1   0  0  0  0   0   0   1   1   1
    g1    1   1   1   0    0    0    0    0    0   1  0   1   1   1   0  0  0  1   1   1   0   0   0
    g12  -1   0   0  -1    0    0    0    1    1   1  1   0   1   1   1  1  0  0   0   1   0   0   1
    g13   0  -1   0   0   -1    0    1    0   -1   1  1   1   0   1   1  0  1  0   1   0   0   1   0
    g23   0   0  -1   0    0   -1   -1   -1    0   1  1   1   1   0   0  1  1  1   0   0   1   0   0
    e1    1   1   0   1    1    0   -1   -1    0   0  0   1   1   0  -1  0  0  1   0   0   1   0   0
    e2    1   0   1   1    0    1    1    0   -1   0  0   1   0   1   0 -1  0  0   1   0   0   1   0
    e3    0   1   1   0    1    1    0    1    1   0  0   0   1   1   0  0 -1  0   0   1   0   0   1
    e01  -1  -1   0   0    0    1    1    1    0   0  1   0   0   1   1  0  0 -1   0   0   0   1   1
    e02  -1   0  -1   0    1    0   -1    0    1   0  1   0   1   0   0  1  0  0  -1   0   1   0   1
    e03   0  -1  -1   1    0    0    0   -1   -1   0  1   1   0   0   0  0  1  0   0  -1   1   1   0
    e11   0   0   1  -1   -1    0    1    1    0   1  0   0   0   1   1  0  0  0   1   1  -1   0   0
    e12   0   1   0  -1    0   -1   -1    0    1   1  0   0   1   0   0  1  0  1   0   1   0  -1   0
    e13   1   0   0   0   -1   -1    0   -1   -1   1  0   1   0   0   0  0  1  1   1   0   0   0  -1


## Intersection numbers for Lemma 22

Intersection products of classes on dP surface of degree 6
as considered in [Lemma 22](https://arxiv.org/abs/?)

```Mathematica
allBX6 = {b12, bp12, bt12};
allGX6 = {g0, g1, g12};
allEX6 = {e1, e2, e01, e02, e11, e12};

(*display table with column and row headers*)
BGE6 = allBX6 ~Join~ allGX6 ~Join~ allEX6;
tab = Table[ BGE6[[i]].M.BGE6[[j]], {i, 1, Length[BGE6]}, {j, 1, Length[BGE6]}];
TableForm@({{x} ~Join~ (str /@ BGE6)} ~Join~ Transpose[{str /@ BGE6} ~Join~ Transpose@tab])
```

Output

    x     b12 bp12 bt12 g0 g1 g12 e1 e2 e01 e02 e11 e12
    b12   -2  -1    0   0  1  -1   1  1 -1  -1   0   0
    bp12  -1  -2    0   1  0  -1   1  1  0   0  -1  -1
    bt12   0   0   -2   0  0   0  -1  1  1  -1   1  -1
    g0     0   1    0   0  1   1   0  0  0   0   1   1
    g1     1   0    0   1  0   1   0  0  1   1   0   0
    g12   -1  -1    0   1  1   0   1  1  0   0   0   0
    e1     1   1   -1   0  0   1  -1  0  1   0   1   0
    e2     1   1    1   0  0   1   0 -1  0   1   0   1
    e01   -1   0    1   0  1   0   1  0 -1   0   0   1
    e02   -1   0   -1   0  1   0   0  1  0  -1   1   0
    e11    0  -1    1   1  0   0   1  0  0   1  -1   0
    e12    0  -1   -1   1  0   0   0  1  1   0   0  -1


## Intersection numbers for Proposition 33

Intersection products of classes on dP surface of degree 4 such that B(X)
has four components (see [Proposition 33](https://arxiv.org/abs/?)).

```Mathematica
BGE4 = BX4 ~Join~ GX4 ~Join~ EX4;
tab = Table[ BGE4[[i]].M.BGE4[[j]], {i, 1, Length[BGE4]}, {j, 1, Length[BGE4]}];
TableForm@({{x} ~Join~ (str/@ BGE4)} ~Join~ Transpose[{str/@ BGE4} ~Join~ Transpose@tab])
```

Output

    x    bp13 bp24 b14 b23 g1 g0 g34 g12 e1 e2 e3 e4
    bp13 -2    0    0   0  0  1  0   0    1  0  1  0
    bp24  0   -2    0   0  0  1  0   0    0  1  0  1
    b14   0    0   -2   0  1  0  0   0    1  0  0  1
    b23   0    0    0  -2  1  0  0   0    0  1  1  0
    g1    0    0    1   1  0  1  1   1    0  0  0  0
    g0    1    1    0   0  1  0  1   1    0  0  0  0
    g34   0    0    0   0  1  1  0   2    0  0  1  1
    g12   0    0    0   0  1  1  2   0    1  1  0  0
    e1    1    0    1   0  0  0  0   1   -1  0  0  0
    e2    0    1    0   1  0  0  0   1    0 -1  0  0
    e3    1    0    0   1  0  0  1   0    0  0 -1  0
    e4    0    1    1   0  0  0  1   0    0  0  0 -1

