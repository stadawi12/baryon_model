(* ::Package:: *)

BeginPackage["oscillatorBracket`"]

oscillatorBracket::usage=
"This contains the oscillator bracket function and all necessary functions to compute matrix representations of operators,
angle, nineJ, talmi, sumArgs, talmiCoeffs, tBlock, stateMatrix, tMatrix, and all the integrals."

Begin["`Private`"]

ClearAll[angle]
angle[system_]:=(
If[system == 1, tan = -Sqrt[(m2*(m1 + m2 + m3))/(m1*m3)], tan = Sqrt[(m1*(m1 + m2 + m3))/(m2*m3)]]; 
Return[ArcTan[tan] + Pi])

(* Nine-J function *)
nineJ[j1_, j2_, j4_, j5_, x_, L_, l_, Lambda_, l1_, l2_]:=(-1)^(2*x)*(2*x + 1)*SixJSymbol[{j1, j4, L}, {l, Lambda, x}]*SixJSymbol[{j2, j5, l}, {j4, x, l2}]*SixJSymbol[{l1, l2, Lambda}, {x, j1, j2}]

(* Talmi coefficient function *)
talmi[n1_, l1_, n2_, l2_, n_, l_, nN_, L_, Lambda_, system_] := (
   fMatrix = {{2, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 2, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 2, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 2, 1}}; 
   ListOfFs = fMatrix.{n1, l1, n2, l2, nN, L, n, l}; 
   aMatrix = {{2, 2, 0, 0, 1, 1, 0, 0}, {0, 0, 2, 2, 0, 0, 1, 1}, {2, 0, 2, 0, 1, 0, 1, 0}, {0, 2, 0, 2, 0, 1, 0, 1}}; 
   (* determine amin, bmin, cmin and dmin *)
   If[ListOfFs[[1]] < ListOfFs[[3]], amin = ListOfFs[[1]], amin = ListOfFs[[3]]]; 
   If[ListOfFs[[1]] < ListOfFs[[4]], bmin = ListOfFs[[1]], bmin = ListOfFs[[4]]]; 
   If[ListOfFs[[2]] < ListOfFs[[3]], cmin = ListOfFs[[2]], cmin = ListOfFs[[3]]]; 
   If[ListOfFs[[2]] < ListOfFs[[4]], dmin = ListOfFs[[2]], dmin = ListOfFs[[4]]]; 
   SummationVariables = {}; 
   (* Obtain summation variables {(a,b,c,d,la,lb,lc,ld)} *)
   For[a = 0, a < Floor[amin/2] + 1, a++, 
   For[b = 0, b < Floor[bmin/2] + 1, b++, 
   For[c = 0, c < Floor[cmin/2] + 1, c++, 
   For[d = 0, d < Floor[dmin/2] + 1, d++, 
   For[la = 0, la < amin + 1, la++, 
   For[lb = 0, lb < bmin + 1, lb++, 
   For[lc = 0, lc < cmin + 1, lc++, 
   For[ld = 0, ld < dmin + 1, ld++, 
   If[aMatrix . {a, b, c, d, la, lb, lc, ld} == ListOfFs, 
   If[MemberQ[Range[Abs[la - lc], la + lc], L], 
   If[MemberQ[Range[Abs[lb - ld], lb + ld], l], 
   If[MemberQ[Range[Abs[la - lb], la + lb], l1], 
   If[MemberQ[Range[Abs[lc - ld], lc + ld], l2], 
   If[EvenQ[la + lc + L], 
   If[EvenQ[lb + ld + l], 
   If[EvenQ[la + lb + l1], 
   If[EvenQ[lc + ld + l2], 
   (* bin for collecting sum terms *)
   AppendTo[SummationVariables, {a, b, c, d, la, lb, lc, ld}]]]]]]]]]]]]]]]]]]; 
   sum = {}; 
   
   coeff = (I^(l1 + l2 + L + l)*Sqrt[n1!*n2!*nN!*n!*(2*(n1 + l1) + 1)!!*(2*(n2 + l2) + 1)!!*(2*(nN + L) + 1)!!*(2*(n + l) + 1)!!])/2^((l1 + l2 + L + l)/4);
   (* calculate sum terms and append them to bin (sum) *)
   For[i = 1, i < Length[SummationVariables] + 1, i++, 
     J = {SummationVariables[[i,5]], SummationVariables[[i,6]], l1, SummationVariables[[i,7]], SummationVariables[[i,8]], l2, L, l, \[CapitalLambda]}; 
     a = SummationVariables[[i,1]];
     b = SummationVariables[[i,2]]; 
     c = SummationVariables[[i,3]]; 
     d = SummationVariables[[i,4]]; 
     la = J[[1]]; 
     lb = J[[2]]; 
     lc = J[[4]]; 
     ld = J[[5]]; 
     y = Intersection[Range[Abs[la - Lambda], la + Lambda], Range[Abs[l - lc], l + lc], Range[Abs[lb - l2], lb + l2]]; 
     (* Append the crazy formula into sum *)
     AppendTo[sum, (-1)^(la + lb + lc)*2^((la + lb + lc + ld)/2)*Sin[angle[system]]^(2*a + la + 2*d + ld)*
       Cos[angle[system]]^(2*b + lb + 2*c + lc)*Sum[nineJ[la, lb, lc, ld, y[[i]], L, l, Lambda, l1, l2], {i, 1, Length[y]}]*(((2*la + 1)*(2*lb + 1)*(2*lc + 1)*(2*ld + 1))/
        (a!*b!*c!*d!*(2*(a + la) + 1)!!*(2*(b + lb) + 1)!!*(2*(c + lc) + 1)!!*(2*(d + ld) + 1)!!))*ClebschGordan[{la, 0}, {lc, 0}, {L, 0}]*ClebschGordan[{lb, 0}, {ld, 0}, {l, 0}]*
       ClebschGordan[{la, 0}, {lb, 0}, {l1, 0}]*ClebschGordan[{lc, 0}, {ld, 0}, {l2, 0}]]]; 
     (* finally output the result which is the coefficient times the sum constituents all added up *)
       t = (-1)^(l + L - Lambda)*Total[sum]*coeff; 
       Return[t])
       
(* OBTAINS THE SUMMATION ARGUMENTS (COMBINATIONS OF |n,l,N,L:\[CapitalLambda]> IN THAT ORDER) THAT ARE ALLOWED FOR A SPECIFIED ENERGY \[Gamma] AND T.O.A.M. L*)
sumArgs[Lambda_, gamma_] := (

nlist = {}; 
llist = {}; 

n = 0; While[2*n <= gamma, AppendTo[nlist, n]; n++]; 
l = 0; While[l <= gamma, AppendTo[llist, l]; l++]; 

sArgs = {}; 

   For[i = 1, i <= Length[nlist], i++, 
   For[j = 1, j <= Length[llist], j++, 
   For[k = 1, k <= Length[nlist], k++, 
   For[m = 1, m <= Length[llist], m++, 
       If[2*nlist[[i]] + llist[[j]] + 2*nlist[[k]] + llist[[m]] == gamma, 
       If[MemberQ[Table[lvals, {lvals, Abs[llist[[j]] - llist[[m]]], llist[[j]] + llist[[m]]}], Lambda], 
         AppendTo[sArgs, {nlist[[i]], llist[[j]], nlist[[k]], llist[[m]]}]]]]]]]; 
         Return[sArgs]
)
 
(* OBTAINS COEFFICIENTS OF ALL STATES (system: 1 or 2) *)
talmiCoeffs[n1_, l1_, n2_, l2_, Lambda_, system_]:= (

coeffs = {}; 
gamma = 2*n1 + l1 + 2*n2 + l2; 
args1 = sumArgs[Lambda, gamma]; 

For[wuba = 1, wuba < Length[args1] + 1, wuba++, 
    nN = args1[[wuba,1]]; 
    lL = args1[[wuba,2]]; 
    nn = args1[[wuba,3]]; 
    ll = args1[[wuba,4]]; 
    coeff2 = talmi[n1, l1, n2, l2, nN, lL, nn, ll, Lambda, system]; 
    If[coeff2 != 0, AppendTo[coeffs, {coeff2, args1[[wuba]]}]]
    ]; 
    
   Return[coeffs]
)
 
(* PRODUCES t-MATRIX BLOCKS, CAN BE USED TO PROVE ORTHOGONALITY OF TALMI COEFFICIENTS *)
tBlock[Lambda_,gamma_,system_]:=(

states=sumArgs[Lambda,gamma];

mat=Table[{states[[i]],states[[j]]},{i,1,Length[states]},{j,1,Length[states]}];
bin={};

For[he=1,he<=Length[states],he++,
For[be=1,be<=Length[states],be++,
	nN=mat[[be,he,1,1]];
	L=mat[[be,he,1,2]];
	n=mat[[be,he,1,3]];
	l=mat[[be,he,1,4]];
	n1=mat[[be,he,2,1]];
	l1=mat[[be,he,2,2]];
	n2=mat[[be,he,2,3]];
	l2=mat[[be,he,2,4]];
	AppendTo[bin,talmi[n1,l1,n2,l2,nN,L,n,l,Lambda,system]*1.]
]];
matrix=ArrayReshape[bin,{Length[states],Length[states]}]
)
(*matrix.Transpose[matrix]//FullSimplify//MatrixForm;*)

(* HAMILTONIAN MATRIX ELEMENT LABELER *)
stateMatrix[Lambda_, gammaMax_]:=(
gammaMinEven=If[EvenQ[Lambda], Lambda, Lambda+1];
gammaMinOdd=If[EvenQ[Lambda], Lambda+1, Lambda];
energies=If[EvenQ[gammaMax],Table[n,{n,gammaMinEven,gammaMax,2}],Table[n,{n,gammaMinOdd,gammaMax,2}]];
states=Table[sumArgs[Lambda,energies[[i]]],{i,1,Length[energies]}];
flatStates=Flatten[states,1];
matrixOfStates=Table[{flatStates[[i]],flatStates[[j]]},{i,1,Length[flatStates]},{j,1,Length[flatStates]}]
)

(* USE THIS TO CONSTRUCT THE t-MATRIX *)
tMatrix[Lambda_,gammaMax_,system_]:=(
gammaMinEven=If[EvenQ[Lambda], Lambda, Lambda+1];
gammaMinOdd=If[EvenQ[Lambda], Lambda+1, Lambda];
energies=If[EvenQ[gammaMax],Table[n,{n,gammaMinEven,gammaMax,2}],Table[n,{n,gammaMinOdd,gammaMax,2}]];
blocks=Map[tBlock[Lambda,#,system]&,energies];
DiagonalMatrix[Hold/@blocks]//ReleaseHold//ArrayFlatten
)

ClearAll[Rcoeff,Lcoeff,ksum,rPsum,rel]
Rcoeff[n_,l_]:=Rcoeff[n,l]=Sqrt[(2*n!)/Gamma[n + l + 3/2]] 
Lcoeff[n_,l_,k_]:=Lcoeff[n,l,k]=((-1)^k*Gamma[n + l + 3/2])/((n - k)!*Gamma[k + l + 3/2]*k!)
ksum[n_,m_,l_,p_]:=ksum[n,m ,l,p]=Sum[Lcoeff[n,l,k] Lcoeff[m,l,p-k-l],{k,Max[0,p-l-m],Min[n,p-l]}]
rPsum[n_,m_,l_,y_]:=rPsum[n,m,l,y]=Sum[Gamma[p+(1/2)*y + 3/2] * ksum[n,m,l,p],{p,l,l+n+m}]

rel[n_,m_,l_,y_]:=rel[n,m,l,y]=1/2*Rcoeff[n,l]*Rcoeff[m,l]*rPsum[n,m,l,y]

End[]

EndPackage[]



