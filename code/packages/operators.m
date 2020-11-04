(* ::Package:: *)

BeginPackage["operators`"]

operators::usage=
"This package stores all operators: kineticOperator[Lambda,gammaMax,omega], polynomialOperator[Lambda,gammaMax,power]"

Begin["`Private`"]

(* KINETIC OPERATOR MATRIX MAKER 1 SIZE FITS ALL *)
kineticOperator[Lambda_,gammaMax_,omega_]:=(
sM=stateMatrix[Lambda,gammaMax];
new={};
For[i=1,i<=Length[sM],i++,
For[j=1,j<=Length[sM],j++,
np=sM[[i,j,1,1]];
lp=sM[[i,j,1,2]];
n=sM[[i,j,2,1]];
l=sM[[i,j,2,2]];
nlp=sM[[i,j,1,3]];
llp=sM[[i,j,1,4]];
nl=sM[[i,j,2,3]];
ll=sM[[i,j,2,4]];
qhoRho=KroneckerDelta[np,n]*KroneckerDelta[lp,l]*KroneckerDelta[nlp,nl]*KroneckerDelta[llp,ll]*omega(2np+lp+3/2);
qhoLam=KroneckerDelta[np,n]*KroneckerDelta[lp,l]*KroneckerDelta[nlp,nl]*KroneckerDelta[llp,ll]*omega(2nl+ll+3/2);
rhoSqrd=KroneckerDelta[lp,l]*KroneckerDelta[nlp,nl]*KroneckerDelta[llp,ll]*1/2 (muRho omega^2)/alpha^2 rel[np,n,l,2];
lamSqrd=KroneckerDelta[np,n]*KroneckerDelta[lp,l]*KroneckerDelta[llp,ll]*1/2 (muLambda omega^2)/beta^2 rel[nlp,nl,ll,2];
AppendTo[new,qhoRho+qhoLam-rhoSqrd-lamSqrd]
]];
new=ArrayReshape[new,{Length[sM],Length[sM]}]
)

polynomialOperator[Lambda_,gammaMax_,power_]:=(
sM=stateMatrix[Lambda,gammaMax];
new={};
For[i=1,i<=Length[sM],i++,
For[j=1,j<=Length[sM],j++,
np=sM[[i,j,1,1]];
lp=sM[[i,j,1,2]];
n=sM[[i,j,2,1]];
l=sM[[i,j,2,2]];
nlp=sM[[i,j,1,3]];
llp=sM[[i,j,1,4]];
nl=sM[[i,j,2,3]];
ll=sM[[i,j,2,4]];
coulomb=KroneckerDelta[lp,l]*KroneckerDelta[nlp,nl]*KroneckerDelta[llp,ll]* rel[np,n,l,power];
AppendTo[new,coulomb]
]];
new=ArrayReshape[new,{Length[sM],Length[sM]}]
)

End[]

EndPackage[]
