(* ::Package:: *)

BeginPackage["mewsAndAlphas`"]

mewsAndAlphas::usage=
	"This package conatins all the mu's and alpha definitions"
	
muRho=(m1*m2)/(m1 + m2);
muLambda=((m1 + m2)*m3)/(m1 + m2 + m3);
mu1=(m2*m3)/(m2 + m3);
mu2=(m1*m3)/(m1 + m3);
{{a1, b1}, {c1, d1}} = {{-(m3/(m2 + m3)), -1}, {(m2*(m1 + m2 + m3))/((m1 + m2)*(m2 + m3)), -(m1/(m1 + m2))}};
{{a2, b2}, {c2, d2}} = {{-(m3/(m1 + m3)), 1}, {((-m1)*(m1 + m2 + m3))/((m1 + m2)*(m1 + m3)), -(m2/(m1 + m2))}};
alpha = Sqrt[muRho*omega];
beta = Sqrt[muLambda*omega];
alpha1 = (alpha*a1)/Cos[theta[1]];
alpha2 = (alpha*a2)/Cos[theta[2]];
alpha11 = (beta*c1)/Sin[theta[1]];
alpha22 = (beta*c2)/Sin[theta[2]];
beta1 = -((alpha*b1)/Sin[theta[1]]);
beta11 = (beta*d1)/Cos[theta[1]];
beta2 = -((alpha*b2)/Sin[theta[2]]);
beta22 = (beta*d2)/Cos[theta[2]];
alpha1 == alpha11;
alpha2 == alpha22;
beta1 == beta11;
beta2 == beta22;

EndPackage[]



