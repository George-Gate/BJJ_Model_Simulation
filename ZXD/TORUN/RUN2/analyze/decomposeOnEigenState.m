J=5;
Ec=-1;
N=10;
N0=10;

generateFockOperators;
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
H=sparse(-J*H1+Ec/8*H2);
generateEigenStates;

load()