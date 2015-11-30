% Generate the energy eigen state for current Hamiltonian
% Only generate eigen state that has total particle number N0
% Result save in V and D
% 
% Requirement:
%   H - Current Hamiltonian
%   N0 - The total particle number of the desired eigen state
%   a1, a2

[V,D]=eig(full(H));
D=diag(D);

% ignore states with total particle number neq to N0
particleNum=full(a1'*a1+a2'*a2);
for j=1:Dim
    exp_n=V(:,j)'*particleNum*V(:,j);
    if (round(exp_n)~=N0)
        D(j)=NaN;
    end
end

% sort eigen states
[D,Index]=sort(D);
V=V(:,Index);
V=V(:,1:N0+1);
D=D(1:N0+1);
clear Index;