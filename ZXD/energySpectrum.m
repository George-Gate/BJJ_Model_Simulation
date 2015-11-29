% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
% This script will calc the energy spectrum of BJJ Hamiltonian, use J as a
% parameter.

%% Parameters
N=10;
J1=0;J2=5;     % J will change from J1 to J2
Ec1=-1;Ec2=-1; % Ec will change from Ec1 to Ec2, linearly as J changes.
spCount=3000;   % The total number of sample point

generateFockOperators;

%% Calc the eigen value & state of Hamiltonian
J=linspace(J1,J2,spCount);
Ec=linspace(Ec1,Ec2,spCount);

% calc constant operators
H1=full(a1'*a2+a2'*a1);
H2=full((a2'*a2-a1'*a1)^2);
particleNum=full(a1'*a1+a2'*a2);
%sqParticleNum=particleNum*particleNum;

energy=zeros(Dim,spCount);   % eigen energy
exp_n=zeros(Dim,spCount);    % expected total particle number of eigen state
%dev_n=zeros(Dim,spCount);    % variance of total particle number

for i=1:spCount
    H=-J(i)*H1+Ec(i)/8*H2;
    % calc eigen value & state
    [V,D]=eig(H);
    % record eigen value
    energy(:,i)=diag(D);
    % analyze eigen state
    for j=1:Dim
        exp_n(j,i)=V(:,j)'*particleNum*V(:,j);
        %dev_n(j,i)=V(:,j)'*sqParticleNum*V(:,j)-exp_n(j,i)^2;
    end
end

% Select the eigen states that have the desired total particle number
% It is checked that all energy eigen states are also an eigen state of
% total particle number.
energyToDraw=energy;
energyToDraw(round(exp_n)<N)=NaN;

%% plot
for i=1:Dim
    scatter(J,energyToDraw(i,:),1,'red');
    hold on;
end
hold off;box on;
xlabel('J');
ylabel('energy');
title(['Energy Spectrum (N=', num2str(N),' Ec from ', num2str(Ec1) ' to ', num2str(Ec2),')']);
