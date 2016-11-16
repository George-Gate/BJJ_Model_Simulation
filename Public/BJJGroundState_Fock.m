% Use Fock representation
% [Encoding] 
%-------------------------------------------------------------------------------
%    [1,0,0,....] -> |0,0>  [0,1,0,...] -> |1,0> ... [0,...,0,1,0,...] -> |N,0>
%    [0,...,0,1,0,...] -> |0,1>  [0,...,0,1,0,...] -> |1,1> ... [...] -> |N-1,1>
%       (N+1 leading 0) 
%    [0,...,0,1,0,...] -> |0,2>  ...    [...] -> |N-2,2>
%       (2N+1 leading 0)
%    ......
%    [0,...,0,1,0,...] -> |0,k>  ...    [...] -> |N-k,k>    
%       ( k*N+(3k-k^2)/2 leading 0)
%-------------------------------------------------------------------------------
% Total Dim: (N+1)*(N+2)/2
% |n1,n2> --> psi( n2*N+(3*n2-n2^2)/2+n1+1 )=1   (n1+n2<=N)
%             psi( nn2k(n1+1,n2+1)  )=1
% psi( k )=1  --> |k2nn(1),k2nn(2)>
%
% Under some circumstance, a1 and a2' may not commute. For example, |k,N-k>,
% if we do a2' first, |k,N-k+1> has N+1 boson which is not included in our
% numerical representation. But if we do a1 first and then do a2', no error
% will occur. To avoid such circumstance, try to start the simulation from 
% N-1 bosons. Or always act annihilation before creation. 
%
%% parameters
J=1000;
Ec=-1;
N=12;     % maximum boson number in the space
N0=12;    % boson number of the ground state

%% init operators
generateFockOperators();
H=sparse(-J*(a1'*a2+a2'*a1)+Ec/8*(a2'*a2-a1'*a1)^2);

%% find ground state
minPsi=makeState('BJJ Ground','psi',nn2k,N0,Dim,H);

psiSCS=makeState('SCS','psi',nn2k,N0,Dim,J);
display(['Fedelity with SCS=',num2str(abs(minPsi'*psiSCS),5)]);

psiSCS'*H*psiSCS


%% plot ground state distribution
plotFockState(minPsi,N0,nn2k);
% display ground state
displayFockState(k2nn,minPsi,'|minPsi> = ');
