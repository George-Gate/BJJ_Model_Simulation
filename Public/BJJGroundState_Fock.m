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
J=0;
Ec=-1;
N=20;     % maximum boson number in the space
N0=18;    % boson number of the ground state
clear minPsi minE;

%% init operators
generateFockOperators();
H=sparse(-J*(a1'*a2+a2'*a1)+Ec/8*(a2'*a2-a1'*a1)^2);

%% find ground state
% find an initial state
% only in N0 boson space
Nindex=zeros(N0+1,1);
for i=0:N0            % find the index of N-boson state
    Nindex(i+1)=nn2k(i+1,N0-i+1);
end
% init minE and minPsi
if (~exist('minE','var') || ~exist('minPsi','var'))
    minPsi=zeros(Dim,1);
    minPsi(Nindex)=randPsi(N0+1);
    minE=real(minPsi'*H*minPsi);
end
% use randPsi to find an init state with low energy
for i=1:100
    psi=minPsi;
    psi(Nindex)=psi(Nindex)+randPsi(N0+1);
    psi=psi/sqrt(psi'*psi);
    tmp=real(psi'*H*psi);
    if (tmp<minE)
        minPsi=psi;
        minE=tmp;
    end
end

% use imaginary time evolution method to get ground state
dt=10/max(1,abs(minE));
tol=1e-12;
Udt=sparse(expm(-H*dt));
while(1)
    for i=10/dt:-1:1
        psi_new=Udt*minPsi;    % evolve dt time
        psi_new=psi_new/norm(psi_new);   % renormalize
    end
    if (1-psi_new'*minPsi < tol)     % check whether reached required precision
        minPsi=psi_new;
        % make output state pure real
        minPsi=minPsi*exp(-1i*angle(max(minPsi)));
        minE=real(minPsi'*H*minPsi); % update minE and exit
        break;
    else
        minPsi=psi_new;      % update new state
    end
end
% display ground state
displayFockState(k2nn,minPsi,'|minPsi> = ');

%% plot ground state distribution
plotFockState(minPsi,N0,nn2k);
