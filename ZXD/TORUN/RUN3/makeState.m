%Author: George-Gate
%Date: 2016/01/15
%--------------------------------------------------------------------------
%Generate required state.
%  Function dependency: randPsi -> randU
%
%  state = makeState( stateName, type, nn2k, varargin )
%
%  stateName - The name of the state, can be one of the following:
%              'BJJ Ground', 'SCS', 'Eigen'
%  type - Describe in what format the result should be return. Can be one
%         of the two:
%              'psi': return a state vector (Column vector)
%              'rho': return a density matrix
%  nn2k - Just pass this array
%  varargin - special arguments required by particular state
%
%[Examples]
%  makeState('SCS', 'psi', nn2k, N0, Dim, J)
%    Return a SCS state, pointing to sign(J)*x direction in the Bloch Sphere.
%       N0  - The state has total particle number N0.
%
%  makeState('BJJ Ground', 'rho', nn2k, N0, Dim, H, tol)
%    Return the ground state of given Halmitonian. If the ground state
%    degenerates, the function will return an arbitary combination of those
%    ground states.
%       N0  - The state has total particle number N0.
%       H   - Hamiltonian of the system.
%       tol - [Optional] The precision of iteration.
%
%  makeState('Eigen', 'psi', nn2k, N0, Dim, H, k, opNum)
%    Return a superposition of specific energy eigen states.
%       N0  - The state has total particle number N0.
%       H   - Hamiltonian of the system.
%       k   - A vector, specify how to combine eigen states. The returned
%             state is sum(k(i)*V(:,i)), in which the sum is over i.
%
%             V(:,i) is the energy eigen state with the i-th smallest
%             energy among all states that have total particle number N0.
%
%             The length of k should not excced N0+1, but can be smaller
%             than N0+1, the function will add zeros at the end of k to
%             make it as long as N0+1.
%
%             The function will normalize k first to ensure the
%             normalization of returned state.
%     opNum - Total particle number operator, opNum=a1'*a1+a2'*a2.
%
function state = makeState( stateName, type, nn2k, varargin )
    % generate state vector
    switch stateName
        case 'BJJ Ground'
            % get parameters
            N0=varargin{1};     
            Dim=varargin{2};
            H=varargin{3};
            if (nargin==7)
                tol=varargin{4};
            else
                tol=1e-12;
            end
            % find an initial state only in N0 boson space
            Nindex=zeros(N0+1,1);
            for i=0:N0            % find the index of N-boson state
                Nindex(i+1)=nn2k(i+1,N0-i+1);
            end
            % init minE and minPsi
            minPsi=zeros(Dim,1);
            minPsi(Nindex)=randPsi(N0+1);
            minE=real(minPsi'*H*minPsi);
                
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
            Udt=sparse(expm(-H*dt));
            counter=0;
            while(1)
                counter=counter+1;
                for i=10/dt:-1:1
                    psi_new=Udt*minPsi;    % evolve dt time
                    psi_new=psi_new/norm(psi_new);   % renormalize
                end
                if (1-psi_new'*minPsi < tol || counter>1000)     % check whether reached required precision
                    if (counter>1000)
                        warning(['Can not reach required tolerance, exit because iterate times exceed 1000.',...
                                 ' Current error: ',num2str(abs(1-psi_new'*minPsi),'%6.2e')]);
                    end
                    minPsi=psi_new;
                    % make output state pure real
                    minPsi=minPsi*exp(-1i*angle(max(minPsi)));
                    %minE=real(minPsi'*H*minPsi); % update minE and exit
                    break;
                else
                    minPsi=psi_new;      % update new state
                end
            end
            state=minPsi;
        case 'SCS'
            N0=varargin{1};     % get parameters
            Dim=varargin{2};
            J=varargin{3};
            state=zeros(Dim,1);
            for i=0:N0
                state(nn2k(i+1,N0-i+1))=sign(J)^i*1/sqrt(factorial(i)*factorial(N0-i));
            end
            state=state*sqrt(factorial(N0))/sqrt(2)^N0;
        case 'Eigen'
            N0=varargin{1};     % get parameters
            Dim=varargin{2};
            H=varargin{3};
            k=varargin{4};
            opNum=varargin{5};
            
            % find eigen vectors & values
            [V,D]=eig(full(H));
            D=diag(D);

            % ignore states with total particle number neq to N0
            for j=1:Dim
                exp_n=V(:,j)'*opNum*V(:,j);
                if (round(exp_n)~=N0)
                    D(j)=NaN;
                end
            end

            % sort eigen states
            [~,Index]=sort(D);
            V=V(:,Index);
            
            % normalize k and get the required state
            if length(k)>N0+1
                k=k(1:N0+1);
            end
            k=k/norm(k);
            state=zeros(Dim,1);
            for i=1:length(k)
                state=state+k(i)*V(:,i);
            end
        otherwise
            error(['Invalid stateName: ''',stateName,'''']);
    end
    
    % reshape state
    if (strcmp(type,'psi'))
    elseif (strcmp(type,'rho'))
        state=sparse(state*state');
    else
        error(['Invalid type: ''',type,'''']);
    end
end

