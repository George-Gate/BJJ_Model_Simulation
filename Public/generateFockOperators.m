%% Description
% Generate following operators in fock representation: 
%     a1, a2
% Generate two mapping:
%     k2nn, nn2k
% Calc the dimension of Hilbert space: 
%     Dim
% Required input: (Should be set in workspace before calling this script)
%     N - maximum number of bosons

Dim=(N+1)*(N+2)/2;
a1=spalloc(Dim,Dim,Dim);
a2=spalloc(Dim,Dim,Dim);
k2nn=zeros(2,Dim);
nn2k=zeros(N+1,N+1);
% make a1 and a2
for n1=0:N
    for n2=N-n1:-1:0
        i=n2*N+(3*n2-n2^2)/2+n1+1;
        j=(n2-1)*N+(5*n2-n2^2-4)/2+n1+1;
        k=n2*N+(3*n2-n2^2)/2+(n1-1)+1;
        if n1~=0 
            a1(k,i)=sqrt(n1); %#ok<SPRIX>
        end
        if n2~=0
            a2(j,i)=sqrt(n2); %#ok<SPRIX>
        end
        % record mapping & inverse mapping
        nn2k(n1+1,n2+1)=i;
        k2nn(1,i)=n1;
        k2nn(2,i)=n2;
    end
end