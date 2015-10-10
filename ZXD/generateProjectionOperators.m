%% Description
% Generate following operators in fock representation: 
%     proj{n1+1,n2+1}: make the projection to state |n1,n2> (n1+n2<=N)
%     projNOON{n}: make the projection to state (|n,0>+|0,n>)/sqrt(2) (n<=N)
% Required input: (Should be set in workspace before calling this script)
%     (Run 'generateFockOperators.m' first.) 
%     N, Dim, nn2k
%
% Encoding: See 'generateFockOperators.m'
%

proj=cell(N+1,N+1);
projNOON=cell(N,1);

% make proj
for n1=0:N
    for n2=N-n1:-1:0
        i=nn2k(n1+1,n2+1);
        proj{n1+1,n2+1}=spalloc(Dim,Dim,1);
        proj{n1+1,n2+1}(i,i)=1;
    end
end

% make projNOON
for n=1:N
    i=nn2k(n+1,1);
    j=nn2k(1,n+1);
    projNOON{n}=spalloc(Dim,Dim,4);
    projNOON{n}(i,i)=0.5;
    projNOON{n}(i,j)=0.5;
    projNOON{n}(j,i)=0.5;
    projNOON{n}(j,j)=0.5;
end
clear n n1 n2 i j;