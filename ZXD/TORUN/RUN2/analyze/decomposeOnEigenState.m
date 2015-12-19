J=5;
Ec=-1;
N=10;
N0=10;

generateFockOperators;
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
H=sparse(-J*H1+Ec/8*H2);
generateEigenStates;

tmp=load('..\mats\RUN2.1\combinedData.mat');
coeff=zeros(tmp.rCount,5);
for i=1:tmp.rCount
    for j=1:5
        coeff(i,j)=tmp.finPsiList(:,i)'*V(:,j);
    end
end

for i=1:tmp.rCount
    coeff(i,:)=coeff(i,:)*exp(-1i*angle(coeff(i,1)));
end

plot(linspace(0,2*pi,100)/pi,abs(coeff(:,1)).^2);hold on;
plot(linspace(0,2*pi,100)/pi,abs(coeff(:,2)).^2);
plot(linspace(0,2*pi,100)/pi,abs(coeff(:,3)).^2);
plot(linspace(0,2*pi,100)/pi,abs(coeff(:,4)).^2);
%scatter(linspace(0,2*pi,100)/pi,abs(tmp.cx0List-tmp.c0xList).^2/2);
legend({'n=0','n=1','n=2','n=3'});
ylabel('Probability');
xlabel('\phi/\pi');
title('N=10, tmax=220, Final State');
% plot(tmp.tmaxList,prob(:,2));hold on;
% plot(tmp.tmaxList,prob(:,3));hold on;
% plot(tmp.tmaxList,prob(:,4));hold on;
% plot(tmp.tmaxList,prob(:,5));hold on;
hold off;