% filename & parameters
path='../mats/RUN2.6/';
iRange=[1,100];

temp=load([path,'i=',num2str(iRange(1)),'.mat']);

% init storage space
len=iRange(2)-iRange(1)+1;
tmaxList  = zeros(len,1);
c0xList   = zeros(len,1);
cx0List   = zeros(len,1);
cx1List   = zeros(len,1);
finPsiList= zeros(temp.results.Dim,len);
rCount=0;

% load data
for i=iRange(1)+1:iRange(2)+1
    rCount=rCount+1;
    tmaxList(rCount)=temp.results.tmax;
    c0xList(rCount)=temp.results.c0x;
    cx0List(rCount)=temp.results.cx0;
    cx1List(rCount)=temp.results.cx1;
    finPsiList(:,rCount)=temp.results.psiList(:,temp.results.rCount);
    
    if (i<=iRange(2))
        temp=load([path,'i=',num2str(i),'.mat']);
    end  
end

% save data
save([path,'combinedData.mat'],'rCount',...
     'tmaxList','c0xList','cx0List','cx1List','finPsiList');
 
 clear temp path iRange len tmaxList c0xList cx0List cx1List;
 clear finPsiList rCount;

