function summer = plotFockState( psi,N,nn2k,ha )
%Use built-in function 'bar' to visualize psi.
%   psi - State vector
%   N   - To indicate which sub-space to visualize. Only plot basis that
%         contain exactly N bosons.
%   nn2k - Just pass in this array.
%   ha   - [Optional] Handle of axis to plot on.
% Return: The total probability that psi has exactly N bosons.

    if (nargin<4)
        ha=gca;
    end
    Nindex=zeros(N+1,1);
    for i=0:N            % find the index of N-boson state
        Nindex(i+1)=nn2k(i+1,N-i+1);
    end
    bar(ha,0:N,abs(psi(Nindex)).^2);
    set(ha,'xlim',[0 N]);
    xlabel('n_1');
    ylabel(['|c(n_1,',num2str(N),'-n_1)|^2']);
    summer=sum(abs(psi(Nindex)).^2);
end

