%Author: George-Gate
%Date: 2015/10/6
%Last Modify Date: 2015/10/10
%--------------------------------------------------------------------------
%Use built-in function 'bar' to visualize psi.
%
%[Usage 1]
%  P=plotFockState( psi,N0,nn2k,ha )
%   psi - State vector size: Dim x 1
%   N0  - To indicate which sub-space to visualize. Only plot basis that
%         contain exactly N0 bosons.
%   nn2k - Just pass in this array.
%   ha   - [Optional] Handle of axis to plot on.
% Return: The probability distribution that being plot.
%[Usage 2]
%  P=plotFockState( rho,N0_Range,nn2k,ha)
%   rho - Density matrix size: Dim x Dim
%   N0_Range - A vector contains two integers. Indicates what sub-space to visualize. Only plot basis that
%         contain N0_Range(1)~N0_Range(2) (included) bosons.
%   nn2k - Just pass in this array.
%   ha   - [Optional] Handle of axis to plot on.
% Return: The probability distribution that being plot.
function P = plotFockState( varargin )
    % set axis handle
    if (nargin<4)
        ha=gca;
    else
        ha=varargin{4};
    end
    if (size(varargin{1},2)==1)
    % if the first argument is a state vector
        psi=varargin{1};
        N0=varargin{2};
        nn2k=varargin{3};
        Nindex=zeros(N0+1,1);
        for i=0:N0            % find the index of N-boson state
            Nindex(i+1)=nn2k(i+1,N0-i+1);
        end
        bar(ha,0:N0,abs(psi(Nindex)).^2);
        set(ha,'xlim',[-0.5 N0+0.5]);
        xlabel('n_1');
        ylabel(['|c(n_1,',num2str(N0),'-n_1)|^2']);
        P=abs(psi(Nindex)).^2;
    elseif (size(varargin{1},1)==size(varargin{1},2))
    % if the first argument is a density matrix
        rho=varargin{1};
        N0_Range=varargin{2};
        if length(N0_Range)==1
            N0_Range(2)=N0_Range(1);
        end
        if (N0_Range(2)<N0_Range(1))
            error('Invalid N0_Range.');
        end
        nn2k=varargin{3};
        proTmp=zeros(N0_Range(2)+1,N0_Range(2)-N0_Range(1)+1);
        % calc projection
        for n=N0_Range(1):N0_Range(2)
            for i=0:n
                proTmp(i+1,n-N0_Range(1)+1)=abs(rho(nn2k(i+1,n-i+1),nn2k(i+1,n-i+1))); 
                %proTmp(i+1,n-N0_Range(1)+1)=abs(trace(proj{i+1,n-i+1}*rho));
            end
        end
        % plot
        if (N0_Range(1)==N0_Range(2))
            bar(ha,0:N0_Range(1),proTmp);
            set(ha,'xlim',[-0.5 N0_Range(1)+0.5]);
            xlabel('n_1');
            ylabel(['|c(n_1,',num2str(N0_Range(1)),'-n_1)|^2']);
        else
            bar3(ha,0:N0_Range(2),proTmp);
            set(gca,'ylim',[-0.5,N0_Range(2)+0.5]);
            % generate tick labels
            tl=cell(N0_Range(2)-N0_Range(1)+1,1);
            for i=N0_Range(1):N0_Range(2)
                tl{i-N0_Range(1)+1}=num2str(i);
            end
            set(gca,'XTickLabel',tl);
            xlabel('n_1+n_2');
            ylabel('n_1');
            zlabel('Probability');
        end
        P=proTmp;
    else
        % The shape of the first argument is invalid.
        error('Syntax error: The shape of the first argument is invalid.');
    end
end

