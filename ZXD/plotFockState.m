%Author: George-Gate
%Date: 2015/10/6
%Last Modify Date: 2015/10/10
%
function P = plotFockState( varargin )
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
%  P=plotFockState( rho,N0,proj,ha)
%   rho - Density matrix size: Dim x Dim
%   N0  - To indicate which sub-space to visualize. Only plot basis that
%         contain exactly N0 bosons.
%   proj - Just pass in this cell matrix.
%   ha   - [Optional] Handle of axis to plot on.
% Return: The probability distribution that being plot.

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
        set(ha,'xlim',[0 N0]);
        xlabel('n_1');
        ylabel(['|c(n_1,',num2str(N0),'-n_1)|^2']);
        P=abs(psi(Nindex)).^2;
    elseif (size(varargin{1},1)==size(varargin{1},2))
    % if the first argument is a density matrix
        rho=varargin{1};
        N0=varargin{2};
        proj=varargin{3};
        proTmp=zeros(N0+1,1);
        % calc projection
        for i=0:N0
            proTmp(i+1)=abs(trace(proj{i+1,N0-i+1}*rho));
        end
        % plot
        bar(ha,0:N0,proTmp);
        set(ha,'xlim',[0 N0]);
        xlabel('n_1');
        ylabel(['|c(n_1,',num2str(N0),'-n_1)|^2']);
        P=proTmp;
    else
        % The shape of the first argument is invalid.
        error('Syntax error: The shape of the first argument is invalid.');
    end
end

