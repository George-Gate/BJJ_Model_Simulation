%Author: George-Gate
%Date: 2016/11/16
%Last Modify Date: 2016/11/17
%--------------------------------------------------------------------------
%plot the Husimi Q function of a multi-boson spin state
%
%[Usage 1]
%  P=plotSpinState( psi,N,M,SCSbase,ha )
%   psi - State vector size: (2*M)^N x 1
%   N   - The number of bosons
%   M   - The dim of spatial state
% SCSbase - One of the output of this function, use [] if do not want this
%           input.
%   ha  - [Optional] Handle of axis to plot on.
% Return: 
%   Q      - The husimi Q function that being plot with coordinate
%           [theta,phi]=meshgrid(linspace(0,pi,50),linspace(0,2*pi,50));
% SCSbase  - The spin coherent state base that used to calc Q function. 
%            A (2*M)^N x 50 x 50 array whose first dim is the same as psi.
%            Can be inputed to this function at the next call for speed up.
%


function [Q,SCSbase] = plotSpinState( varargin )
    % set axis handle
    if (nargin<5)
        ha=gca;
    else
        ha=varargin{5};
    end
    if (size(varargin{1},2)==1)
    % if the first argument is a state vector
        % read args 
        psi=varargin{1};
        N=varargin{2};
        M=varargin{3};
        SCSbase=varargin{4};
        
        % check SCSbase
        gridSize=[50 50];
        if (  ~(size(SCSbase,1)==(2*M)^N && size(SCSbase,2)==gridSize(1) && size(SCSbase,3)==gridSize(2)  )     )
            % generate SCSbase
            SCSbase=zeros([(2*M)^N,gridSize]);
            
            psi_x0=zeros(M,1);psi_x0(1)=1;   % spatial state for single boson
            
            % generate grid
            [theta,phi]=meshgrid(linspace(0,pi,gridSize(1)),linspace(0,2*pi,gridSize(2)));
            % calc coherent state
            for i=1:size(phi,1)
                for j=1:size(phi,2);
                    th=theta(i,j);
                    ph=phi(i,j);
                    % generate spin coherent state |th,ph>_SCS
                    psi_a0=kron([cos(th/2)*exp(-0.5i*ph);sin(th/2)*exp(0.5i*ph)],psi_x0);
                    psi_a=psi_a0;
                    for k=2:N
                        psi_a=kron(psi_a0,psi_a);
                    end
                    SCSbase(:,i,j)=psi_a;
                end
            end
        end 
        
        % generate grid
        if (~exist('theta','var'))
            [theta,phi]=meshgrid(linspace(0,pi,gridSize(1)),linspace(0,2*pi,gridSize(2)));
        end
        % calc Husimi Q function
        Q=phi;
        psi=psi';       % transpose in advanced for speed up
        for i=1:size(phi,1)
            for j=1:size(phi,2);
                Q(i,j)=abs(psi*SCSbase(:,i,j))^2/pi;
            end
        end
        psi=psi';       % transpose back
        
        X=sin(theta).*cos(phi);
        Y=sin(theta).*sin(phi);
        Z=cos(theta);
        % plot grid line
        [gth,gph]=meshgrid(linspace(0,pi,20),linspace(0,2*pi,20));
        h=mesh(ha,sin(gth).*cos(gph),sin(gth).*sin(gph),cos(gth));hold on;
        set(h,'FaceAlpha',0,'EdgeColor',[0.7 0.7 0.7]);
        % plot Q function
        h=surf(ha,X,Y,Z,Q);hold off;
        %h=surf(ha,theta/pi,phi/pi,Q);
        h.AlphaData=Q/max(Q(:));
        set(h,'EdgeAlpha',0,'FaceColor','interp','FaceAlpha','interp');
        set(ha,'CLim',[min(Q(:)),max(Q(:))]);
        
        % set view
            %set(ha,'CameraPosition',[0.5 -1 0.5],'CameraTarget',[-0.5 1 -0.5]);
            set(ha,'CameraPosition',[13 -10 6],'CameraTarget',[0 0 0]);
        
        xlabel('x (Husimi Q Function)');
        ylabel('y');
        zlabel('z');
    else
        % The shape of the first argument is invalid.
        error('Syntax error: The shape of the first argument is invalid.');
    end
end

