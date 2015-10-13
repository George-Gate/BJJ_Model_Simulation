% Author: LCY
% Date: 2015/10/12
% Last edited by George-Gate 2015/10/13
%--------------------------------------------------------------------------
% This is the function that calculate the QFI of a given state psi
% based on a pseudo unitary transformation U=exp(-i\omega S_z t)
%
% QFI = QFI(state,a1,a2)
%
%  state -    State vector or density matrix
%  a1/a2 -    Annihilation operator of two modes
 

function QFI = QFI(state,a1,a2)
    H=(a2'*a2-a1'*a1)/2;  % H=S_z, the generator of U
    omega=1;
    dt=0.01;              % A short time for transformation
    if (isvector(state))
        psi_th=expm(-1i*H*omega*dt)*state;
        psi_th_d=-1i*H*psi_th;
        QFI=real(4*(psi_th_d'*psi_th_d-(psi_th_d'*psi_th)^2));
    else
        rho_th=expm(-1i*H*omega*dt)*state*expm(1i*H*omega*dt);
        rho_th_d=-1i*H*rho_th;
        QFI=4*real(trace(rho_th_d'*rho_th_d)-(trace(rho_th_d'))^2);
    end
end

