function psi_t=BJJ_SEQ(t,psi,tmax,relaxT,H1,H2)
    J = (1-t/(tmax-relaxT))*(t<tmax-relaxT);
    Ec=  -t/(tmax-relaxT)*(t<tmax-relaxT)-(t>=tmax-relaxT);
    H=sparse(-J*H1+Ec/8*H2);
    psi_t=-1i*H*psi;
end