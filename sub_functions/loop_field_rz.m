function [Br,Bz]=loop_field_rz(a,z_loops,rho,z)
%LOOP_FIELD_RZ  Field [T] of coaxial circular loops (radius a, at axial
%positions z_loops, 1 A each, same direction, axis = z) at points (rho,z)
%(column vectors, m). Exact (complete elliptic integrals).
mu0=4*pi*1e-7;
rho=rho(:); z=z(:); z_loops=z_loops(:)';
dz=z-z_loops;                                   %M x N
s=a^2+rho.^2+dz.^2;
alpha2=s-2*a*rho; beta2=s+2*a*rho; beta=sqrt(beta2);
m=1-alpha2./beta2;
[K,E]=ellipke(min(m,1-1e-15));
C=mu0/pi;
Bz=sum(C./(2*alpha2.*beta).*((a^2-rho.^2-dz.^2).*E+alpha2.*K),2);
rr=max(rho,1e-12);
Br=sum(C*dz./(2*alpha2.*beta.*rr).*(s.*E-alpha2.*K),2);
Br(rho<1e-9)=0;
end
