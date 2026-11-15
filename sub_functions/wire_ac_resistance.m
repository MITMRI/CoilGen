function out=wire_ac_resistance(path_xyz,wire,freqs,temperature_C)
%Resistance [Ohm] of a wound coil against frequency, with skin effect and
%proximity losses from the field of all other turns.
%
%path_xyz: wire path (3xN, m). freqs: frequencies (Hz). temperature_C: copper
%temperature. wire: struct with
%   .type        'solid' or 'litz'
%   .radius      conductor radius (solid) or strand radius (litz), m
%   .num_strands litz only
%   .bundle_radius litz only, m
%   .lay_factor  litz only, added length of the twisted strands (e.g. 1.03)
%
%Solid wire: exact Bessel solutions for a round conductor, one for its own
%current (skin) and one in a uniform transverse field (proximity), with the
%local field taken from Biot-Savart along the path. Litz: the same proximity
%solution per strand, in the external field plus the bundle's own field.
%Returns .R_dc, .R_ac (per frequency), .H_rms_per_A (A/m per A along the wire).

rho=1.724e-8*(1+0.00393*(temperature_C-20)); %annealed copper, 100% IACS
mu0=4*pi*1e-7;
sigma=1/rho;

step=vecnorm(diff(path_xyz,1,2));
wire_length=sum(step);

%local transverse field at sample points along the wire, for 1 A
A=path_xyz(:,1:end-1)'; B=path_xyz(:,2:end)';
num_samples=min(12000,size(A,1));
sample_inds=round(linspace(1,size(A,1),num_samples));
mid=(A(sample_inds,:)+B(sample_inds,:))/2;
tangent=(B(sample_inds,:)-A(sample_inds,:)); tangent=tangent./vecnorm(tangent,2,2);
H=zeros(num_samples,3);
block=400;
for i0=1:block:num_samples
ii=i0:min(i0+block-1,num_samples);
H(ii,:)=field_without_own_segments(A,B,mid(ii,:),sample_inds(ii))/mu0;
end
H_transverse2=sum(H.^2,2)-sum(H.*tangent,2).^2;
out.H_rms_per_A=sqrt(mean(H_transverse2));
mean_H2=mean(H_transverse2);

switch wire.type
case 'solid'
area=pi*wire.radius^2;
out.R_dc=rho*wire_length/area;
skin=arrayfun(@(f) skin_factor(wire.radius,sigma,f),freqs);
prox=arrayfun(@(f) proximity_loss(wire.radius,sigma,f),freqs);
out.R_ac=out.R_dc*skin+2*wire_length*prox*mean_H2;
case 'litz'
area=wire.num_strands*pi*wire.radius^2;
out.R_dc=rho*wire_length*wire.lay_factor/area;
H_bundle2=1/(8*pi^2*wire.bundle_radius^2); %mean H^2 inside the bundle for 1 A
skin=arrayfun(@(f) skin_factor(wire.radius,sigma,f),freqs);
prox=arrayfun(@(f) proximity_loss(wire.radius,sigma,f),freqs);
out.R_ac=out.R_dc*skin+2*wire_length*wire.lay_factor*wire.num_strands*prox*(mean_H2+H_bundle2);
end
out.freqs=freqs;
end

function Fr=skin_factor(a,sigma,f)
%R_ac/R_dc of a round wire from its own current
if f==0, Fr=1; return; end
delta=sqrt(1/(pi*f*4e-7*pi*sigma));
x=(1-1i)*a/delta;
Fr=real(x/2*besselj(0,x,1)/besselj(1,x,1));
end

function g=proximity_loss(a,sigma,f)
%mean loss per length [W/m] of a round wire in a transverse field of 1 A/m peak
if f==0, g=0; return; end
mu0=4*pi*1e-7; w=2*pi*f;
delta=sqrt(2/(w*mu0*sigma));
x=(1-1i)*a/delta;
Q=x*besselj(0,x,1)/besselj(1,x,1)-1;
c=(1-Q)/(1+Q);
g=pi/2*w*mu0*a^2*abs(real(1i*(1+c)*conj(1-c)));
end

function b=field_without_own_segments(A,B,P,own)
%Biot-Savart of the straight segments A->B at points P for 1 A, leaving out
%the segments next to each point along the wire
mu0_over_4pi=1e-7;
L=B-A; len=vecnorm(L,2,2); Lh=L./len;
b=zeros(size(P));
for k=1:size(P,1)
r1=P(k,:)-A; r2=P(k,:)-B;
n1=vecnorm(r1,2,2); n2=vecnorm(r2,2,2);
c=cross(L,r1,2); c2=sum(c.^2,2); c2(c2<1e-24)=inf;
f=mu0_over_4pi*(sum(Lh.*r1,2)./n1-sum(Lh.*r2,2)./n2).*len./c2;
near=max(1,own(k)-2):min(size(A,1),own(k)+2);
f(near)=0;
b(k,:)=sum(f.*c,1);
end
end
