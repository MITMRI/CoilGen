function r=evaluate_rf_solenoid(g,spec)
%EVALUATE_RF_SOLENOID  B1 field, inductance, resistance at the Larmor frequency
%and the capacitor values of an RF solenoid built by build_rf_solenoid_path.
w=2*pi*spec.f0;
r.f0=spec.f0;

%--- field at 1 A: litz + strip (5 filaments of 1/5 A) ---
A=g.litz(:,1:end-1)'; B=g.litz(:,2:end)'; I=ones(size(A,1),1);
for i=1:numel(g.strip_filaments)
f=g.strip_filaments{i};
%strip runs from bottom to top (current flows up it back to the feed)
A=[A; f(:,1:end-1)']; B=[B; f(:,2:end)']; I=[I; ones(size(f,2)-1,1)/numel(g.strip_filaments)]; %#ok<AGROW>
end
B0=biot_savart(A,B,I,[0 0 0]);
r.B1_centre_T_per_A=hypot(B0(2),B0(3));
d=[0.060 0.080 0.100 0.110 0.120];
r.sphere_d=d; r.spread=zeros(size(d)); r.dev_min=r.spread; r.dev_max=r.spread;
for i=1:numel(d)
P=sphere_points(d(i)/2);
Bp=biot_savart(A,B,I,P);
b=hypot(Bp(:,2),Bp(:,3));                       %perpendicular to B0 (x)
r.spread(i)=(max(b)-min(b))/(max(b)+min(b));
r.dev_min(i)=min(b)/r.B1_centre_T_per_A-1; r.dev_max(i)=max(b)/r.B1_centre_T_per_A-1;
end

%--- inductance of the closed loop: litz as a round 2.5 mm conductor, strip as
%its centreline with the same radius (a small part of the total) ---
%closed path: litz (top -> bottom), then the strip (bottom -> top), closed through the feed
loop=[g.litz g.strip_centre];
r.L_H=wire_inductance(loop,spec.litz_od/2,'closed',true,'segment_length',0.003);

%--- resistance at f0, 20 C and 40 C ---
litz_model=struct('type','litz','radius',spec.litz_strand_radius,'num_strands',spec.litz_strands, ...
    'bundle_radius',spec.litz_od/2*0.9,'lay_factor',spec.litz_lay_factor);
ac=wire_ac_resistance(g.litz,litz_model,spec.f0,20);
r.R_litz_dc=ac.R_dc; r.R_litz_ac=ac.R_ac;
rho=1.724e-8; mu0=4*pi*1e-7;
delta=sqrt(rho/(pi*spec.f0*mu0));
wst=g.options.strip_width_mm*1e-3;
r.R_strip=rho*g.strip_length/(2*wst*delta)*1.3;  %both faces, +30 % for edge crowding
if isfield(g,'return_type') && strcmp(g.return_type,'bore')
r.R_strip=r.R_litz_ac/g.litz_length*g.strip_length;   %litz return, same loss per metre
end
r.R_ac=r.R_litz_ac+r.R_strip;
r.Q_unloaded=w*r.L_H/r.R_ac;
r.B1_uT_per_sqrtW=r.B1_centre_T_per_A/sqrt(r.R_ac)*1e6;  %peak B1 for 1 W, coil loss only (I peak = sqrt(2 W/R): x sqrt 2)
r.B1_uT_per_sqrtW=r.B1_uT_per_sqrtW*sqrt(2);

%--- tuning: num_breaks equal capacitors in series ---
C_total=1/(w^2*r.L_H);
r.C_total_F=C_total;
r.C_break_F=spec.num_breaks*C_total;
r.C_each_F=r.C_break_F/2;                         %two 1206 in parallel per break
r.V_break_per_A=1/(w*r.C_break_F);               %V peak across a break per A peak
end

function P=sphere_points(R)
n=3000; k=(0:n-1)'+0.5;
pol=acos(1-2*k/n); az=pi*(1+sqrt(5))*k;
S=R*[cos(az).*sin(pol) sin(az).*sin(pol) cos(pol)];
[gx,gy,gz]=ndgrid(linspace(-1,1,15)*R);
In=[gx(:) gy(:) gz(:)]; In=In(vecnorm(In,2,2)<0.999*R,:);
P=[S;In];
end
