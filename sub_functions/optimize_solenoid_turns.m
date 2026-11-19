function out=optimize_solenoid_turns(N,a,z_max,p_min,roi_radius)
%OPTIMIZE_SOLENOID_TURNS  Axial positions of N circular turns (radius a, m) of
%a solenoid along z, symmetric about z = 0, that make the RF field B1 most
%uniform over a sphere of radius roi_radius. B0 is along x, so B1 is the part
%of the field perpendicular to x: |B_perp|^2 = Bz^2 + (Br*sin(phi))^2. Over
%phi it lies between |Bz| (phi = 0) and sqrt(Bz^2+Br^2) (phi = 90 deg), so
%both are constrained on a (rho,z) grid.
%
%Minimises the peak-to-peak spread (max-min)/(max+min) of |B_perp| over the
%sphere, with turn centres at |z| <= z_max and at least p_min apart (m).
%Returns .z (1xN, m), .spread, .B_centre_per_A (T/A).

[pr,pz]=quarter_disc_points(roi_radius,0.002);
odd=mod(N,2)==1;
nh=floor(N/2);
if nh==0, error('N must be at least 2'); end
first_min=p_min/2; if odd, first_min=p_min; end
if z_max<first_min+(nh-1)*p_min, error('%d turns do not fit',N); end

%linear constraints on h (positive half): h1>=first_min, h(i+1)-h(i)>=p_min, h(end)<=z_max
Aineq=zeros(nh+1,nh+2); bineq=zeros(nh+1,1);
Aineq(1,1)=-1; bineq(1)=-first_min;
for i=1:nh-1, Aineq(i+1,i)=1; Aineq(i+1,i+1)=-1; bineq(i+1)=-p_min; end
Aineq(nh+1,nh)=1; bineq(nh+1)=z_max;
full_z=@(h) [-fliplr(h(:)') zeros(1,odd) h(:)'];

scale=4*pi*1e-7*N/(2*a);                         %rough field, keeps values ~1
field=@(h) b_perp_bounds(a,full_z(h),pr,pz)/scale;

starts={};
starts{end+1}=linspace(first_min,z_max,nh);                     %uniform
u=linspace(0,1,nh); starts{end+1}=first_min+(z_max-first_min)*sin(u*pi/2); %end weighted
starts{end+1}=first_min+(z_max-first_min)*u.^0.7;
opts=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxFunctionEvaluations',4e4,'MaxIterations',800,'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-8);
best=inf;
for s=1:numel(starts)
h0=starts{s};
if any(diff(h0)<p_min), h0=first_min+(0:nh-1)*(z_max-first_min)/max(nh-1,1); end
b=field(h0);
x0=[h0(:);min(b);max(b)];
obj=@(x) (x(end)-x(end-1))/(x(end)+x(end-1));
nonl=@(x) deal([field(x(1:nh))-x(end); x(end-1)-field(x(1:nh))],[]);
lb=[zeros(nh,1);0;0]; ub=[z_max*ones(nh,1);inf;inf];
[x,f]=fmincon(obj,x0,Aineq,bineq,[],[],lb,ub,nonl,opts);
if f<best, best=f; xb=x; end
end
h=sort(xb(1:nh))';
out.z=full_z(h);
out.spread=best;
[~,Bz0]=loop_field_rz(a,out.z,0,0);
out.B_centre_per_A=abs(Bz0);
end

function b=b_perp_bounds(a,z,pr,pz)
[Br,Bz]=loop_field_rz(a,z,pr,pz);
b=[abs(Bz);sqrt(Bz.^2+Br.^2)];
end

function [pr,pz]=quarter_disc_points(R,h)
%(rho,z) points in the quarter disc rho,z >= 0 of radius R, plus its rim
[g1,g2]=ndgrid(0:h:R);
keep=g1.^2+g2.^2<=R^2;
t=linspace(0,pi/2,ceil(R*pi/2/h)+1)';
pr=[g1(keep);R*cos(t)]; pz=[g2(keep);R*sin(t)];
end
