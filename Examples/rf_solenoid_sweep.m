%Sweep of the 12 cm RF solenoid: for each number of turns N and coil length,
%the turn positions that make B1 most uniform over the optimisation sphere
%(optimize_solenoid_turns, ring model), and how uniform that is over spheres
%of other sizes. Litz length per turn is fixed by the radius, so N sets how
%much litz is used. Writes designs/rf_12cm/sweep.txt and sweep.mat.
%
%Magnet: 72 mT Halbach, B0 along x, bore along z -> solenoid along z.

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
out=fullfile(pwd,'designs','rf_12cm'); if ~exist(out,'dir'), mkdir(out); end

spec=rf_12cm_spec();
Ns=4:30;
lengths=[0.100 0.120 0.140 0.160];
eval_d=[0.100 0.110 0.116 0.120];

rows=[];
designs={};
for L=lengths
z_max=L/2-spec.litz_od/2;
for N=Ns
try
o=optimize_solenoid_turns(N,spec.a,z_max,spec.pitch_min,spec.roi_radius);
catch
continue
end
sp=arrayfun(@(d) ring_spread(spec.a,o.z,d/2),eval_d);
litz=N*2*pi*spec.a;
eff=o.B_centre_per_A/sqrt(litz);                 %T/A per sqrt(m of litz): ~ B1/sqrt(W)
rows(end+1,:)=[L N litz o.B_centre_per_A*1e6 eff*1e6 sp*100]; %#ok<SAGROW>
designs{end+1}=struct('L',L,'N',N,'z',o.z); %#ok<SAGROW>
fprintf('L %3.0f  N %2d  litz %5.2f m  spread(%%) %s\n',L*1e3,N,litz,sprintf('%6.2f',sp*100));
end
end
save(fullfile(out,'sweep.mat'),'rows','designs','spec','eval_d');

fid=fopen(fullfile(out,'sweep.txt'),'w');
fprintf(fid,'RF solenoid sweep (ring model), turns at radius %.2f mm, pitch >= %.1f mm, optimised over a %.0f mm sphere\n', ...
    spec.a*1e3,spec.pitch_min*1e3,2*spec.roi_radius*1e3);
fprintf(fid,'B1 spread = (max-min)/(max+min) of |B perpendicular to B0| over the sphere, in %%\n\n');
fprintf(fid,'%6s %3s %7s %10s %13s %s\n','L_mm','N','litz_m','uT_per_A','uT/A/sqrt(m)',sprintf('  d%-5.0f',eval_d*1e3));
for i=1:size(rows,1)
fprintf(fid,'%6.0f %3d %7.2f %10.2f %13.2f %s\n',rows(i,1)*1e3,rows(i,2),rows(i,3),rows(i,4),rows(i,5),sprintf('%8.2f',rows(i,6:end)));
end
fclose(fid);

function s=ring_spread(a,z,R)
[g1,g2]=ndgrid(0:0.001:R);
keep=g1.^2+g2.^2<=R^2;
t=linspace(0,pi/2,200)';
pr=[g1(keep);R*cos(t)]; pz=[g2(keep);R*sin(t)];
[Br,Bz]=loop_field_rz(a,z,[pr;pr],[pz;-pz]);
b=[abs(Bz);sqrt(Bz.^2+Br.^2)];
s=(max(b)-min(b))/(max(b)+min(b));
end
