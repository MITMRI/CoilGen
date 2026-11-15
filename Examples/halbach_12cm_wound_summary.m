%Summary of the wire layouts of halbach_12cm_wound_set.m: efficiency per watt,
%linearity of the real wire on spheres of several diameters, and the places
%where neighbouring wires are closer than one wire diameter. Written to
%layout_summary.txt in the layout folder (WOUND_CLEARANCE_MM picks the folder).

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
clearance_mm=str2double(getenv('WOUND_CLEARANCE_MM')); if isnan(clearance_mm), clearance_mm=1; end
if clearance_mm==1, sub='layouts'; else, sub=['layouts_' strrep(sprintf('%gmm',clearance_mm),'.','p')]; end
folder=fullfile(pwd,'Wound_gradient_coils','halbach_12cm',sub);
files=dir(fullfile(folder,'*.mat'));
diameters=[60 80 100 110 120];
n=3000; k=(0:n-1)'+0.5; po=acos(1-2*k/n); az=pi*(1+sqrt(5))*k;
u=[cos(az).*sin(po) sin(az).*sin(po) cos(po)];

fid=fopen(fullfile(folder,'layout_summary.txt'),'w');
hdr=sprintf('%-22s %4s %5s %6s %6s %6s %5s %7s','layout','lev','wire','eff','R20','R80','L','G/sqrtW');
hdr=[hdr sprintf('  max nonlinearity [%%] on sphere of diameter [mm]:') sprintf(' %5d',diameters) sprintf('  %6s %6s','min cc','<OD')];
units=sprintf('%-22s %4s %5s %6s %6s %6s %5s %7s','','','[m]','mT/m/A','[ohm]','[ohm]','[uH]','20C');
for out=[1 fid], fprintf(out,'%s\n%s\n',hdr,units); end
for f=1:numel(files)
r=load(fullfile(folder,files(f).name)); r=r.res;
wp=r.wire_path; A=wp(:,1:end-1)'; B=wp(:,2:end)';
nl=zeros(size(diameters));
for d=1:numel(diameters)
P=u*diameters(d)/2e3; b=biot_savart(A,B,ones(size(A,1),1),P);
X=[P*r.spec.gradient_dir(:) ones(n,1)]; c=X\b(:,3);
nl(d)=max(abs(b(:,3)-X*c))/abs(c(1))/(diameters(d)/2e3)*100;
end
od=r.wire.outer_diameter;
cl=r.spacing.closest_mm; too_close=sum(cl<od*1e3)*0.8;   %mm of wire closer than OD
eff=r.field.efficiency_mT_per_m_per_A;
line=sprintf('%-22s %4d %5.1f %6.3f %6.3f %6.3f %5.0f %7.2f',erase(files(f).name,'.mat'),r.spec.levels,r.wire_length, ...
    eff,r.R20_ohm,r.R80_ohm,r.inductance_H*1e6,eff/sqrt(r.R20_ohm));
line=[line sprintf('%51s','') sprintf(' %5.2f',nl) sprintf('  %6.2f %6.0f',r.spacing.min_mm,too_close)];
for out=[1 fid], fprintf(out,'%s\n',line); end
%where the wires collide
bad=find(cl<od*1e3);
if ~isempty(bad)
p=r.spacing.points(bad,:);
ang=atan2d(p(:,3),p(:,2));
fprintf('    collisions: x %s mm, angle %s deg\n',mat2str(unique(round(p(:,1)*1e3/5)*5))',mat2str(unique(round(ang/10)*10)'));
end
end
fprintf(fid,'\nmin cc = closest centre-to-centre distance of wire parts >30 mm apart along the wire [mm];\n<OD = length of wire closer than one wire diameter to another part [mm].\n');
fclose(fid);
