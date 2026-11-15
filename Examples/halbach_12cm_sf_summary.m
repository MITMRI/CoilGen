%Best length for each linearity limit, from the stream function sweep
%(halbach_12cm_sf_sweep.m), estimated for one layer of round wire.
%Wire count is set by the densest place: neighbouring wires are one pitch
%apart where the current density is highest.
%   efficiency = 1/dpsi,  wire length = int|K| dA / dpsi,  dpsi = Kmax*pitch
%   gradient per sqrt(watt) = 1/sqrt(ohm_per_m * dpsi * int|K| dA)

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
folder=fullfile(pwd,'Wound_gradient_coils','halbach_12cm','sf_sweep');
inch=0.0254; per_kft=1/304.8; clearance=0.001;
wires=struct('name',{'AWG14','AWG16'},'od',{0.0682*inch,0.0545*inch}, ...
    'ohm_per_m',{2.524*per_kft,4.019*per_kft});
names={'Z','Y','X'}; lengths=140:20:340;
limits=[0.05 0.1 0.2 0.5 1 2 5];   %max nonlinearity on the 120 mm sphere, %

fid=fopen(fullfile(folder,'..','sf_summary.txt'),'w');
for out=[1 fid]
fprintf(out,'Best gradient per sqrt(W) at 20 C, 1 mm between wires, max nonlinearity on the 120 mm sphere\n');
fprintf(out,'(stream function estimate; mT/m per sqrt(W); length in mm; turns = contour levels)\n\n');
for w=1:numel(wires)
pitch=wires(w).od+clearance;
fprintf(out,'%s  pitch %.2f mm\n',wires(w).name,pitch*1e3);
fprintf(out,'%-3s','');
fprintf(out,'  %7s%%: G/sqrtW  L  lev  mT/m/A',num2str(limits(1)));
for l=limits(2:end), fprintf(out,' | %4s%%: G/sqW   L lev  eff',num2str(l)); end
fprintf(out,'\n');
for c=1:3
fprintf(out,'%-3s',names{c});
for lim=limits
best=struct('fom',0);
for L=lengths
r=load(fullfile(folder,sprintf('%s_L%03d.mat',names{c},L))); r=r.res;
ok=r.nonlin_max_pct<=lim;
dpsi=r.K_max*pitch;
levels=floor(r.sf_range./dpsi);
dpsi=r.sf_range./max(levels,1);
fom=1./(wires(w).ohm_per_m*dpsi.*r.K_int);  %(T/m)^2/W
fom(~ok)=0;
[m,k]=max(fom);
if m>best.fom, best=struct('fom',m,'L',L,'k',k,'lev',levels(k),'eff',1e3/dpsi(k)); end
end
if best.fom>0
fprintf(out,' | %6.2f %3d %3d %5.2f',sqrt(best.fom)*1e3,best.L,best.lev,best.eff);
else
fprintf(out,' | %6s %3s %3s %5s','-','','','');
end
end
fprintf(out,'\n');
end
fprintf(out,'\n');
end
end
fclose(fid);
