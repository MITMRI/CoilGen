%Real wire paths of the RF solenoid candidates from rf_solenoid_sweep.m: helix
%of rings + jogs, return strip, capacitor breaks; field, L, R at the Larmor
%frequency and capacitor values. Writes designs/rf_12cm/candidates.txt/.mat.

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
out=fullfile(pwd,'designs','rf_12cm');
S=load(fullfile(out,'sweep.mat')); spec=rf_12cm_spec();
Ns=[12 14 16 18 20 22 28];
L=0.160;

res={};
fid=fopen(fullfile(out,'candidates.txt'),'w');
hdr=sprintf('%3s %6s %6s %6s %7s %7s %7s %6s %8s %8s %8s %8s | %s\n','N','len','litz','L_uH','Rdc','Rac','Q','uT/sW', ...
    'C_brk_pF','C_ea_pF','Vbrk/A','chunk%','B1 spread % over sphere d = 80/100/110/116/120 mm, vol% of 120 mm within 5/10 %');
fprintf(fid,'RF solenoid candidates, f0 = %.4f MHz, litz %dx48 AWG, %.2f mm OD, turn radius %.2f mm\n', ...
    spec.f0/1e6,spec.litz_strands,spec.litz_od*1e3,spec.a*1e3);
fprintf(fid,'len = coil length over the litz; litz in m (incl. leads); R in ohm at 20 C; uT/sW = peak B1 per sqrt(W) of coil loss (unloaded)\n');
fprintf(fid,'C_brk = capacitance per break (%d in series), C_ea = each of the two 1206 in parallel; Vbrk/A = V peak per A peak\n',spec.num_breaks);
fprintf(fid,'chunk%% = largest deviation of a chunk length from equal\n\n%s',hdr);
fprintf('%s',hdr);
for N=Ns
i=find(cellfun(@(d) d.N==N && abs(d.L-L)<1e-9,S.designs),1);
z=S.designs{i}.z;
g=build_rf_solenoid_path(z,spec,'return_type','strip','seam_deg',90,'pocket_halfwidth_mm',0.85);
r=evaluate_rf_solenoid(g,spec);
r.N=N; r.z=z; r.g=g; r.coil_length=max(z)-min(z)+spec.litz_od;
res{end+1}=r; %#ok<SAGROW>
ch=max(abs(g.chunk_length/mean(g.chunk_length)-1))*100;
line=sprintf('%3d %6.1f %6.2f %6.2f %7.3f %7.3f %7.0f %6.1f %8.0f %8.0f %8.1f %8.1f | %s\n',N,r.coil_length*1e3, ...
    g.litz_length,r.L_H*1e6,r.R_litz_dc,r.R_ac,r.Q_unloaded,r.B1_uT_per_sqrtW,r.C_break_F*1e12,r.C_each_F*1e12, ...
    r.V_break_per_A,ch,[sprintf('%6.2f',r.spread*100) sprintf('  | %5.1f %5.1f',r.vol_within_5*100,r.vol_within_10*100)]);
fprintf(fid,'%s',line); fprintf('%s',line);
end
fclose(fid);
save(fullfile(out,'candidates.mat'),'res','spec');
