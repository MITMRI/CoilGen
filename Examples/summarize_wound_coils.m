%Summary of the wire gauge study of halbach_wound_gradient_set.m, printed and
%written to coil_summary.txt next to the results.

if ~exist('skip_cd','var')
if ispc
cd('..\');
else
cd('../');
end
end
result_folder=fullfile(pwd,'Wound_gradient_coils','halbach_320mm_bore');
coil_names={'Z','Y','X'}; targets=[28 12 12]; %mT/m
wire_names={'AWG12','AWG14','AWG16','AWG18','LITZ16'};
ramps=[500e-6 1e-3];

fid=fopen(fullfile(result_folder,'coil_summary.txt'),'w');
for out_id=[1 fid]
fprintf(out_id,'Current, voltage and power for the target gradient, copper at 80 C\n');
fprintf(out_id,'%-3s %-7s %4s %6s %6s %6s %6s %6s %7s %7s %6s %6s %7s %7s %7s\n','','wire','turn','wire','pitch','clos.','eff','nonl','R20','R80','L','I','V_R','V_L','P');
fprintf(out_id,'%-3s %-7s %4s %6s %6s %6s %6s %6s %7s %7s %6s %6s %7s %7s %7s\n','','','lev','[m]','[mm]','[mm]','mT/m/A','[%]','[ohm]','[ohm]','[mH]','[A]','[V]','1ms[V]','[W]');
for c=1:3
for w=1:numel(wire_names)
file=fullfile(result_folder,'gauge_study',[coil_names{c} '_' wire_names{w} '.mat']);
if ~exist(file,'file'), continue; end
r=load(file); r=r.res; f=r.field;
I=targets(c)/f.efficiency_mT_per_m_per_A;
fprintf(out_id,'%-3s %-7s %4d %6.1f %6.2f %6.2f %6.3f %6.2f %7.3f %7.3f %6.3f %6.2f %7.1f %7.1f %7.0f\n',coil_names{c},wire_names{w},r.spec.levels, ...
    r.wire_length,r.wire.outer_diameter*1e3+1,r.spacing.min_mm,f.efficiency_mT_per_m_per_A,f.nonlinearity_max_percent, ...
    r.R20_ohm,r.R80_ohm,r.inductance_H*1e3,I,I*r.R80_ohm,r.inductance_H*I/ramps(2),I^2*r.R80_ohm);
end
fprintf(out_id,'\n');
end
fprintf(out_id,'AC resistance, solid wire model at 80 C: R(f)/R_dc\n');
fprintf(out_id,'%-3s %-7s %8s %8s %8s %8s\n','','wire','100 Hz','1 kHz','3 kHz','10 kHz');
for c=1:3
for w=1:numel(wire_names)
file=fullfile(result_folder,'gauge_study',[coil_names{c} '_' wire_names{w} '.mat']);
if ~exist(file,'file'), continue; end
r=load(file); r=r.res;
fprintf(out_id,'%-3s %-7s %8.4f %8.4f %8.4f %8.4f\n',coil_names{c},wire_names{w},r.ac.R_ac/r.ac.R_dc);
end
end
end
fclose(fid);
