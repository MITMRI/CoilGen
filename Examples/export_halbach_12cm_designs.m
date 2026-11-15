%Export the chosen 12 cm Halbach gradient coils (halbach_12cm_wound_set.m,
%0.8 mm between wires) to designs/12cm for CAD: the wire centreline of each
%coil in the magnet frame as CSV, a design summary and the full results
%(run view_halbach_12cm_coils.m afterwards for the figures). Terminals that are
%not at the +z end are moved there (move_wire_terminals). See designs/12cm/README.md.
%
%Magnet frame: x = B0, z = bore, y = z cross x.  From the CoilGen frame
%(B0 = z_cg, bore = x_cg): x = z_cg, y = -y_cg, z = x_cg.

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
src=fullfile(pwd,'Wound_gradient_coils','halbach_12cm','layouts_0p8mm');
out=fullfile(pwd,'designs','12cm');
for d={'','wire_paths','results','figures'}, if ~exist(fullfile(out,d{1}),'dir'), mkdir(fullfile(out,d{1})); end, end

%       coil  layout file              former inner/outer radius [mm]
coils={'Z','Z_L330_AWG14_lin300',66,70; ...
       'Y','Y_L324_AWG14_lin100',71,75; ...
       'X','X_L318_AWG14_lin020',76,81};

fid=fopen(fullfile(out,'design_summary.txt'),'w');
fprintf(fid,'12 cm Halbach gradient set: 14 AWG heavy build magnet wire, 0.8 mm plastic between wires\n');
fprintf(fid,'Magnet frame: x = B0, z = bore. All lengths in mm. Copper at 20 C unless noted.\n\n');
for c=1:size(coils,1)
r=load(fullfile(src,[coils{c,2} '.mat'])); r=r.res;
%terminals at the +z (bore) end of the former, the outermost turn there
bore=[1;0;0];   %+z magnet = +x CoilGen
if max(bore'*r.wire_path(:,[1 end]))<max(bore'*r.wire_path)-0.001
r=move_wire_terminals(r,bore,r.wire.outer_diameter+0.0008);
end
res=r; save(fullfile(out,'results',[coils{c,1} '_coil.mat']),'res','-v7.3');
v=r.wire_path;
p=[v(3,:);-v(2,:);v(1,:)]'*1e3;                     %magnet frame, mm
rad=hypot(p(:,1),p(:,2));
theta=mod(atan2d(p(:,2),p(:,1)),360);
od=r.wire.outer_diameter*1e3;
depth=coils{c,4}-rad;                               %wire centre below the former surface
crossing=rad<r.spec.radius*1e3-0.5;
T=table(p(:,1),p(:,2),p(:,3),rad,theta,depth,crossing, ...
    'VariableNames',{'x_mm','y_mm','z_mm','r_mm','theta_deg','centre_depth_mm','crossing'});
writetable(T,fullfile(out,'wire_paths',[coils{c,1} '_wire_path.csv']));
f=r.field;
gap=norm(p(1,:)-p(end,:));
runs=find(diff([0;crossing])==1);
fprintf(fid,'%s coil (%s)\n',coils{c,1},coils{c,2});
fprintf(fid,'  former: inner radius %g, outer radius %g (ID %g / OD %g), length %g\n', ...
    coils{c,3},coils{c,4},2*coils{c,3},2*coils{c,4},r.spec.length*1e3);
fprintf(fid,'  wire: %s, max outer diameter %.3f, bare %.3f, pitch min %.2f\n',r.wire.name,od,r.wire.bare_radius*2e3,od+0.8);
fprintf(fid,'  wire centre radius %.2f (diameter %.2f), groove depth %.2f (wire flush with the surface)\n', ...
    r.spec.radius*1e3,r.spec.radius*2e3,od);
fprintf(fid,'  crossings: %d, wire centre %.2f deeper there (groove depth %.2f, %.2f plastic left inside)\n', ...
    numel(runs),-r.spec.bridge_lift*1e3,max(depth)+od/2,coils{c,4}-coils{c,3}-(max(depth)+od/2));
fprintf(fid,'  wire extent along z: %.1f to %.1f\n',min(p(:,3)),max(p(:,3)));
fprintf(fid,'  terminals: start (%.1f, %.1f, %.1f), end (%.1f, %.1f, %.1f), %.1f apart',p(1,:),p(end,:),gap);
if isfield(r,'terminals'), fprintf(fid,' (moved to the +z end; the old cut is closed by a %.1f mm jumper)',r.terminals.jumper_length_m*1e3); end
fprintf(fid,'\n');
fprintf(fid,'  turns (contour levels) %d, wire length %.1f m (+ leads)\n',r.spec.levels,r.wire_length);
fprintf(fid,'  efficiency %.3f mT/m/A, %.2f mT/m per sqrt(W)\n',f.efficiency_mT_per_m_per_A,f.efficiency_mT_per_m_per_A/sqrt(r.R20_ohm));
fprintf(fid,'  R %.3f ohm at 20 C, %.3f ohm at 80 C, L %.0f uH\n',r.R20_ohm,r.R80_ohm,r.inductance_H*1e6);
fprintf(fid,'  closest approach of two wire parts %.2f centre to centre\n',r.spacing.min_mm);
fprintf(fid,'  max nonlinearity on the 120 mm sphere %.2f %%\n\n',f.nonlinearity_max_percent);
end
fclose(fid);
copyfile(fullfile(src,'layout_summary.txt'),fullfile(out,'all_layouts_0p8mm.txt'));
copyfile(fullfile(src,'..','sf_summary.txt'),fullfile(out,'length_study_1mm.txt'));
type(fullfile(out,'design_summary.txt'));
