%Wire layouts for the 12 cm Halbach gradient set from the stream functions of
%halbach_12cm_sf_sweep.m: one layer of round wire in grooves on the outside of
%the printed formers, 1 mm (or WOUND_CLEARANCE_MM) between wires, 120 mm spherical imaging volume.
%   Z (bore) on the 66-70 mm former, Y on 71-75 mm, X (along B0) on 76-80 mm
%The wire centre sits half a wire diameter below the outer former surface; the
%wires that cross over the turns go further into the former.
%CoilGen frame: B0 = z_cg, bore = x_cg (x_magnet = z_cg, y_magnet = -y_cg,
%z_magnet = x_cg).
%
%For each coil, length, wire and linearity limit, the power weight with the most
%gradient per watt within the limit is chosen, and the number of contour levels
%is the largest that keeps neighbouring wires one pitch apart at the densest
%place. Set WOUND_JOB to "Z,340,AWG14,0.20" (or a ';' separated list).

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
threads=getenv('WOUND_THREADS'); if ~isempty(threads), maxNumCompThreads(str2double(threads)); end

sweep_folder=fullfile(pwd,'Wound_gradient_coils','halbach_12cm','sf_sweep');
%plastic between neighbouring wires: WOUND_CLEARANCE_MM (default 1 mm -> 'layouts')
clearance_mm=str2double(getenv('WOUND_CLEARANCE_MM')); if isnan(clearance_mm), clearance_mm=1; end
clearance=clearance_mm/1e3; pot_offset=0.5;
result_folder=fullfile(pwd,'Wound_gradient_coils','halbach_12cm',layout_folder_name(clearance_mm));
if ~exist(result_folder,'dir'), mkdir(result_folder); end

%       name shape dir      former outer radius  lengths [mm]
coils={'Z', 'x', [1;0;0], 0.070, [340 330]; ...
       'Y', 'y', [0;1;0], 0.075, [300 324]; ...
       'X', 'z', [0;0;1], 0.080, 318};   %maximum length: Z 330, Y 324, X 318 mm
limits=[0.05 0.2 1 2 3 5];  %max nonlinearity of the stream function on the 120 mm sphere, %

%copper magnet wire, heavy build (MWS, NEMA MW1000): bare / max outer diameter
inch=0.0254; per_kft=1/304.8;
wires(1)=solid_wire('AWG14',0.0641*inch,0.0682*inch,2.524*per_kft);
wires(2)=solid_wire('AWG16',0.0508*inch,0.0545*inch,4.019*per_kft);
sweep_od=wires(1).outer_diameter;   %the sweep meshes sit at the 14 AWG wire centre

job=getenv('WOUND_JOB');
for c=1:size(coils,1)
for len=coils{c,5}
for w=1:numel(wires)
for lim=limits
tag=sprintf('%s,%d,%s,%.2f',coils{c,1},len,wires(w).name,lim);
if ~isempty(job) && ~any(strcmp(tag,strsplit(job,';'))), continue; end
name=sprintf('%s_L%03d_%s_lin%03d',coils{c,1},len,wires(w).name,round(lim*100));
result_file=fullfile(result_folder,[name '.mat']);
sw=load(fullfile(sweep_folder,sprintf('%s_L%03d.mat',coils{c,1},len))); sw=sw.res;
pitch=wires(w).outer_diameter+clearance;
dpsi=sw.K_max*pitch; levels=floor(sw.sf_range./dpsi-1+2*pot_offset);
levels=levels-mod(levels,2);       %even: no contour on the zero line
fom=levels.^2./(wires(w).ohm_per_m_20C*sw.K_int.*sw.sf_range); fom(sw.nonlin_max_pct>lim)=0;
[~,k]=max(fom);
%stream function on the mesh at this wire's centre radius
radius=coils{c,4}-wires(w).outer_diameter/2;
coil_mesh=sw.mesh;
coil_mesh.vertices(2:3,:)=coil_mesh.vertices(2:3,:)*radius/(coils{c,4}-sweep_od/2);
stream_function=sw.sf(:,k);
t=load(fullfile(pwd,'target_fields','halbach_12cm_targets.mat')); t=t.targets;
target_field=struct('coords',t.coords,'b',[zeros(2,size(t.coords,2)); t.(['g' coils{c,2}])]);
target_field.weights=ones(1,size(t.coords,2)); target_field.target_field_group_inds=ones(1,size(t.coords,2));
sf_file=['halbach_12cm_' name '.mat'];
save(fullfile(pwd,'Pre_Optimized_Solutions',sf_file),'coil_mesh','stream_function','target_field');
spec=struct('name',coils{c,1},'shape',coils{c,2},'gradient_dir',coils{c,3}, ...
    'length',len/1e3,'radius',radius,'roi_radius',0.060,'levels',levels(k), ...
    'sf_source_file',sf_file,'pot_offset_factor',pot_offset,'interconnection_cut_width',0.010, ...
    'bridge_lift',-(wires(w).outer_diameter+0.0003),'make_sweep',false, ...
    'mu',sw.mu(k),'sf_nonlin_max_pct',sw.nonlin_max_pct(k),'former_outer_radius',coils{c,4});
t0=tic;
res=design_wound_coil(spec,wires(w),result_file);
delete(fullfile(pwd,'Pre_Optimized_Solutions',sf_file));
f=res.field;
fprintf('DONE %-22s %3.0fs lev %3d: %5.1f m, %.3f mT/m/A, nonlin %.2f %% (sf %.2f %%), R20 %.3f ohm, L %.0f uH, closest %.2f mm, %.2f mT/m/sqrtW\n', ...
    name,toc(t0),spec.levels,res.wire_length,f.efficiency_mT_per_m_per_A,f.nonlinearity_max_percent,spec.sf_nonlin_max_pct, ...
    res.R20_ohm,res.inductance_H*1e6,res.spacing.min_mm,f.efficiency_mT_per_m_per_A/sqrt(res.R20_ohm));
end
end
end
end

function w=solid_wire(name,bare_diameter,outer_diameter,ohm_per_m)
w=struct('name',name,'type','solid','outer_diameter',outer_diameter, ...
    'bare_radius',bare_diameter/2,'outer_radius_for_L',bare_diameter/2, ...
    'ohm_per_m_20C',ohm_per_m,'strand_radius',[],'num_strands',[],'lay_factor',[]);
end

function name=layout_folder_name(clearance_mm)
if clearance_mm==1, name='layouts'; else, name=['layouts_' strrep(sprintf('%gmm',clearance_mm),'.','p')]; end
end
