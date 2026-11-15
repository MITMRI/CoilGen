%This script designs a 3-axis gradient set for a Halbach magnet, wound from
%round wire in one layer in 3D printed formers, and compares wire gauges:
%32 cm clear inside the formers, 30 cm spherical imaging volume, 1 mm between
%neighbouring wires.
%
%Frame of the magnet: B0 along x (Halbach), bore along z. CoilGen optimizes the
%z-component of the field, so the CoilGen frame has B0 = z_cg, bore = x_cg:
%   x_magnet = z_cg,  y_magnet = -y_cg,  z_magnet = x_cg
%
%Each coil is saved as soon as it is finished; finished coils are not
%recomputed. Set the environment variable WOUND_JOB to "Z,AWG16" to run a
%single coil in its own MATLAB process.

if ~exist('skip_cd','var')
if ispc
cd('..\');
else
cd('../');
end
end
addpath(fullfile(pwd,'sub_functions'),fullfile(pwd,'verify'));

result_folder=fullfile(pwd,'Wound_gradient_coils','halbach_320mm_bore');
if ~exist(fullfile(result_folder,'gauge_study'),'dir'), mkdir(fullfile(result_folder,'gauge_study')); end

%       name shape  gradient dir  length[m] radius[m] tikonov   levels at 16 AWG
coils={'Z',  'x',  [1;0;0],      0.90,     0.163,    31622.7766, 114; ...
       'Y',  'y',  [0;1;0],      0.60,     0.171,    1e4,        66; ...
       'X',  'z',  [0;0;1],      0.60,     0.179,    1e3,        66};

clearance=0.001;      %between neighbouring wires, m
roi_radius=0.15;      %30 cm sphere
bridge_lift=0.0035;   %the wires that cross over the turns

%copper magnet wire, heavy build (MWS data, NEMA MW1000-2020): nominal bare
%diameter and resistance, maximum overall diameter. inch and ohm/1000 ft.
inch=0.0254; per_kft=1/304.8;
wires(1)=solid_wire('AWG12',0.0808*inch,0.0851*inch,1.589*per_kft);
wires(2)=solid_wire('AWG14',0.0641*inch,0.0682*inch,2.524*per_kft);
wires(3)=solid_wire('AWG16',0.0508*inch,0.0545*inch,4.019*per_kft);
wires(4)=solid_wire('AWG18',0.0403*inch,0.0437*inch,6.386*per_kft);
%16 AWG litz, 3x54/38 unserved single build, 0.071 inch overall (Remington
%Industries): 162 strands of 38 AWG, 648.2 ohm/1000 ft per strand, 3 % lay
lay_factor=1.03;
wires(5)=struct('name','LITZ16','type','litz','outer_diameter',0.071*inch, ...
    'bare_radius',0.0040*inch/2,'outer_radius_for_L',0.071*inch/2, ...
    'ohm_per_m_20C',648.2*per_kft/162*lay_factor, ...
    'strand_radius',0.0040*inch/2,'num_strands',162,'lay_factor',lay_factor);

base_pitch=wires(3).outer_diameter+clearance;

job=getenv('WOUND_JOB');
for c=1:size(coils,1)
for w=1:numel(wires)
tag=[coils{c,1} ',' wires(w).name];
if ~isempty(job) && ~strcmp(job,tag), continue; end
spec=struct('name',coils{c,1},'shape',coils{c,2},'gradient_dir',coils{c,3}, ...
    'length',coils{c,4},'radius',coils{c,5},'tikonov',coils{c,6},'roi_radius',roi_radius, ...
    'bridge_lift',bridge_lift,'make_sweep',false);
pitch=wires(w).outer_diameter+clearance;
spec.levels=round(coils{c,7}*base_pitch/pitch);
result_file=fullfile(result_folder,'gauge_study',[coils{c,1} '_' wires(w).name '.mat']);
res=design_wound_coil(spec,wires(w),result_file);
f=res.field;
fprintf('DONE %-8s levels %3d: %6.1f m, %.3f mT/m/A, nonlin %.2f %%, L %.3f mH, R80 %.3f ohm, closest %.2f mm\n', ...
    strrep(tag,',',' '),spec.levels,res.wire_length,f.efficiency_mT_per_m_per_A, ...
    f.nonlinearity_max_percent,res.inductance_H*1e3,res.R80_ohm,res.spacing.min_mm);
end
end

function w=solid_wire(name,bare_diameter,outer_diameter,ohm_per_m)
w=struct('name',name,'type','solid','outer_diameter',outer_diameter, ...
    'bare_radius',bare_diameter/2,'outer_radius_for_L',bare_diameter/2, ...
    'ohm_per_m_20C',ohm_per_m,'strand_radius',[],'num_strands',[],'lay_factor',[]);
end
