function res=design_wound_coil(spec,wire,result_file)
%Design one single-layer gradient coil wound from round wire in a printed
%former, then evaluate the wire path that comes out of CoilGen.
%
%spec:  .name .shape .gradient_dir .length .radius .tikonov .levels .roi_radius
%       .bridge_lift (height of the wires that cross over the turns, m)
%       .make_sweep  (true: also write the swept wire body as STL)
%wire:  .name .type .outer_diameter .bare_radius .ohm_per_m_20C and, for litz,
%       .strand_radius .num_strands .lay_factor (see wire_ac_resistance)
%result_file: .mat file for the result; nothing is recomputed if it exists.
%
%CoilGen frame: B0 = z, bore along x. The field is checked with Biot-Savart on
%points that the optimization did not use, because a stream function can
%overfit its target points.

if exist(result_file,'file'), res=load(result_file); res=res.res; return; end

num_cross_points=13;
angle=linspace(0,2*pi,num_cross_points);
cross_section=[sin(angle);cos(angle)]*wire.outer_diameter/2;
num_axial=max(30,round(spec.length/0.016));
output_folder=fileparts(result_file);

coil_out=CoilGen(...
    'field_shape_function',spec.shape,...
    'coil_mesh_file','create cylinder mesh',...
    'cylinder_mesh_parameter_list',[spec.length spec.radius 64 num_axial 0 1 0 pi/2],...
    'surface_is_cylinder_flag',true,...
    'target_mesh_file','sphere_radius150mm.stl',...
    'target_region_resolution',5,...
    'levels',spec.levels,...
    'pot_offset_factor',0.25,...
    'tikonov_reg_factor',spec.tikonov,...
    'min_loop_signifcance',0.5,...
    'interconnection_cut_width',0.012,...
    'normal_shift_length',spec.bridge_lift,...
    'normal_shift_smooth_factors',[3 3 3],...
    'smooth_flag',true,...
    'smooth_factor',1,...
    'skip_postprocessing',false,...
    'skip_sweep',~spec.make_sweep,...
    'save_stl_flag',spec.make_sweep,...
    'cross_sectional_points',cross_section,...
    'conductor_cross_section_width',wire.outer_diameter,...
    'conductor_cross_section_height',wire.outer_diameter,...
    'plot_flag',false,...
    'skip_inductance_calculation',true,...
    'output_directory',output_folder);

%CoilGen removes these from the path
addpath(fullfile(pwd,'sub_functions'),fullfile(pwd,'verify'));

wire_path=coil_out.coil_parts(1).wire_path.v;
res.spec=spec; res.wire=wire;
res.wire_path=wire_path;
res.shift_array=coil_out.coil_parts(1).shift_array;
res.wire_length=sum(vecnorm(diff(wire_path,1,2)));

res.field=validate_field(wire_path,spec.gradient_dir,spec.roi_radius);
res.spacing=wire_spacing(wire_path,wire.outer_diameter+0.001);

%inductance: uniform current in the round conductor (for litz: the bundle)
res.inductance_H=wire_inductance(wire_path,wire.outer_radius_for_L);

%resistance
res.R20_ohm=wire.ohm_per_m_20C*res.wire_length;
res.R_of_T=@(T) res.R20_ohm*(1+0.00393*(T-20));
res.R60_ohm=res.R_of_T(60);
res.R80_ohm=res.R_of_T(80);
res=rmfield(res,'R_of_T');
res.freqs=[100 1e3 3e3 1e4];
ac_model=struct('type',wire.type,'radius',wire.bare_radius);
if strcmp(wire.type,'litz')
ac_model.radius=wire.strand_radius;
ac_model.num_strands=wire.num_strands;
ac_model.lay_factor=wire.lay_factor;
ac_model.bundle_radius=wire.outer_diameter/2;
end
res.ac=wire_ac_resistance(wire_path,ac_model,res.freqs,80);

%keep the result in the form the CoilGen plots use
res.coil_layouts=struct('out',strip_coilgen_result(coil_out));
save(result_file,'res','-v7.3');
end

function r=validate_field(wire_path,gradient_dir,roi_radius)
%field of the wire at 1 A on a dense set of points on and inside the sphere
num_surface=4000;
n=(0:num_surface-1)'+0.5;
polar=acos(1-2*n/num_surface); azimuth=pi*(1+sqrt(5))*n;
surface_points=roi_radius*[cos(azimuth).*sin(polar) sin(azimuth).*sin(polar) cos(polar)];
[gx,gy,gz]=ndgrid(linspace(-1,1,17)*roi_radius);
inner=[gx(:) gy(:) gz(:)]; inner=inner(vecnorm(inner,2,2)<0.999*roi_radius,:);
points=[surface_points;inner];
A=wire_path(:,1:end-1)'; B=wire_path(:,2:end)';
current=ones(size(A,1),1);
b=biot_savart(A,B,current,points);
along=points*gradient_dir(:);
design=[along ones(numel(along),1)];
fit=design\b(:,3);
residual=b(:,3)-design*fit;
step=2e-3;
centre=biot_savart(A,B,current,[gradient_dir(:)';-gradient_dir(:)']*step);
r.efficiency_mT_per_m_per_A=abs(centre(1,3)-centre(2,3))/(2*step)*1e3;
r.fit_efficiency_mT_per_m_per_A=abs(fit(1))*1e3;
r.nonlinearity_max_percent=max(abs(residual))/abs(fit(1))/roi_radius*100;
r.nonlinearity_rms_percent=rms(residual)/abs(fit(1))/roi_radius*100;
r.max_position_error_mm=max(abs(residual))/abs(fit(1))*1e3;
r.max_field_mT_per_A=max(vecnorm(b,2,2))*1e3;
end

function r=wire_spacing(wire_path,nominal_pitch)
%closest approach between parts of the wire more than 30 mm apart along it
s=[0 cumsum(vecnorm(diff(wire_path,1,2)))];
sample_step=0.0008;
num_samples=ceil(s(end)/sample_step);
[s_unique,unique_inds]=unique(s);
s_new=linspace(0,s(end),num_samples)';
points=interp1(s_unique,wire_path(:,unique_inds)',s_new);
neighbours=rangesearch(points,points,nominal_pitch*1.5);
closest=inf(num_samples,1);
for i=1:num_samples
j=neighbours{i};
j=j(abs(s_new(j)-s_new(i))>0.03);
if ~isempty(j), closest(i)=min(vecnorm(points(j,:)-points(i,:),2,2)); end
end
valid=isfinite(closest);
r.min_mm=min(closest)*1e3;
r.percentile1_mm=prctile(closest(valid),1)*1e3;
r.wire_below_pitch_m=sum(closest<nominal_pitch*0.999)*sample_step;
r.closest_mm=closest*1e3;
r.points=points;
end
