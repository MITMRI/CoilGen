function res=evaluate_wound_wire(res)
%Field, wire spacing, inductance and resistance of a wound wire path.
%res: .wire_path (3xN, m, CoilGen frame) .spec (.gradient_dir .roi_radius)
%     .wire (see design_wound_coil). The results are added to res.

res.wire_length=sum(vecnorm(diff(res.wire_path,1,2)));

res.field=validate_field(res.wire_path,res.spec.gradient_dir,res.spec.roi_radius);
res.spacing=wire_spacing(res.wire_path,res.wire.outer_diameter+0.001);

%inductance: uniform current in the round conductor (for litz: the bundle)
res.inductance_H=wire_inductance(res.wire_path,res.wire.outer_radius_for_L);

%resistance
res.R20_ohm=res.wire.ohm_per_m_20C*res.wire_length;
res.R_of_T=@(T) res.R20_ohm*(1+0.00393*(T-20));
res.R60_ohm=res.R_of_T(60);
res.R80_ohm=res.R_of_T(80);
res=rmfield(res,'R_of_T');
res.freqs=[100 1e3 3e3 1e4];
ac_model=struct('type',res.wire.type,'radius',res.wire.bare_radius);
if strcmp(res.wire.type,'litz')
ac_model.radius=res.wire.strand_radius;
ac_model.num_strands=res.wire.num_strands;
ac_model.lay_factor=res.wire.lay_factor;
ac_model.bundle_radius=res.wire.outer_diameter/2;
end
res.ac=wire_ac_resistance(res.wire_path,ac_model,res.freqs,80);
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

