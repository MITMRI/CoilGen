function res=stack_crossings(res,lift,ramp)
%Crossings with the wires touching, for thin printed formers ("A+").
%
%CoilGen moves the wire that passes under other turns (the lead) inward by
%normal_shift_length. Here the lead is put directly under the turns instead,
%enamel on enamel, and the short pieces of turn that pass over a lead are
%lifted by `lift` (m) above their normal position, so the lead needs less
%depth: lead centre = wire radius + lift - wire diameter - 0.05 mm.
%The lift ramps to zero over `ramp` (m) along the turn on each side.
%The wire is evaluated again afterwards (evaluate_wound_wire).
%
%res: result of design_wound_coil (CoilGen frame, bore along x)

v=res.wire_path; od=res.wire.outer_diameter; R=res.spec.radius;
rad=vecnorm(v(2:3,:)); phi=atan2(v(3,:),v(2,:));
old_shift=max(R-rad,0);
is_lead=old_shift>0.5e-3;
nominal=abs(res.spec.bridge_lift);                %CoilGen's full shift
frac=min(old_shift/(0.95*nominal),1);              %0..1, keeps CoilGen's ramps
margin=0.05e-3;                                    %so the wires just touch
new_rad=R-frac*(od+margin-lift);

%turn pieces over a lead: distance to the lead on the unrolled surface
uv=[v(1,:); phi*R];
lead_uv=uv(:,is_lead & frac>0.3);
d=inf(1,size(v,2));
if ~isempty(lead_uv)
%both sides of the phi wrap
cand=[lead_uv lead_uv+[0;2*pi*R] lead_uv-[0;2*pi*R]];
[~,dist]=knnsearch(cand',uv');
d=dist';
end
full=od+1e-3;                                      %full lift this close to a lead
up=lift*min(max((full+ramp-d)/ramp,0),1);
%not the wire running into or out of a lead (within 15 mm along the wire)
s=[0 cumsum(vecnorm(diff(v,1,2)))];
near_lead=false(size(s));
lead_s=s(is_lead);
if ~isempty(lead_s)
[~,ds]=knnsearch(lead_s',s'); near_lead=ds'<0.015;
end
up(is_lead|near_lead)=0;
new_rad=new_rad+up;

scale=new_rad./rad;
res.wire_path(2:3,:)=v(2:3,:).*scale;
res.crossing_stack=struct('lift_m',lift,'ramp_m',ramp,'is_lead',is_lead,'lift_profile_m',up);
res=evaluate_wound_wire(res);
end
