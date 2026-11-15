function res=move_wire_terminals(res,end_dir,lead_gap)
%Move the two terminals of a wound coil to the end of the former.
%
%The wire of a CoilGen coil is one series loop opened at a single cut, so its
%two ends lie next to each other. Here the two ends are joined by a straight
%jumper, and the wire is opened instead where it reaches furthest along
%end_dir (the outermost turn at that end of the former), leaving lead_gap (m)
%between the new terminals. The jumper must not run into other turns; that is
%checked. The wire is then evaluated again (evaluate_wound_wire).
%
%res: result of design_wound_coil   end_dir: unit vector, CoilGen frame
%(e.g. [1;0;0] = +bore end)   lead_gap: distance between the terminals, m

v=res.wire_path;
jumper=v(:,end)-v(:,1);
t=linspace(0,1,50);
seg=v(:,1)+jumper.*t;
s=[0 cumsum(vecnorm(diff(v,1,2)))];
far=s>0.03 & s<s(end)-0.03;
clearance=min(min(pdist2(seg',v(:,far)')));
if clearance<res.wire.outer_diameter
error('move_wire_terminals: the jumper passes %.2f mm from another turn',clearance*1e3);
end
%close the loop, then open it at the point furthest along end_dir
closed=[v v(:,1)];
sc=[0 cumsum(vecnorm(diff(closed,1,2)))];
[~,i]=max(end_dir(:)'*v);
keep_before=sc<sc(i)-lead_gap/2;
keep_after=sc>sc(i)+lead_gap/2;
keep_after(end)=false;          %the closing point equals the first point
res.wire_path=[closed(:,keep_after) closed(:,keep_before)];
res.shift_array=[];
res.terminals=struct('moved',true,'jumper_length_m',norm(jumper),'jumper_clearance_m',clearance, ...
    'lead_gap_m',lead_gap,'end_dir',end_dir(:));
res=evaluate_wound_wire(res);
end
