%Export the chosen 12 cm RF solenoid to designs/rf_12cm for CAD: litz
%centreline (CSV), turn positions, return and the two winding capacitors, a
%summary and figures. Run rf_solenoid_sweep.m first. Set N before running
%to export another turn count (default 12). See designs/rf_12cm/README.md.
%
%Magnet frame: x = B0, z = bore, theta from +x toward +y. Lengths in mm in the files.

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
out=fullfile(pwd,'designs','rf_12cm');
for d={'figures'}, if ~exist(fullfile(out,d{1}),'dir'), mkdir(fullfile(out,d{1})); end, end
if ~exist('N','var'), N=18; end
S=load(fullfile(out,'sweep.mat')); spec=rf_12cm_spec();
i=find(cellfun(@(d) d.N==N && abs(d.L-0.160)<1e-9,S.designs),1);
z=S.designs{i}.z;

%return: copper foil under the groove floor, under the seam (the 120 mm sphere
%touches the inner wall, so nothing can go in the bore); seam and return at
%theta = 90 deg, where the return's near field is mostly along B0
ret_opts={'return_type','strip','seam_deg',90,'pocket_halfwidth_mm',0.85}; %caps stacked: 1.6 wide, fits the groove
%widen the jogs until neighbouring wires keep the minimum pitch everywhere
for jog_deg=30:5:90
g=build_rf_solenoid_path(z,spec,'jog_deg',jog_deg,ret_opts{:});
[dmin,where]=min_wire_distance(g.litz,0.015);
if dmin>=spec.pitch_min-1e-5, break, end
end
r=evaluate_rf_solenoid(g,spec);
r.N=N; r.z=z; r.g=g; r.coil_length=max(z)-min(z)+spec.litz_od; r.min_wire_distance=dmin;
res=r; save(fullfile(out,sprintf('rf_solenoid_N%d.mat',N)),'res','spec');
bore=strcmp(g.return_type,'bore');
d360=@(t) mod(rad2deg(t),360);

%--- litz centreline (coil, then the return) ---
p=[g.litz g.strip_centre]'*1e3;
rad=hypot(p(:,1),p(:,2)); th=mod(atan2d(p(:,2),p(:,1)),360);
names={'ring','jog','lead'};
part=[names(g.part) repmat({'return'},1,size(g.strip_centre,2))]';
T=table(p(:,1),p(:,2),p(:,3),rad,th,spec.former_od/2*1e3-rad,part, ...
    'VariableNames',{'x_mm','y_mm','z_mm','r_mm','theta_deg','depth_below_OD_mm','part'});
writetable(T,fullfile(out,'litz_path.csv'));
writetable(table((1:N)',sort(z,'descend')'*1e3,'VariableNames',{'turn','z_mm'}),fullfile(out,'turns.csv'));

%--- capacitors in the winding (the loop's 3rd equal share is the port's tune/match network) ---
nb=spec.num_breaks; ic=2:nb; nc=numel(ic);
C=r.C_break_F*1e12;
Tb=table((1:nc)',d360(g.break_theta(ic))',g.break_z(ic)'*1e3,spec.a*1e3*cos(g.break_theta(ic))', ...
    spec.a*1e3*sin(g.break_theta(ic))',g.break_s(ic)'*1e3,repmat(C,nc,1),repmat(C/2,nc,1), ...
    'VariableNames',{'cap','theta_deg','z_mm','x_mm','y_mm','litz_arc_from_port_mm','C_total_pF','C_each_1206_pF'});
writetable(Tb,fullfile(out,'caps.csv'));
if exist(fullfile(out,'breaks.csv'),'file'), delete(fullfile(out,'breaks.csv')); end

%--- summary ---
fid=fopen(fullfile(out,'design_summary.txt'),'w');
w=2*pi*spec.f0;
fprintf(fid,'12 cm Halbach RF solenoid, %d turns. B0 = %.0f mT along x, f0 = %.4f MHz. Bore = z.\n',N,spec.B0*1e3,spec.f0/1e6);
fprintf(fid,'Magnet frame: x = B0, z = bore, theta from +x toward +y. Lengths in mm.\n\n');
fprintf(fid,'Former: ID %.0f, OD %.0f. Litz %d x 48 AWG, %.2f OD, in a groove on the outside, flush:\n', ...
    spec.former_id*1e3,spec.former_od*1e3,spec.litz_strands,spec.litz_od*1e3);
fprintf(fid,'  turn (centre) radius %.2f, groove depth %.2f (width ~%.1f), %.2f of plastic under the groove\n', ...
    spec.a*1e3,spec.litz_od*1e3,spec.litz_od*1e3+0.1,(spec.former_od-spec.former_id)/2*1e3-spec.litz_od*1e3);
fprintf(fid,'  min pitch %.2f (0.8 plastic between grooves); closest wires in the path %.2f\n',spec.pitch_min*1e3,dmin*1e3);
fprintf(fid,'Turns (z of the ring centre, wound from +z to -z, counter-clockwise seen from +z):\n ');
fprintf(fid,' %.2f',sort(z,'descend')*1e3); fprintf(fid,'\n');
pz=-diff(sort(z,'descend'))*1e3;
fprintf(fid,'  pitches: '); fprintf(fid,' %.2f',pz); fprintf(fid,'\n');
fprintf(fid,'  coil length over the litz %.1f; leads to %.1f / %.1f\n',r.coil_length*1e3,g.z_top*1e3,g.z_bottom*1e3);
fprintf(fid,'Jogs: each turn steps to the next across the seam at theta = %.0f, over %.0f deg (%.1f of arc), smooth (raised cosine) in z\n', ...
    d360(g.seam_theta),g.options.jog_deg,deg2rad(g.options.jog_deg)*spec.a*1e3);
fprintf(fid,'  top lead: the litz starts at theta %.1f, z %.1f and runs axially to the first ring\n',d360(g.feed_theta(2)),g.z_top*1e3);
fprintf(fid,'  bottom lead: the litz leaves the last ring at theta %.1f and runs axially to z %.1f\n',d360(g.strip_theta_end),g.z_bottom*1e3);
if bore
rr=g.options.return_radius*1e3;
fprintf(fid,'Return: the same litz, inside the bore, centre at r %.2f, in a groove %.2f deep in the inner wall\n',rr,rr+spec.litz_od*500-spec.former_id*500);
fprintf(fid,'  (sticks %.2f into the bore; %.2f of plastic between it and the outer grooves), %.0f long:\n', ...
    spec.former_id*500-(rr-spec.litz_od*500),(spec.former_od/2-spec.litz_od)*1e3-(rr+spec.litz_od*500),g.strip_length*1e3);
fprintf(fid,'  hole through the wall at theta %.1f, z %.1f (end of the bottom lead); inside: along z = %.1f to theta %.0f,\n', ...
    d360(g.strip_theta_end),g.z_bottom*1e3,g.z_bottom*1e3,d360(g.seam_theta));
fprintf(fid,'  axially at theta %.0f from z %.1f to %.1f (under the seam), along z = %.1f to theta %.1f, out through a hole there\n', ...
    d360(g.seam_theta),g.z_bottom*1e3,g.z_top*1e3,g.z_top*1e3,d360(g.feed_theta(1)));
else
fprintf(fid,'Return: copper foil %.0f wide (2 oz, >= 0.07 thick) insulated with Kapton, %.0f long, in a slot whose floor is %.2f below the OD\n', ...
    g.options.strip_width_mm,g.strip_length*1e3,spec.litz_od*1e3+0.15);
fprintf(fid,'  (0.15 under the groove floor, so the jogs over it stay flush; %.2f of plastic left inside)\n', ...
    (spec.former_od-spec.former_id)/2*1e3-spec.litz_od*1e3-0.15);
fprintf(fid,'  bottom tab along z = %.1f from theta %.1f (end of the bottom lead) to %.0f; axial at theta %.0f from z %.1f to %.1f\n', ...
    g.z_bottom*1e3,d360(g.strip_theta_end),d360(g.seam_theta),d360(g.seam_theta),g.z_bottom*1e3,g.z_top*1e3);
fprintf(fid,'  (under the jogs); top tab along z = %.1f from theta %.0f to %.1f (port)\n',g.z_top*1e3,d360(g.seam_theta),d360(g.feed_theta(1)));
end
fprintf(fid,'Port (to the tune/match board): at the +z end, z %.1f, between the litz start (theta %.1f) and the return (theta %.1f)\n', ...
    g.z_top*1e3,d360(g.feed_theta(2)),d360(g.feed_theta(1)));
fprintf(fid,'Capacitors in the winding (%d, each 2 x 1206 in parallel, stacked, in the groove widened to %.1f over %.1f along the wire, ~%.1f deep):\n', ...
    nc,2*g.options.pocket_halfwidth_mm,2*g.options.pocket_halflength_mm,spec.litz_od*1e3);
for k=1:nc
fprintf(fid,'  C%d: theta %6.2f deg, z %7.2f\n',k,d360(g.break_theta(ic(k))),g.break_z(ic(k))*1e3);
end
fprintf(fid,'  chunk lengths (port -> C1 -> C2 -> return -> port): %s m\n',sprintf(' %.3f',g.chunk_length));
fprintf(fid,'\nElectrical (unloaded, 20 C, at f0):\n');
fprintf(fid,'  L %.2f uH, R %.3f ohm (litz dc %.3f, litz ac %.3f, return %.3f), Q %.0f\n', ...
    r.L_H*1e6,r.R_ac,r.R_litz_dc,r.R_litz_ac,r.R_strip,r.Q_unloaded);
fprintf(fid,'  litz length %.2f m coil + %.2f m return\n',g.litz_length,g.strip_length);
fprintf(fid,'  B1 (perpendicular to B0) at the centre %.1f uT per A; %.0f uT per sqrt(W) of coil loss\n', ...
    r.B1_centre_T_per_A*1e6,r.B1_uT_per_sqrtW);
fprintf(fid,'  tuning: %.0f pF in series around the loop in total, split in %d equal shares of %.0f pF:\n',r.C_total_F*1e12,nb,C);
fprintf(fid,'    C1..C%d in the winding = %.0f pF each (2 x %.0f pF 1206), the remaining share is the port network (tune + match ~ %.0f pF)\n',nc,C,C/2,C);
fprintf(fid,'    %.1f V peak across each per A peak; with 2 caps in parallel each carries half the current\n',r.V_break_per_A);
fprintf(fid,'\nB1 of the real path (coil + return): %.1f %% of the 120 mm sphere volume within +-5 %% of the centre, %.1f %% within +-10 %%\n', ...
    r.vol_within_5*100,r.vol_within_10*100);
fprintf(fid,'Uniformity over a sphere:\n');
fprintf(fid,'  diameter    spread(max-min)/(max+min)   min..max vs centre\n');
for k=1:numel(r.sphere_d)
fprintf(fid,'  %5.0f        %6.2f %%                  %+6.2f .. %+6.2f %%\n',r.sphere_d(k)*1e3,r.spread(k)*100,r.dev_min(k)*100,r.dev_max(k)*100);
end
fclose(fid);
type(fullfile(out,'design_summary.txt'))

%--- figures ---
f=figure('Visible','off','Position',[0 0 1400 700]);
a1=spec.a*1e3;
s_arc=@(t) t*a1;                                  %unrolled: arc length at the turn radius
rel=@(t) mod(t-g.seam_theta+pi,2*pi)-pi;          %angle from the seam
tt=rel(g.theta);
br=find(abs(diff(tt))>pi);
xs=s_arc(tt); xs(br+1)=nan;                       %do not join across the cut
hold on
patch([-pi pi pi -pi]*a1,[-1 -1 1 1]*spec.roi_radius*1e3,[0.85 0.93 1],'EdgeColor','none');
text(-pi*a1+5,-spec.roi_radius*1e3+12,sprintf('%.0f mm sphere (axial extent)',2*spec.roi_radius*1e3),'Color',[0.3 0.4 0.7]);
rt=rel(atan2(g.strip_centre(2,:),g.strip_centre(1,:)));
plot(s_arc(rt),g.strip_centre(3,:)*1e3,'k-','LineWidth',3);
plot(xs,g.z*1e3,'Color',[0.8 0.45 0.1],'LineWidth',2);
bt=rel(g.break_theta(ic));
plot(s_arc(bt),g.break_z(ic)*1e3,'s','MarkerSize',12,'MarkerFaceColor','y','MarkerEdgeColor','k');
for k=1:nc, text(s_arc(bt(k))+4,g.break_z(ic(k))*1e3+3,sprintf('C%d',k),'FontWeight','bold'); end
pf=mean(rel(g.feed_theta));
plot(s_arc(pf),g.z_top*1e3,'o','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');
text(s_arc(pf)+4,g.z_top*1e3+4,'port','FontWeight','bold');
axis equal; xlim([-pi pi]*a1); grid on
xlabel(sprintf('arc length at r = %.2f mm from the seam (theta = %.0f deg), mm',a1,d360(g.seam_theta))); ylabel('z (bore), mm');
title(sprintf('RF solenoid, %d turns of %.1f mm litz, unrolled (orange: litz, black: return foil under the seam, yellow: caps, green: port)',N,spec.litz_od*1e3));
exportgraphics(f,fullfile(out,'figures','unrolled.png'),'Resolution',150); close(f);

%B1 maps in the planes through the axis: x-z (along B0) and y-z, wires to scale
A=g.litz(:,1:end-1)'; B=g.litz(:,2:end)'; I=ones(size(A,1),1);
for k=1:numel(g.strip_filaments)
fl=g.strip_filaments{k}; A=[A;fl(:,1:end-1)']; B=[B;fl(:,2:end)']; I=[I;ones(size(fl,2)-1,1)/numel(g.strip_filaments)]; %#ok<AGROW>
end
allw=[g.litz g.strip_centre];
allw=densify(allw,0.0003);
[u,v]=meshgrid(-0.068:0.0005:0.068,-0.09:0.0005:0.09);
f=figure('Visible','off','Position',[0 0 1300 850]);
names_pl={'x-z plane (y = 0)','y-z plane (x = 0)'};
rw=spec.litz_od/2*1e3;
for k=1:2
P=zeros(numel(u),3); P(:,3)=v(:);
if k==1, P(:,1)=u(:); else, P(:,2)=u(:); end
Bp=biot_savart(A,B,I,P);
b=reshape(hypot(Bp(:,2),Bp(:,3))/r.B1_centre_T_per_A,size(u));
[~,dw]=knnsearch(allw',P);
b(reshape(dw,size(u))<spec.litz_od/2)=nan;        %inside a wire
b(abs(u)>spec.former_od/2)=nan;
subplot(1,2,k);
contourf(u*1e3,v*1e3,(b-1)*100,[-20 -10 -5 -2 -1 1 2 5 10 20],'ShowText','on'); hold on
caxis([-10 10]); colormap(gca,turbo); colorbar
t=linspace(0,2*pi,200); plot(55*cos(t),55*sin(t),'k:',60*cos(t),60*sin(t),'k--');
%former wall
for sgn=[-1 1]
patch(sgn*[spec.former_id spec.former_od spec.former_od spec.former_id]*500,[-1 -1 1 1]*g.z_top*1e3*1.08, ...
    [0.6 0.6 0.6],'FaceAlpha',0.25,'EdgeColor','none');
end
%wire cross-sections: where the path crosses the plane (circles), or runs in it (bars)
c=k; o=3-k;                                       %coordinate normal to the plane, in-plane one
W=allw*1e3; 
sg=sign(W(o,:)+1e-9);                            %points exactly on the plane count as one side
cr=find(sg(1:end-1)~=sg(2:end) & abs(W(o,1:end-1)-W(o,2:end))<5);
for j=cr
s=W(o,j)/(W(o,j)-W(o,j+1)); q=W(:,j)+s*(W(:,j+1)-W(:,j));
rectangle('Position',[q(c)-rw q(3)-rw 2*rw 2*rw],'Curvature',[1 1],'FaceColor',[0.85 0.5 0.15],'EdgeColor','k');
end
Wr=densify(g.strip_centre,0.0003)*1e3;           %only the return can lie in a plane
inpl=abs(Wr(o,:))<1e-3 & abs(diff([Wr(3,1) Wr(3,:)]))>0;
if any(inpl)
qi=Wr(:,inpl);
for sgn=unique(sign(qi(c,:)))
qq=qi(:,sign(qi(c,:))==sgn);
wr=rw; if ~bore, wr=0.15; end                 %foil seen edge-on
rectangle('Position',[mean(qq(c,:))-wr min(qq(3,:))-wr 2*wr max(qq(3,:))-min(qq(3,:))+2*wr],'FaceColor','k','EdgeColor','k');
end
end
axis equal; xlim([-68 68]); ylim([-90 90]); xlabel('mm'); ylabel('z, mm');
title({['B1 deviation from centre, %, ' names_pl{k}],'dashed 120 mm, dotted 110 mm sphere; litz to scale, grey = former'});
end
exportgraphics(f,fullfile(out,'figures','B1_maps.png'),'Resolution',200); close(f);
fprintf('closest wires %.2f mm (jog %d deg)\n',dmin*1e3,jog_deg);

function Q=densify(P,h)
%points along the polyline P at most h apart
s=[0 cumsum(vecnorm(diff(P,1,2)))];
keep=[true diff(s)>0]; P=P(:,keep); s=s(keep);
Q=interp1(s,P',0:h:s(end))';
end

function [dmin,where]=min_wire_distance(P,skip)
%closest approach of the path to itself, ignoring points within skip (m) of
%arc length of each other
s=[0 cumsum(vecnorm(diff(P,1,2)))];
dmin=inf; where=0;
for i=1:numel(s)
far=abs(s-s(i))>skip;
d=sqrt(sum((P(:,far)-P(:,i)).^2,1));
[m,j]=min(d);
if m<dmin, dmin=m; where=i; end
end
end
