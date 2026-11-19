function g=build_rf_solenoid_path(z,spec,varargin)
%BUILD_RF_SOLENOID_PATH  Litz path of a variable-pitch RF solenoid along z,
%with its return strip and capacitor breaks (all in m, magnet frame:
%x = B0, z = bore, theta from +x toward +y).
%
%Each turn k is a ring at z(k). The wire goes counter-clockwise (seen from +z),
%starts at the +z end and steps down to the next ring in a "jog" that spans
%jog_deg of arc centred on theta = 0 (the seam). After the last (-z) ring the
%wire drops to below the coil and joins a flat copper return strip that runs
%back to the +z end along theta = 0, in a shallow slot under the groove floor
%(under the jogs). Running it under the seam cancels the net axial current of
%the jogs locally.
%
%Breaks: the loop (litz + strip) is cut into num_breaks pieces of equal length.
%Break 1 is at the +z end between the litz start and the top of the strip
%(the feed/match point). The others are moved along the wire to the nearest
%place on a ring, away from the seam, with room for the capacitor pocket.
%
%Options: 'jog_deg' (30), 'lead_mm' (4, axial lead past the end rings),
%'feed_gap_mm' (6, arc between the strip and the litz start at the feed),
%'pocket_halfwidth_mm' (1.8, half width of the capacitor pocket across the wire),
%'pocket_halflength_mm' (4.0, half length along the wire), 'strip_width_mm' (6),
%'return_type' ('strip': foil under the seam; 'bore': litz on the inside wall of
%the former, reached through holes past the end turns), 'return_radius' (m, centre
%of the bore return, default 0.060), 'seam_deg' (0: angle of the seam and the
%return; the whole coil is rotated about z by this).

p=inputParser;
addParameter(p,'jog_deg',30); addParameter(p,'lead_mm',4);
addParameter(p,'feed_gap_mm',6); addParameter(p,'pocket_halfwidth_mm',1.8);
addParameter(p,'pocket_halflength_mm',4.0); addParameter(p,'strip_width_mm',6);
addParameter(p,'step_deg',1);
addParameter(p,'return_type','strip'); addParameter(p,'return_radius',0.060);
addParameter(p,'seam_deg',0);
parse(p,varargin{:}); o=p.Results;

a=spec.a;
z=sort(z(:)','descend');                         %wind from +z to -z
N=numel(z);
jog=deg2rad(o.jog_deg); dstep=deg2rad(o.step_deg);
delta=o.feed_gap_mm*1e-3/a;                      %arc between strip tab and litz at the feed
lead=o.lead_mm*1e-3;
t_start=jog/2;                                   %first ring starts and last ring ends clear of the seam
t_end=2*pi-jog/2;

th=[]; zz=[]; part=[];                           %part: 1 ring, 2 jog, 3 lead
%start lead: from the feed (above the top ring) down to the top ring
zt=z(1)+lead;
add(t_start*[1 1],[zt z(1)],3);
for k=1:N
t0=jog/2;
t1=2*pi-jog/2;
t=[t0:dstep:t1 t1]; t=unique(t);
add(t,z(k)*ones(size(t)),1);
if k<N                                           %jog to the next ring across the seam
t=linspace(2*pi-jog/2,2*pi+jog/2,max(3,ceil(jog/dstep)+1));
s=(t-t(1))/(t(end)-t(1)); s=(1-cos(pi*s))/2;   %smooth step
add(t,z(k)+(z(k+1)-z(k))*s,2);
end
end
%end lead: from the bottom ring down to below the coil
zb=z(N)-lead;
add([t_end t_end],[z(N) zb],3);
th=th(:)'; zz=zz(:)'; part=part(:)';
%remove repeated points
keep=[true hypot(diff(th)*a,diff(zz))>1e-9];
th=th(keep); zz=zz(keep); part=part(keep);
litz=[a*cos(th);a*sin(th);zz];

%return strip: flat copper under the groove floor. From the end of the litz
%(theta t_end, z zb) a tab runs at zb to theta 0, the strip runs up theta 0 to
%zt, and a tab at zt runs to theta t_start-delta, next to the litz start (feed).
r_strip=spec.former_od/2-spec.litz_od-0.0001;
if strcmp(o.return_type,'bore'), r_strip=o.return_radius; end
t_tab_top=t_start-delta;
nz=max(2,ceil((zt-zb)/0.002)+1);
tb=linspace(t_end-2*pi,0,max(2,ceil((2*pi-t_end)*a/0.002)+1));
tt=linspace(0,t_tab_top,max(2,ceil(t_tab_top*a/0.002)+1));
st=[tb zeros(1,nz) tt]; sz=[zb*ones(size(tb)) linspace(zb,zt,nz) zt*ones(size(tt))];
keep=[true hypot(diff(st)*r_strip,diff(sz))>1e-9]; st=st(keep); sz=sz(keep);
g.strip_centre=[r_strip*cos(st);r_strip*sin(st);sz];
if strcmp(o.return_type,'bore')
%radial pieces through the wall at both ends, litz as one filament
g.strip_centre=[[a*cos(st(1));a*sin(st(1));sz(1)] g.strip_centre [a*cos(st(end));a*sin(st(end));sz(end)]];
end
%modelled as 5 parallel filaments across its width, each carrying 1/5
%(offset sideways: in theta along the axial run, in z along the tabs)
w=o.strip_width_mm*1e-3; off=linspace(-w/2,w/2,5)*(4/5);
g.strip_filaments=cell(1,5);
axial=[false abs(diff(sz))>1e-9];
for i=1:5
phi=st+~axial.*0+axial.*off(i)/r_strip;
zf=sz+(~axial).*off(i).*sign(sz);
g.strip_filaments{i}=[r_strip*cos(phi);r_strip*sin(phi);zf];
end
if strcmp(o.return_type,'bore'), g.strip_filaments={g.strip_centre}; end
g.return_type=o.return_type;

g.litz=litz; g.theta=th; g.z=zz; g.part=part; g.z_rings=z;
g.litz_length=sum(vecnorm(diff(litz,1,2)));
g.strip_length=sum(vecnorm(diff(g.strip_centre,1,2)));
g.z_top=zt; g.z_bottom=zb; g.feed_theta=[t_tab_top t_start]; g.strip_theta_end=t_end;
g.options=o;

%capacitor breaks
s=[0 cumsum(vecnorm(diff(litz,1,2)))];
loop_len=g.litz_length+g.strip_length;
nb=spec.num_breaks;
g.break_s=zeros(1,nb); g.break_xyz=zeros(3,nb); g.break_theta=zeros(1,nb);
g.break_z=zeros(1,nb); g.break_ok=true(1,nb);
%break 1: the feed, between the top strip tab and the litz start
tf=t_start-delta/2;
g.break_s(1)=0; g.break_theta(1)=tf; g.break_z(1)=zt;
g.break_xyz(:,1)=[a*cos(tf);a*sin(tf);zt];
%room needed around a pocket: neighbouring ring centres at least this far in z
need_z=o.pocket_halfwidth_mm*1e-3+(spec.pitch_min-spec.litz_od)+spec.litz_od/2;
ring_gap=@(zr) min(abs(z(abs(z-zr)>1e-9)-zr));
cand=find(part==1);
ok=false(size(cand));
pocket_arc=o.pocket_halflength_mm*1e-3/a;
for i=1:numel(cand)
j=cand(i); tt=mod(th(j),2*pi);
far_from_seam=tt>jog/2+pocket_arc+deg2rad(3) && tt<2*pi-jog/2-pocket_arc-deg2rad(3);
ok(i)=far_from_seam && ring_gap(zz(j))>=need_z;
end
for b=2:nb
target=(b-1)/nb*loop_len;                      %arc length from the feed along the litz
[~,i]=min(abs(s(cand(ok))-target)+1e9*(~ok(ok)));
okc=cand(ok); j=okc(i);
g.break_s(b)=s(j); g.break_theta(b)=mod(th(j),2*pi); g.break_z(b)=zz(j);
g.break_xyz(:,b)=litz(:,j);
end
%chunk lengths (loop: feed -> breaks -> end of litz -> strip -> feed)
edges=[g.break_s(2:end) g.litz_length];
g.chunk_length=diff([0 edges]); g.chunk_length(end)=g.chunk_length(end)+g.strip_length;
g.pocket_need_z=need_z;

%rotate the whole coil so the seam and the return sit at seam_deg
if o.seam_deg~=0
c=cosd(o.seam_deg); sn=sind(o.seam_deg); Rz=[c -sn 0; sn c 0; 0 0 1];
g.litz=Rz*g.litz; g.strip_centre=Rz*g.strip_centre; g.break_xyz=Rz*g.break_xyz;
for i=1:numel(g.strip_filaments), g.strip_filaments{i}=Rz*g.strip_filaments{i}; end
r0=deg2rad(o.seam_deg);
g.theta=g.theta+r0; g.break_theta=mod(g.break_theta+r0,2*pi);
g.feed_theta=g.feed_theta+r0; g.strip_theta_end=g.strip_theta_end+r0;
end
g.seam_theta=deg2rad(o.seam_deg);

function add(t,zv,pt)
th=[th t(:)']; zz=[zz zv(:)']; part=[part pt*ones(1,numel(t))];
end
end
