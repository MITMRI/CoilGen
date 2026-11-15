%Close-up of one crossing of the Z coil in designs/12cm: where the wire leaves
%the inside of a lobe it passes under the turns in an open slot, touching them;
%the turns over it are 0.4 mm proud (stack_crossings).
%Left: 3D, wire at its real diameter, former surface see-through.
%Right: section along the lead (plane theta = 0): the turns cut across, the
%lead running under them, and the former wall.

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
folder=fullfile(pwd,'designs','12cm');
headless=~usejava('desktop');
if headless, set(0,'DefaultFigureVisible','off'); end
T=readtable(fullfile(folder,'wire_paths','Z_wire_path.csv'));
P=[T.x_mm T.y_mm T.z_mm]; th=atan2d(T.y_mm,T.x_mm); cr=T.crossing==1;
od=1.732; R_in=66; R_out=70;

%the longest crossing
runs=[find(diff([0;cr])==1) find(diff([cr;0])==-1)];
[~,k]=max(runs(:,2)-runs(:,1)); lead=runs(k,1):runs(k,2);
zc=median(P(lead,3)); win_z=zc+[-9 9]; win_th=[-9 9]/R_out*180/pi;  %+-9 mm

in=P(:,3)>win_z(1)-3 & P(:,3)<win_z(2)+3 & th>win_th(1)-3 & th<win_th(2)+3;
f=figure('Color','w','Position',[50 50 1700 700]);
t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); hold(ax,'on');
%wire tubes, piece by piece
pieces=split_runs(find(in));
for i=1:numel(pieces)
idx=pieces{i};
col=[0.85 0.45 0.15]; 
draw_tube(ax,P(idx,:),od/2,col,cr(idx));
end
%former outer surface (see-through) and inner surface
[TH,Z]=meshgrid(linspace(win_th(1),win_th(2),30),linspace(win_z(1),win_z(2),30));
surf(ax,R_out*cosd(TH),R_out*sind(TH),Z,'FaceColor',[0.6 0.7 0.9],'FaceAlpha',0.25,'EdgeColor','none');
surf(ax,R_in*cosd(TH),R_in*sind(TH),Z,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none');
axis(ax,'equal'); view(ax,[115 30]); camlight(ax,'headlight'); lighting(ax,'gouraud');
xlim(ax,[R_in-1 R_out+1]); ylim(ax,R_out*sind(win_th)); zlim(ax,win_z);
xlabel(ax,'x (radial here) [mm]'); ylabel(ax,'y [mm]'); zlabel(ax,'z, bore [mm]');
title(ax,{'Z coil crossing, 3D (18 x 18 mm)','orange: turns in their grooves (lifted 0.4 mm over the lead), red: lead under them','blue: outer former surface, grey: inside of the former'});

%section in the plane theta = 0 (r-z)
ax2=nexttile(t); hold(ax2,'on');
fill(ax2,[win_z(1) win_z(2) win_z(2) win_z(1)],[R_in R_in R_out R_out],[0.85 0.88 0.95],'EdgeColor','none');
r=hypot(P(:,1),P(:,2));
%turns cut by the plane
s=find(sign(th(1:end-1))~=sign(th(2:end)) & abs(th(1:end-1))<20 & ~cr(1:end-1));
turn_r=[];
for i=s'
a=th(i)/(th(i)-th(i+1)); q=P(i,:)+a*(P(i+1,:)-P(i,:));
if q(3)<win_z(1)||q(3)>win_z(2), continue; end
turn_r(end+1)=hypot(q(1),q(2));
rectangle(ax2,'Position',[q(3)-od/2 hypot(q(1),q(2))-od/2 od od],'Curvature',1,'FaceColor',[0.85 0.45 0.15],'EdgeColor','k');
%groove
fill(ax2,q(3)+[-1 1 1 -1]*od/2,[R_out R_out R_out-od R_out-od]+0.001,[1 1 1],'EdgeColor','none','FaceAlpha',0);
end
%lead along the plane
L=lead(abs(th(lead))<2 & P(lead,3)>win_z(1)-1 & P(lead,3)<win_z(2)+1);
plot(ax2,P(L,3),r(L)+od/2,'k',P(L,3),r(L)-od/2,'k');
fill(ax2,[P(L,3);flipud(P(L,3))],[r(L)+od/2;flipud(r(L)-od/2)],[0.8 0.15 0.1],'EdgeColor','none');
yline(ax2,R_out,'b-','outer surface','LabelHorizontalAlignment','left');
yline(ax2,R_in,'-','inside of former','Color',[0.4 0.4 0.4],'LabelHorizontalAlignment','left');
axis(ax2,'equal'); xlim(ax2,win_z); ylim(ax2,[R_in-1.5 R_out+1.5]);
xlabel(ax2,'z, bore [mm]'); ylabel(ax2,'radius [mm]');
title(ax2,{'Cut along the lead (side view)', ...
    sprintf('turns over the lead: top %.2f mm proud of the surface; lead centre %.2f mm below the surface', max(turn_r)+od/2-R_out,R_out-median(r(L))), ...
    sprintf('gap between lead and turns %.2f mm (touching), plastic under the lead %.2f mm',max(turn_r)-median(r(L))-od,median(r(L))-od/2-R_in)});
if headless, exportgraphics(f,fullfile(folder,'figures','Z_crossing_closeup.png'),'Resolution',140); end

function c=split_runs(idx)
br=[0;find(diff(idx)>1);numel(idx)];
c=arrayfun(@(i) idx(br(i)+1:br(i+1)),1:numel(br)-1,'UniformOutput',false);
end

function draw_tube(ax,p,rad,col,is_cross)
if size(p,1)<2, return; end
n=16; a=linspace(0,2*pi,n);
t=gradient(p')'; t=t./vecnorm(t,2,2);
radial=[p(:,1:2) zeros(size(p,1),1)]; radial=radial./vecnorm(radial,2,2);
u=cross(t,radial,2); u=u./vecnorm(u,2,2); w=cross(t,u,2);
X=p(:,1)+rad*(u(:,1)*cos(a)+w(:,1)*sin(a));
Y=p(:,2)+rad*(u(:,2)*cos(a)+w(:,2)*sin(a));
Z=p(:,3)+rad*(u(:,3)*cos(a)+w(:,3)*sin(a));
C=repmat(reshape(col,1,1,3),size(X,1),n,1);
C(is_cross,:,1)=0.8; C(is_cross,:,2)=0.15; C(is_cross,:,3)=0.1;
surf(ax,X,Y,Z,C,'EdgeColor','none');
end
