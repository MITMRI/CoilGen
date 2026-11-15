%Pictures of the 12 cm Halbach gradient set in designs/12cm (written by
%export_halbach_12cm_designs.m) to judge
%how dense the winding is: CoilGen's own plots, each coil unrolled with a 1:1
%zoom of its densest area (wire drawn at its real diameter), the coils nested,
%and a winding summary (turns, crossovers, tightest bends).
%Outside the MATLAB desktop the figures are written as PNG files.
%
%CoilGen frame: B0 = z, bore along x.  Magnet frame: x_magnet = z_cg (B0),
%y_magnet = -y_cg, z_magnet = x_cg (bore).

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'plotting'),fullfile(pwd,'sub_functions'));

folder=fullfile(pwd,'designs','12cm','results');   %from export_halbach_12cm_designs.m
figure_folder=fullfile(pwd,'designs','12cm','figures');
if ~exist(figure_folder,'dir'), mkdir(figure_folder); end
headless=~usejava('desktop');
if headless, set(0,'DefaultFigureVisible','off'); end

choice={'Z','Z_coil','Z gradient (bore axis), 330 mm'; ...
        'Y','Y_coil','Y gradient, 324 mm'; ...
        'X','X_coil','X gradient (B0 axis), 318 mm'};
colors=[0.80 0.25 0.15; 0.15 0.50 0.25; 0.15 0.35 0.75];
zoom_size=[50 30];   %mm, bore x arc

close all;
res=cell(1,3);
fid=fopen(fullfile(figure_folder,'winding_summary.txt'),'w');
for out=[1 fid]
fprintf(out,'%-22s %5s %6s %7s %10s %11s %12s\n','coil','turns','wire','former','crossings','tightest','bends <5 mm');
fprintf(out,'%-22s %5s %6s %7s %10s %11s %12s\n','','','[m]','r [mm]','','bend [mm]','');
end
for c=1:3
r=load(fullfile(folder,[choice{c,2} '.mat'])); r=r.res; res{c}=r;
%CoilGen's plots
first=numel(findall(0,'Type','figure'));
plot_2D_contours_with_sf(r.coil_layouts,1,choice{c,3});
plot_3D_contours_with_sf(r.coil_layouts,1,choice{c,3});
if headless
figs=findall(0,'Type','figure'); figs=sort([figs.Number]);
names={'contours_2D','contours_3D'};
for k=1:numel(figs)-first
exportgraphics(figure(figs(first+k)),fullfile(figure_folder,sprintf('%s_%s.png',choice{c,1},names{k})),'Resolution',130);
end
end
%winding numbers
wp=r.wire_path;
shifted=is_crossing(r);
crossings=sum(diff([0 shifted])==1);
s=[0 cumsum(vecnorm(diff(wp,1,2)))]; [s,iu]=unique(s);
q=interp1(s,wp(:,iu)',0:0.001:s(end))';       %1 mm steps
q=movmean(q,5,2);                               %below the wire diameter
a=q(:,1:end-4)-q(:,3:end-2); b=q(:,5:end)-q(:,3:end-2); cc=q(:,5:end)-q(:,1:end-4);
bend=vecnorm(a).*vecnorm(b).*vecnorm(cc)./(2*vecnorm(cross(a,b)));  %circumradius over 4 mm
tight=bend<0.005; tight_runs=find(diff([0 tight])==1);
num_tight=numel(tight_runs);
for out=[1 fid]
fprintf(out,'%-22s %5d %6.1f %7.2f %10d %11.1f %12d\n',choice{c,2},r.spec.levels,r.wire_length, ...
    r.spec.radius*1e3,crossings,prctile(bend,0.5)*1e3,num_tight);
end
end
fprintf(fid,['\nturns = contour levels; crossings = places where the wire passes under other turns\n' ...
    '(into the former); tightest bend = 0.5 percentile of the bend radius over 4 mm of wire;\nbends <5 mm = places bent tighter than 5 mm radius.\n']);
fclose(fid);

%unrolled, true scale, with a 1:1 zoom of the densest area
for c=1:3
r=res{c}; v=r.wire_path; od=r.wire.outer_diameter*1e3;
p=[v(3,:);-v(2,:);v(1,:)]*1e3;                 %magnet frame, mm
z=p(3,:); radius=r.spec.radius*1e3;
arc=atan2(p(2,:),p(1,:))*radius;
brk=[abs(diff(arc))>100 false]; arc(brk)=NaN;
%densest area: most wire inside the zoom window
pts=[z(~isnan(arc))' arc(~isnan(arc))'];
[gz,ga]=meshgrid(min(z):5:max(z)-zoom_size(1),-pi*radius:5:pi*radius-zoom_size(2));
count=zeros(size(gz));
for k=1:numel(gz)
count(k)=sum(pts(:,1)>gz(k)&pts(:,1)<gz(k)+zoom_size(1)&pts(:,2)>ga(k)&pts(:,2)<ga(k)+zoom_size(2));
end
[~,k]=max(count(:)); win=[gz(k) ga(k)];

f=figure('Color','w','Units','centimeters','Position',[2 2 34 17]);
t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); hold(ax,'on');
plot(ax,z,arc,'Color',colors(c,:),'LineWidth',0.5);
sh=is_crossing(r); zs=z; zs(~sh)=NaN;
plot(ax,zs,arc,'k','LineWidth',1.2);
rectangle(ax,'Position',[win zoom_size],'EdgeColor','k','LineWidth',1);
plot(ax,z([1 end]),atan2([p(2,1) p(2,end)],[p(1,1) p(1,end)])*radius,'ko','MarkerFaceColor','y','MarkerSize',8);
axis(ax,'equal'); axis(ax,'tight'); grid(ax,'on');
xlabel(ax,'z along bore [mm]'); ylabel(ax,'arc around former from +x (B0) [mm]');
title(ax,sprintf('%s unrolled: %d turns, %.1f m of 14 AWG (black: crossings, yellow: terminals)',choice{c,3},r.spec.levels,r.wire_length));
%zoom, wire at real diameter
ax2=nexttile(t); hold(ax2,'on');
in=z>win(1)-5&z<win(1)+zoom_size(1)+5&arc>win(2)-5&arc<win(2)+zoom_size(2)+5;
zz=z; aa=arc; zz(~in)=NaN; aa(~in)=NaN;
h=plot(ax2,zz,aa,'Color',colors(c,:));
axis(ax2,'equal'); xlim(ax2,win(1)+[0 zoom_size(1)]); ylim(ax2,win(2)+[0 zoom_size(2)]);
grid(ax2,'on'); box(ax2,'on');
drawnow;
pos=getpixelposition(ax2); pts_per_mm=pos(3)/zoom_size(1)*72/get(0,'ScreenPixelsPerInch');
h.LineWidth=od*pts_per_mm;
plot(ax2,zz,aa,'k-','LineWidth',0.3);
xlabel(ax2,'z [mm]'); ylabel(ax2,'arc [mm]');
title(ax2,sprintf('densest %dx%d mm, wire %.2f mm + 0.8 mm plastic',zoom_size,od));
if headless, exportgraphics(f,fullfile(figure_folder,sprintf('%s_unrolled.png',choice{c,1})),'Resolution',150); end
end

%all three nested, one quadrant cut away, with the 120 mm sphere
nested=figure('Color','w','Position',[50 50 1500 650]);
t=tiledlayout(nested,1,2,'TileSpacing','compact','Padding','compact');
for view_ind=1:2
nexttile(t); hold on;
for c=1:3
v=res{c}.wire_path; p=[v(3,:);-v(2,:);v(1,:)]*1e3;
if view_ind==2, p(:,p(1,:)<0&p(2,:)>0)=NaN; end
plot3(p(3,:),p(1,:),p(2,:),'Color',colors(c,:),'LineWidth',0.5);
end
[sx,sy,sz]=sphere(40);
surf(60*sz,60*sx,60*sy,'FaceColor',[0.95 0.75 0.2],'EdgeColor','none','FaceAlpha',0.6);
axis equal tight; grid on; view(-40,22);
xlabel('z, bore [mm]'); ylabel('x, B0 [mm]'); zlabel('y [mm]');
end
title(t,'Z (red) inside Y (green) inside X (blue), with the 120 mm imaging sphere');
if headless, exportgraphics(nested,fullfile(figure_folder,'nested_coils.png'),'Resolution',120); end

function shifted=is_crossing(r)
%wire parts that go into the former to pass under other turns
shifted=vecnorm(r.wire_path(2:3,:))<r.spec.radius-0.0005;
end
