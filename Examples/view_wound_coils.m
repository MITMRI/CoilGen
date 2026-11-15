%This script shows the wound gradient coils of halbach_wound_gradient_set.m in
%3D with CoilGen's plotting functions, plus all three nested in the magnet
%frame. Outside the MATLAB desktop the figures are written as PNG files.
%
%CoilGen frame (as in the plots): B0 = z, bore along x.
%Magnet frame: x_magnet = z, y_magnet = -y, z_magnet = x.

if ~exist('skip_cd','var')
if ispc
cd('..\');
else
cd('../');
end
end
addpath(fullfile(pwd,'plotting'));

result_folder=fullfile(pwd,'Wound_gradient_coils','halbach_320mm_bore');
figure_folder=fullfile(result_folder,'figures');
if ~exist(figure_folder,'dir'), mkdir(figure_folder); end
headless=~usejava('desktop');
if headless, set(0,'DefaultFigureVisible','off'); end

coil_names={'Z','Y','X'};
coil_titles={'Z gradient (bore axis)','Y gradient','X gradient (B0 axis)'};
colors=[0.80 0.25 0.15; 0.15 0.50 0.25; 0.15 0.35 0.75];

close all;
paths=cell(1,3);
for coil_ind=1:numel(coil_names)
loaded=load(fullfile(result_folder,'gauge_study',[coil_names{coil_ind} '_AWG16.mat']));
paths{coil_ind}=loaded.res.wire_path;
coil_layouts=loaded.res.coil_layouts;
first_figure=numel(findall(0,'Type','figure'));
plot_coil_track_with_resulting_bfield(coil_layouts,1,coil_titles{coil_ind});
plot_2D_contours_with_sf(coil_layouts,1,coil_titles{coil_ind});
plot_various_error_metrics(coil_layouts,1,coil_titles{coil_ind});
plot_resulting_gradient(coil_layouts,1,coil_titles{coil_ind});
if headless
figs=findall(0,'Type','figure');
figs=sort([figs.Number]);
for k=figs(first_figure+1:end)
exportgraphics(figure(k),fullfile(figure_folder,sprintf('%s_view_%d.png',coil_names{coil_ind},k-first_figure)),'Resolution',130);
end
end
end

%each coil: wire in 3D and the winding pattern unrolled, magnet frame
shapes=figure('Color','w','Position',[50 50 1500 1300]);
layout=tiledlayout(shapes,3,2,'TileSpacing','compact','Padding','compact');
for coil_ind=1:3
v=paths{coil_ind};
p=[v(3,:);-v(2,:);v(1,:)];
nexttile(layout);
plot3(p(3,:)*1e3,p(1,:)*1e3,p(2,:)*1e3,'Color',colors(coil_ind,:),'LineWidth',0.4);
axis equal tight; grid on; view(35,22);
xlabel('z, bore [mm]'); ylabel('x, B0 [mm]'); zlabel('y [mm]');
title(coil_titles{coil_ind});
nexttile(layout);
arc=atan2(p(2,:),p(1,:))*hypot(p(1,1),p(2,1))*1e3;
z=p(3,:)*1e3;
arc_break=[abs(diff(arc))>100 false];
arc(arc_break)=NaN;
plot(z,arc,'Color',colors(coil_ind,:),'LineWidth',0.4);
axis equal tight; grid on;
xlabel('z, bore [mm]'); ylabel('arc around the former [mm]');
title([coil_names{coil_ind} ' unrolled']);
end
if headless, exportgraphics(shapes,fullfile(figure_folder,'coil_shapes.png'),'Resolution',120); end

%all three coils nested, one quadrant cut away to see inside
nested=figure('Color','w','Position',[50 50 1500 650]);
layout=tiledlayout(nested,1,2,'TileSpacing','compact','Padding','compact');
for view_ind=1:2
nexttile(layout); hold on;
for coil_ind=1:3
v=paths{coil_ind};
p=[v(3,:);-v(2,:);v(1,:)]; %magnet frame
if view_ind==2
p(:,p(1,:)<0&p(2,:)>0)=NaN; %cut away one quadrant
end
plot3(p(3,:)*1e3,p(1,:)*1e3,p(2,:)*1e3,'Color',colors(coil_ind,:),'LineWidth',0.5);
end
[sx,sy,sz]=sphere(40);
surf(150*sz,150*sx,150*sy,'FaceColor',[0.95 0.75 0.2],'EdgeColor','none','FaceAlpha',0.6);
axis equal tight; grid on; view(-40,22);
xlabel('z, bore [mm]'); ylabel('x, B0 [mm]'); zlabel('y [mm]');
end
title(layout,'Z (red) inside Y (green) inside X (blue), with the 30 cm imaging sphere');
if headless, exportgraphics(nested,fullfile(figure_folder,'nested_coils.png'),'Resolution',120); end
