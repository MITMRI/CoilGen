%Length study for the 3-axis gradient set of the 12 cm Halbach system:
%one layer of round wire in grooves on the outside of three printed formers,
%120 mm spherical imaging volume.
%   Z (bore) on the 66-70 mm former, Y on 71-75 mm, X (along B0) on 76-80 mm
%Frame of the magnet: B0 along x, bore along z. CoilGen frame: B0 = z_cg,
%bore = x_cg:  x_magnet = z_cg, y_magnet = -y_cg, z_magnet = x_cg.
%
%For each coil and length, the stream function is solved for a list of
%power weights (sweep_sf_tradeoff). Set SF_JOB to "Z,0.200" to run one job.

if ~exist('skip_cd','var')
if ispc, cd('..\'); else, cd('../'); end
end
addpath(fullfile(pwd,'sub_functions'));
maxNumCompThreads(str2double(getenv_default('SF_THREADS','3')));

result_folder=fullfile(pwd,'Wound_gradient_coils','halbach_12cm','sf_sweep');
roi_radius=0.060;
wire_od_14=0.0682*0.0254;
%       name shape dir      former outer radius
coils={'Z', 'x', [1;0;0], 0.070; ...
       'Y', 'y', [0;1;0], 0.075; ...
       'X', 'z', [0;0;1], 0.080};
lengths=[0.14:0.02:0.34 0.318 0.324 0.330];   %last three: the maximum length of each former

%target: 1 T/m on the surface of the 120 mm sphere
target_file=fullfile(pwd,'target_fields','halbach_12cm_targets.mat');
if ~exist(target_file,'file')
n=2500; k=(0:n-1)'+0.5; polar=acos(1-2*k/n); az=pi*(1+sqrt(5))*k;
p=roi_radius*[cos(az).*sin(polar) sin(az).*sin(polar) cos(polar)]';
targets=struct('coords',p,'gx',p(1,:),'gy',p(2,:),'gz',p(3,:));
save(target_file,'targets');
end

job=getenv('SF_JOB');
for c=1:size(coils,1)
for len=lengths
tag=sprintf('%s,%.3f',coils{c,1},len);
if ~isempty(job) && ~any(strcmp(tag,strsplit(job,';'))), continue; end
spec=struct('shape',coils{c,2},'gradient_dir',coils{c,3},'length',len, ...
    'radius',coils{c,4}-wire_od_14/2,'roi_radius',roi_radius,'mu',logspace(-14,-7,36), ...
    'num_circ',80,'axial_step',0.006,'target_file','halbach_12cm_targets.mat', ...
    'target_name',['g' coils{c,2}]);
file=fullfile(result_folder,sprintf('%s_L%03d.mat',coils{c,1},round(len*1e3)));
t0=tic; r=sweep_sf_tradeoff(spec,file);
fprintf('DONE %s %.0fs  nonlin %%: %s\n',tag,toc(t0),mat2str(round(r.nonlin_max_pct(1:4:end),2)));
end
end

function v=getenv_default(name,default)
v=getenv(name); if isempty(v), v=default; end
end
