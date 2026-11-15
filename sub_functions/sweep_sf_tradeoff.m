function res=sweep_sf_tradeoff(spec,result_file)
%Linearity versus power of the stream function of one cylindrical gradient
%coil, for a list of regularization weights.
%
%CoilGen sets up the mesh, the sensitivity matrix and the resistance matrix.
%The regularization here is the dissipated power x'*R*x (CoilGen's own
%Tikhonov term is |R*x|^2), solved for all weights at once with one SVD.
%Each solution is checked on points that the optimization did not use.
%
%spec: .shape .gradient_dir .length .radius .roi_radius .mu (list)
%      .num_circ .axial_step  .target_file (in target_fields)  .target_name
%result_file: .mat file; nothing is recomputed if it exists.
%
%CoilGen frame: B0 = z, bore along x.

if exist(result_file,'file'), res=load(result_file); res=res.res; return; end

num_axial=max(20,round(spec.length/spec.axial_step));
coil_out=CoilGen(...
    'field_shape_function',spec.shape,...
    'coil_mesh_file','create cylinder mesh',...
    'cylinder_mesh_parameter_list',[spec.length spec.radius spec.num_circ num_axial 0 1 0 pi/2],...
    'surface_is_cylinder_flag',true,...
    'target_field_definition_file',spec.target_file,...
    'target_field_definition_field_name',spec.target_name,...
    'levels',10,...
    'tikonov_reg_factor',1,...
    'skip_postprocessing',true,...
    'skip_inductance_calculation',true,...
    'skip_sweep',true,...
    'save_stl_flag',false,...
    'plot_flag',false,...
    'conductor_thickness',1,...                 %R in ohm*m (per 1/sigma), only relative
    'specific_conductivity_conductor',1,...
    'output_directory',tempdir);
addpath(fullfile(pwd,'sub_functions'));

part=coil_out.coil_parts(1);
mesh=part.coil_mesh;
num_v=size(mesh.vertices,2);
S=squeeze(part.sensitivity_matrix(3,:,:));       %Bz at targets per unit sf
R=part.resistance_matrix;
b=coil_out.target_field.b(3,:)';

%stream function on each cylinder end edge is one free value
[T,~]=boundary_reduction(mesh,num_v);
Sr=S*T; Rr=T'*R*T; Rr=(Rr+Rr')/2;
L=chol(Rr+eye(size(Rr))*1e-12*trace(Rr)/size(Rr,1),'lower');
A=Sr/L';                                          %y = L'*x
[U,sv,V]=svd(A,'econ'); sv=diag(sv);
Ub=U'*b;

%validation points: two spheres that the optimization did not use
val=validation_points(spec.roi_radius);
faces=mesh.faces; verts=mesh.vertices;
[quad_pts,quad_w,quad_face]=face_quadrature(verts,faces,3);
face_area=accumarray(quad_face,quad_w);

num_mu=numel(spec.mu);
res.spec=spec; res.mesh=struct('vertices',verts,'faces',faces);
res.mu=spec.mu(:)';
f={'gradient_T_per_m','nonlin_max_pct','nonlin_rms_pct','nonlin_max_pct_inner','pos_err_mm', ...
   'power_rel','K_max','K_p99','K_int','sf_range','fit_rms_target_pct'};
for i=1:numel(f), res.(f{i})=nan(1,num_mu); end
res.sf=zeros(num_v,num_mu);
Y=V*((sv./(sv.^2+res.mu)).*Ub);                  %all weights at once
X=T*(L'\Y);
Kall=zeros(3,size(faces,2),num_mu);
for k=1:num_mu, Kall(:,:,k)=face_current_density(verts,faces,mesh.n,X(:,k)); end
Bz=sheet_field_z(quad_pts,quad_w,quad_face,Kall,val.points);   %Nval x num_mu
outer=val.radius==spec.roi_radius; inner=~outer;
along=val.points*spec.gradient_dir(:);
for k=1:num_mu
fit=[along ones(numel(along),1)]\Bz(:,k);
g=fit(1);
x=X(:,k)/g; Kf=Kall(:,:,k)/g;
resid=Bz(:,k)/g-along-fit(2)/g;
Kn=vecnorm(Kf,2,1);
res.sf(:,k)=x;
res.gradient_T_per_m(k)=g;
res.nonlin_max_pct(k)=max(abs(resid(outer)))/spec.roi_radius*100;
res.nonlin_rms_pct(k)=rms(resid(outer))/spec.roi_radius*100;
res.nonlin_max_pct_inner(k)=max(abs(resid(inner)))/max(val.radius(inner))*100;
res.pos_err_mm(k)=max(abs(resid(outer)))*1e3;
res.power_rel(k)=x'*R*x;
res.K_max(k)=max(Kn);
res.K_p99(k)=prctile(Kn,99);
res.K_int(k)=sum(Kn(:).*face_area(:));
res.sf_range(k)=max(x)-min(x);
bt=S*x; res.fit_rms_target_pct(k)=rms(bt-b)/spec.roi_radius*100;
end
save(result_file,'res','-v7.3');
end

function [T,keep]=boundary_reduction(mesh,num_v)
%x = T*z: every boundary loop shares one value
bnd=cellfun(@(b) unique(b(:)),mesh.boundary,'UniformOutput',false);
is_b=false(num_v,1); is_b(cat(1,bnd{:}))=true;
keep=find(~is_b);
cols=numel(keep)+numel(bnd);
rows=[keep;cat(1,bnd{:})];
c=[(1:numel(keep))';repelem(numel(keep)+(1:numel(bnd))',cellfun(@numel,bnd))];
T=sparse(rows,c,1,num_v,cols);
T=full(T);
end

function val=validation_points(r)
n1=3000; n2=1500;
val.points=[fib(n1)*r; fib(n2)*r*5/6];
val.radius=[ones(n1,1)*r; ones(n2,1)*r*5/6];
end

function p=fib(n)
k=(0:n-1)'+0.5; polar=acos(1-2*k/n); az=pi*(1+sqrt(5))*k;
p=[cos(az).*sin(polar) sin(az).*sin(polar) cos(polar)];
end

function K=face_current_density(verts,faces,vn,x)
%surface current density (A/m) of a piecewise linear stream function: K = grad(sf) x n
p1=verts(:,faces(1,:)); p2=verts(:,faces(2,:)); p3=verts(:,faces(3,:));
e1=p2-p1; e2=p3-p1;
n=cross(e1,e2); a2=vecnorm(n); n=n./a2;
%outward normals as in the mesh
vn_f=vn(:,faces(1,:))+vn(:,faces(2,:))+vn(:,faces(3,:));
s=sign(sum(n.*vn_f)); n=n.*s;
s1=x(faces(1,:))'; s2=x(faces(2,:))'; s3=x(faces(3,:))';
%gradient of the linear interpolant
g=(cross(n,p3-p2).*s1+cross(n,p1-p3).*s2+cross(n,p2-p1).*s3)./a2;
K=cross(g,n);
end

function [pts,w,fid]=face_quadrature(verts,faces,m)
%m x m subdivision of each triangle, one point per sub-triangle
[a,b]=ndgrid(0:m-1); keep=a+b<=m-1; a=a(keep); b=b(keep);
up=[(a+1/3) (b+1/3)]/m;
[a2,b2]=ndgrid(0:m-2); keep=a2+b2<=m-2; a2=a2(keep); b2=b2(keep);
down=[(a2+2/3) (b2+2/3)]/m;
bary=[up;down];
p1=verts(:,faces(1,:))'; p2=verts(:,faces(2,:))'; p3=verts(:,faces(3,:))';
area=vecnorm(cross(p2-p1,p3-p1,2),2,2)/2;
nf=size(faces,2); nq=size(bary,1);
pts=zeros(nf*nq,3); w=zeros(nf*nq,1); fid=zeros(nf*nq,1);
for q=1:nq
idx=(1:nf)+(q-1)*nf;
pts(idx,:)=p1+bary(q,1)*(p2-p1)+bary(q,2)*(p3-p1);
w(idx)=area/nq; fid(idx)=1:nf;
end
end

function Bz=sheet_field_z(pts,w,fid,K,P)
%Bz of the current sheet for several current distributions K (3 x faces x n)
mu0=4*pi*1e-7;
n=size(K,3);
Jx=reshape(K(1,fid,:),[],n).*w; Jy=reshape(K(2,fid,:),[],n).*w;
Bz=zeros(size(P,1),n);
chunk=max(1,floor(2e6/size(pts,1)));
for i0=1:chunk:size(P,1)
ii=i0:min(i0+chunk-1,size(P,1));
rx=P(ii,1)-pts(:,1)'; ry=P(ii,2)-pts(:,2)'; rz=P(ii,3)-pts(:,3)';
r3=(rx.^2+ry.^2+rz.^2).^1.5;
Bz(ii,:)=(ry./r3)*Jx-(rx./r3)*Jy;
end
Bz=Bz*mu0/(4*pi);
end
