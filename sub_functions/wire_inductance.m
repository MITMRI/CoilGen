function L=wire_inductance(path_xyz,wire_radius,varargin)
%Low-frequency self inductance [H] of a round wire of radius wire_radius [m]
%with uniform current, following the polyline path_xyz (3xN, in m).
%
%Method: Dengler's curve integral (Adv. Electromagnetics 5(1), 2016),
%   L = mu0/(4 pi) * (double Neumann integral over |s-s'| > a/2) + mu0*len/(8 pi)
%where the last term holds the internal inductance of the wire. The path is
%resampled to segments of length h. Pairs further apart than near_factor*h use
%the midpoint rule, closer pairs use Gauss-Legendre quadrature, and a segment
%with itself or its neighbour along the wire is integrated in closed form.
%Options: 'segment_length' (m), 'near_factor', 'closed' (true for a closed loop).

p=inputParser;
addParameter(p,'segment_length',0.004,@isnumeric);
addParameter(p,'near_factor',4,@isnumeric);
addParameter(p,'closed',false,@islogical);
parse(p,varargin{:});
h_target=p.Results.segment_length;
near_factor=p.Results.near_factor;
closed=p.Results.closed;

%resample the path to segments of equal arc length
if closed && norm(path_xyz(:,1)-path_xyz(:,end))>1e-9
path_xyz=[path_xyz path_xyz(:,1)];
end
step=vecnorm(diff(path_xyz,1,2));
path_xyz=path_xyz(:,[true step>1e-9]);
s=[0 cumsum(vecnorm(diff(path_xyz,1,2)))];
wire_length=s(end);
num_seg=max(8,ceil(wire_length/h_target));
s_new=linspace(0,wire_length,num_seg+1);
nodes=interp1(s,path_xyz',s_new)';
A=nodes(:,1:end-1); B=nodes(:,2:end);
mid=(A+B)/2; dl=B-A;
h=mean(vecnorm(dl)); %nominal; chords of a curved path are slightly shorter
b=wire_radius/2;     %excluded |s-s'| of the Dengler integral

%gauss-legendre rule on [0,1]
[gx,gw]=gauss_rule(6);
[ga,gb]=ndgrid(gx,gx); wg=gw(:)*gw(:)'; ga=ga(:)'; gb=gb(:)'; wg=wg(:)';

near_dist=near_factor*h;
total=0;
chunk=max(1,floor(4e7/num_seg/8)); %about 40 MB per temporary
near_i=cell(0); near_j=cell(0);
for i0=1:chunk:num_seg
ii=i0:min(i0+chunk-1,num_seg);
dx=mid(1,ii)'-mid(1,:); dy=mid(2,ii)'-mid(2,:); dz=mid(3,ii)'-mid(3,:);
D=sqrt(dx.^2+dy.^2+dz.^2);
dots=dl(1,ii)'*dl(1,:)+dl(2,ii)'*dl(2,:)+dl(3,ii)'*dl(3,:);
is_near=D<near_dist;
kernel=dots./D; kernel(is_near)=0;
total=total+sum(kernel,'all');
[r,c]=find(is_near);
near_i{end+1}=ii(r)'; near_j{end+1}=c; %#ok<AGROW>
end
near_i=vertcat(near_i{:}); near_j=vertcat(near_j{:});

%pairs close to each other
sep=abs(near_i-near_j);
if closed, sep=min(sep,num_seg-sep); end
reg=sep>=2;
total=total+near_pairs(A,dl,near_i(reg),near_j(reg),ga,gb,wg);

%a segment with itself: the integral over |s-s'|>b of 1/|s-s'| on a straight piece
len=vecnorm(dl);
total=total+sum(2*(len.*log(len/b)-len+b));

%neighbours along the wire: nearly collinear, integral of 1/(u+v) on the square
%without the triangle u+v<b, scaled by the cosine of the corner angle
nb=sep==1;
ci=near_i(nb); cj=near_j(nb);
len_mean=(len(ci)+len(cj))/2;
cos_ang=sum(dl(:,ci).*dl(:,cj),1)./(len(ci).*len(cj));
total=total+sum(cos_ang.*(2*len_mean*log(2)-b));

L=1e-7*total+4*pi*1e-7*wire_length/(8*pi);
end

function total=near_pairs(A,dl,i_ind,j_ind,ga,gb,wg)
total=0;
block=2e5;
for k0=1:block:numel(i_ind)
k=k0:min(k0+block-1,numel(i_ind));
Ai=A(:,i_ind(k)); Aj=A(:,j_ind(k)); di=dl(:,i_ind(k)); dj=dl(:,j_ind(k));
ex=(Ai(1,:)'-Aj(1,:)')+di(1,:)'*ga-dj(1,:)'*gb;
ey=(Ai(2,:)'-Aj(2,:)')+di(2,:)'*ga-dj(2,:)'*gb;
ez=(Ai(3,:)'-Aj(3,:)')+di(3,:)'*ga-dj(3,:)'*gb;
inv_r=sum(wg./sqrt(ex.^2+ey.^2+ez.^2),2);
total=total+sum(sum(di.*dj,1)'.*inv_r);
end
end

function [x,w]=gauss_rule(n)
%Golub-Welsch, mapped to [0,1]
beta=(1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));
[x,ord]=sort((diag(D)+1)/2);
w=(V(1,ord)'.^2);
end
