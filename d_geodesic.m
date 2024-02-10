function dist=d_geodesic(u0,v0,u1,v1,df_du_total,df_dv_total)
df_du_1=df_du_total(:);
df_dv_1=df_dv_total(:);
nx=size(df_du_total,1)
ny=size(df_du_total,2)
nd=20;% number of segments
%Legendre-Gauss Quadrature Weights and Nodes
[absc,wght]=lgwt(nd,-1,1);
ulam=@(t) u0+(u1-u0)*t';
vlam=@(t) v0+(v1-v0)*t';
t=0.5*(absc+1);
ult=ulam(t);% u(lambda)
vlt=vlam(t);% v(lambda)
% interpolation u_lambda and v_lambda
df_du=interp2(reshape(u1,ny,nx),reshape(v1,ny,nx),reshape(df_du_1,ny,nx),ult(:),vlt(:));
df_dv=interp2(reshape(u1,ny,nx),reshape(v1,ny,nx),reshape(df_dv_1,ny,nx),ult(:),vlt(:));

% Computing distance of each segment and sum
g11=1+reshape(df_du,nx*ny,nd).^2;
g12=reshape(df_du,nx*ny,nd).*reshape(df_dv,nx*ny,nd);
g22=1+reshape(df_dv,nx*ny,nd).^2;

% Integration 
dist=sqrt((u1-u0).^2.*g11+2*(u1-u0).*(v1-v0).*g12+(v1-v0).^2.*g22)*0.5*wght;
dist=sum(dist,2);
end

