/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Create uniform tetrahedral mesh
    */

#include <vtkSmartPointer.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>




int main(int argc, char **argv)
{
    std::cout << "tetmesh_uniform. Version 0.9"<< std::endl;

    const double h0 = 0.1; //distance between points in the initial distribution
    const double dptol = .001; //point displacement tolerance
    const double ttol = .1; //tolerance for retriangulation
    const double Fscale = 1.2; //scale factor to make sure forces are repulsive
    const double deltat = .2; //point update step
    const double geps = .001 * h0; //threshold to identify interior points
//    const double deps = sqrt(eps) * h0; //step for numerical gradient


    return 0;
}

/*
 function make_mesh
fd=inline('-0.3+abs(0.7-sqrt(sum(p.^2,2)))');
[p,t] = distmesh2d(fd,0.1,[-1,-1;1,1],[])


function [p,t]=distmesh2d(fd,h0,bbox,pfix,varargin)
dptol=.001; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0; deps=sqrt(eps)*h0;
% 1. Create initial distribution in bounding box (equilateral triangles)
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;

%fh = @huniform;

% Shift even rows
p=[x(:),y(:)];

% List of node coordinates
% 2. Remove points outside the region, apply the rejection method
p=p(feval(fd,p,varargin{:})<geps,:);

% Keep only d<0 points
%r0=1./feval(fh,p,varargin{:}).^2;
r0 = 1;

% Probability to keep point
%p=[pfix; p(rand(size(p,1),1)<r0./max(r0),:)];
p=[pfix; p];

% Rejection method
N=size(p,1);

% Number of points N
pold=inf;

while 1
    % 3. Retriangulation by the Delaunay algorithm
    if max(sqrt(sum((p-pold).^2,2))/h0)>ttol
        pold=p;
        t=delaunayn(p);
        pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
        t=t(feval(fd,pmid,varargin{:})<-geps,:);
        
        % 4. Describe each bar by a unique pair of nodes
        bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
        bars=unique(sort(bars,2),'rows');
        
        % 5. Graphical output of the current mesh
        trimesh(t,p(:,1),p(:,2),zeros(N,1))
        view(2),axis equal,axis off,drawnow
    end
    
    % For first iteration
    %     Any large movement?
    %     Save current positions
    %     List of triangles
    %     Compute centroids
    %     Keep interior triangles
    % Interior bars duplicated
    % Bars as node pairsfunction d=dcircle(p,xc,yc,r)
    %d=sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;
    d = fd(p);
    
    % 6. Move mesh points based on bar lengths L and forces F
    barvec=p(bars(:,1),:)-p(bars(:,2),:);
    % List of bar vectors
    L=sqrt(sum(barvec.^2,2));
    % L = Bar lengths
    %hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
    %L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));
    L0=Fscale*sqrt(sum(L.^2)/numel(L));
    % L0 = Desired lengths
    
    F=max(L0-L,0);
    % Bar forces (scalars)
    Fvec=F./L*[1,1].*barvec;
    
    % Bar forces (x,y components)
    Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:size(pfix,1),:)=0;
    
    % Force = 0 at fixed points
    p=p+deltat*Ftot;
    
    % Update node positions
    % 7. Bring outside points back to the boundary
    d=feval(fd,p,varargin{:}); ix=d>0;
    % Find points outside (d>0)
    dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
    dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; % gradient
    %p(ix,:)=p(ix,:)-[d(ix).*dgradx,d(ix).*dgrady];
    dgrad2=dgradx.^2+dgrady.^2;
    p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];
    
    % Project back to boundary
    % 8. Termination criterion: All interior nodes move less than dptol (scaled)
    if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
end

%function h=huniform(p,varargin)
%h=ones(size(p,1),1);


function d=dcircle(p,xc,yc,r)
d=sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;*/
