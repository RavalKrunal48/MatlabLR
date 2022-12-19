
fprintf('assembling system\n'); % add a linebreak since assembly proccedure prints progress
fprintf('\n'); 

t = cputime;
tic;
%%%%%   ASSEMBLE THE 'STIFFNESS' MATRIX A  %%%%%

N = size(lr.knots,1);
nel = size(lr.elements,1);

A = sparse(N,N);
b = zeros(N,1);

gauss_n
%%% pre-evaluate bezier functions
xg = GaussLegendre(gauss_n(1));
yg = GaussLegendre(gauss_n(2));

nel
fprintf('(  0%%)');
% for all elements
for el=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el/nel*100)); % print progress to screen

  el_du = lr.elements(el,3) - lr.elements(el,1);
  el_dv = lr.elements(el,4) - lr.elements(el,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  ug = (xg+1)/2.0*el_du + lr.elements(el,1);
  vg = (yg+1)/2.0*el_dv + lr.elements(el,2);

  ind = lr.support{el};
  sup = numel(ind);

  % initialize element matrices
  Ak  = zeros(sup);

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),
      N     = lr.computeBasis(ug(gauss_i),vg(gauss_j), 1);
      x     = lr.point(ug(gauss_i), vg(gauss_j), 1);
      Jt    = x(:,2:3);
      x     = x(:,1);
      dNdu  = N(2:3,:);
      N     = N(1,:);
      dNdx  = inv(Jt) * dNdu;

      detJw = det(Jt)*wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      %%% Test function syntax:
      Ak = Ak + dNdx'*dNdx* detJw;  % N_{i,j}^k  N_{i,j}^l  diffusion term,  for all (k,l)

      %%% right-hand side
      fVal    = Problem.Force(x(1), x(2));
      b(ind)  = b(ind) + N'*fVal * detJw;

    end
  end
  % end gauss points

  A(ind, ind)  = A(ind, ind) + Ak;
end
% end element loop

time_assemble        = cputime - t;
time_assemble_wall   = toc;

bodyForce = b;
b = zeros(size(b));

%%% set boundary conditions
disp 'setting boundary conditions'
% neumannBndryCond;

traction = b;
b = zeros(size(b));

edges     = [];
edgVal    = [];

disp('%%%%%%%%%%%%%%%% Extracting Neumann Boundary Edges %%%%%%%%%%%%%%%%')
% lr.plot('enumeration');
% edgesall  = lr.getEdge(0)

for iedge = 1:4
  edge_el{iedge} = lr.getEdge(iedge,'elements');
  edge{iedge}  = lr.getEdge(iedge)
end  

nmnn_ind = [];
iedge = 1
for iel = 1:length(edge_el{iedge})
  el = edge_el{iedge}(iel);
  el_dv = lr.elements(el,4) - lr.elements(el,2);
  [yg wyg] = GaussLegendre(gauss_n(2));
  vg = (yg+1)/2.0*el_dv + lr.elements(el,2);
  %% all functions on this element and on this edge
  % ind = intersect(lr.support{el},edge{iedge})
  ind = lr.support{el};
  nmnn_ind = [nmnn_ind; edge{iedge}];

  for gauss_j=1:gauss_n(2)
    N     = lr.computeBasis(0,vg(gauss_j), 1);
    x     = lr.point(0, vg(gauss_j), 1);
    Jt    = x(:,2:3);
    x     = x(:,1);
    dNdu  = N(2:3,:);
    N     = N(1,:);
    dNdx  = inv(Jt) * dNdu;
    detJw = det(Jt)*wyg(gauss_j)*el_dv / 2.0;
    
    h = Exact_solution.grad_u(x(1), x(2))*(-1);
    
    b(ind) = b(ind) + N'*h(1) * detJw;
  end
end
% return

iedge = 2
for iel = 1:length(edge_el{iedge})
  el = edge_el{iedge}(iel);
  el_dv = lr.elements(el,4) - lr.elements(el,2);
  [yg wyg] = GaussLegendre(gauss_n(2));
  vg = (yg+1)/2.0*el_dv + lr.elements(el,2);
  %% all functions on this element and on this edge
  ind = lr.support{el};
  nmnn_ind = [nmnn_ind; edge{iedge}];

  for gauss_j=1:gauss_n(2)
    N     = lr.computeBasis(1,vg(gauss_j), 1);
    x     = lr.point(1, vg(gauss_j), 1);
    Jt    = x(:,2:3);
    x     = x(:,1);
    dNdu  = N(2:3,:);
    N     = N(1,:);
    dNdx  = inv(Jt) * dNdu;
    detJw = det(Jt)*wyg(gauss_j)*el_dv / 2.0;
    
    h = Exact_solution.grad_u(x(1), x(2))*(1);
    
    b(ind) = b(ind) + N'*h(1) * detJw;
  end
end
b
% return

% iedge = 4
% for iel = 1:length(edge_el{iedge})
%   el = edge_el{iedge}(iel);
%   el_du = lr.elements(el,3) - lr.elements(el,1);
%   [xg wxg] = GaussLegendre(gauss_n(1));
%   ug = (xg+1)/2.0*el_du + lr.elements(el,1);
%   %% all functions on this element and on this edge
%   ind = lr.support{el};
%   nmnn_ind = [nmnn_ind ind];

%   for gauss_j=1:gauss_n(2)
%     N     = lr.computeBasis(ug(gauss_i),1, 1);
%     x     = lr.point(ug(gauss_i),1, 1);
%     Jt    = x(:,2:3);
%     x     = x(:,1);
%     dNdu  = N(2:3,:);
%     N     = N(1,:);
%     dNdx  = inv(Jt) * dNdu;
%     detJw = det(Jt)*wxg(gauss_i)* el_du/ 2.0;

  
%     h = Exact_solution.grad_u(x(1), x(2))*(1);
    
%     if lr.elements(el,1) >= 0.5    
%       b(ind) = b(ind) + N'*h(2) * detJw;
%     else
%       b(ind) = b(ind) + N'*h(1) * detJw;
%     end
%   end
% end

traction = traction + b
nmnn_ind = unique(nmnn_ind)
% 
disp('%%%%%%%%%%%%%%%% Extracting Dirichlet Boundary Edges %%%%%%%%%%%%%%%%')
% Dirichlet BC at edge-3
edgVal3 = zeros(size(edge{3}))
edge{3}
edge{4}
edgVal = zeros(size([edge{3}; edge{4}]))

% strip down DOFs appearing on multple edges (i.e. corners)
[edges i] = unique([edge{3}; edge{4}])
edgVal    = edgVal(i);

%%% put all boundary conditions into traction-vector
if numel(edges)>0
  traction = traction - A(:,edges)*edgVal;
end
traction
% return


%%% remove boundary DOFs from the system
A(:,edges)             = [];
A(edges,:)             = [];
traction(edges)        = [];
bodyForce(edges)       = [];

%%% always include these terms
dF = @(u) A;
F  = @(u) dF(u)*u - bodyForce;

%%% add boundary conditions 
if Problem.Static
  F = @(u) F(u) - traction;
end


% iedge = 3
% for iel = 1:length(edge_el{iedge})
%   el = edge_el{iedge}(iel);
%   el_du = lr.elements(el,3) - lr.elements(el,1);
%   [xg wxg] = GaussLegendre(gauss_n(1));
%   ug = (xg+1)/2.0*el_du + lr.elements(el,1);
%   %% all functions on this element and on this edge
%   ind = lr.support{el};

%   for gauss_j=1:gauss_n(2)
%     N     = lr.computeBasis(ug(gauss_i),0, 1);
%     x     = lr.point(ug(gauss_i),0, 1);
%     Jt    = x(:,2:3);
%     x     = x(:,1);
%     dNdu  = N(2:3,:);
%     N     = N(1,:);
%     dNdx  = inv(Jt) * dNdu;
%     detJw = det(Jt)*wxg(gauss_i)* el_du/ 2.0;
    
%     h = Exact_solution.grad_u(x(1), x(2))*(-1);
    
%     % b(ind) = b(ind) + N(ind)'*h(1) * detJw
%     b(ind) = b(ind) + N'*h(2) * detJw;
%   end
% end