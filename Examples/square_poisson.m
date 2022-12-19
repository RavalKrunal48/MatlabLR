clc; clear;
p = [2,2];
controlpoints = [0,0; .5,  0; 1,0;
                 0,.5; .5,.5; 1,.5;
                 0,1;  .5,1; 1,1]';
knot1 = [0,0,0,1,1,1];
knot2 = [0,0,0,1,1,1];
lr = LRSplineSurface(p, knot1, knot2, controlpoints);

% perform tensor product refinement
% lr.refine();
% lr.refine();
basis_index = [];
% perform local refinement
for i=0
  element_index = lr.getElementContaining(0.25, 0.25);
  basis_index1   = lr.support{element_index};
  element_index = lr.getElementContaining(0.75, 0.25);
  basis_index2   = lr.support{element_index};
  element_index = lr.getElementContaining(0.25, 0.75);
  basis_index3   = lr.support{element_index};
  element_index = lr.getElementContaining(0.75, 0.75);
  basis_index4   = lr.support{element_index};
  basis_index   = [basis_index basis_index1 basis_index2 basis_index3 basis_index4] 
  lr.refine(basis_index, 'basis');

end
figure;
% lr.plot('basis');
% lr.localRaiseOrder(basis_index, 'basis')

% plot mesh for debugging purposes
% figure;
% lr.plot('enumerate');
print(lr)
lr.plot();

nel  = size(lr.elements,1); % number of elements
nfun = size(lr.knots,1);    % number of basis functions

A = zeros(nfun, nfun);      % poisson system matrix (stiffness matrix)
b = zeros(nfun, 1);         % right-hand side (load vector)

% 3-point gauss rule
ng = 3;                          % number of gauss points
xg = [-sqrt(3/5), 0, sqrt(3/5)]; % points
wg = [5,8,5]/9;                  % weights

for e=1:nel % for all elements
  
  % fetch element parametric size (u0,v0)x(u1,v1)
  u0 = lr.elements(e, 1);
  v0 = lr.elements(e, 2);
  u1 = lr.elements(e, 3);
  v1 = lr.elements(e, 4);

  I = lr.support{e}; % global index of all functions with support on this element
  n = numel(I);      % number of active functions on this element
  Ak = zeros(n);     % local stiffness matrix
  bk = zeros(n,1);   % local load vector

  m = 2;
  for i=1:ng % for all gauss points
    for j=1:ng
      u  = u0 + (xg(i)+1)/2 * (u1-u0); % parametric evaluation point (u,v)
      v  = v0 + (xg(j)+1)/2 * (v1-v0);
      N  = lr.computeBasis(u,v,1); % three rows of [N; dN/du; dN/dv]
      x  = lr.point(u,v,1);
      Jt = x(:,2:3); % transpose jacobian matrix: [dx/du, dy/du; dx/dv, dy/dv]
      x  = x(:,1);   % physical point where (u,v) maps to

      dNdu = N(2:3,:);       % parametric derivatives [dN/du, dN/dv]
      dNdx = inv(Jt) * dNdu; % physical derivatives   [dN/dx, dN/dy]

      for k=1:n % for all basis functions with support here
        for l=1:n
          Ak(k,l) = Ak(k,l) + dNdx(:,k)' * dNdx(:,l);
        end
%         bk(k) = bk(k) + N(1,k)*1;
      end
      %%% Optimization note (replace nested for-loops above):
%       Ak = dNdx' * dNdx
      bk = N(1,:)'*(-m.^2.*pi.^2.*sin(m*pi*u).*sin(m*pi*v)-m.^2.*pi.^2.*sin(m*pi*u).*sin(m*pi*v));

      % assemble contributions into global system matrix
      detJw = (u1-u0)*(v1-v0)/4*det(Jt)*wg(i)*wg(j); % weights and mapping contribution
      A(I,I) = A(I,I) + Ak*detJw;
      b(I)   = b(I)   + bk*detJw;
    end
  end
end

% add boundary conditions
e = lr.getEdge(0)
A(e,:) = 0;
A(:,e) = 0;
A(e,e) = eye(numel(e));
b(e)   = 0;


% solve system
u = A \ b;

% plot solution
figure;
lr.surf(u);
