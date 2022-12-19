clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Square geometry with known exact solution (no penetration boundary conditions)
%    u   =  sin(2*pi*y);
%
%   +--------+
%   |        |
%   | Omega  |
%   |        |
%   +--------+
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Problem = struct(...
'Title'             ,  'Starting',  ...
'Subtitle'          ,  'sinus',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'id',       ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [1,1],      ...
'H_Max'             ,  1/8,        ...
'H_Min'             ,  1/8,        ...
'Reynolds'          ,  1,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Latex'             ,  true,       ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10],    ...
'Ref_Strategy'      ,  'href',    ...
'Mark_Param'        ,   .1 ); 


Convergence_rates = struct( ...
'uniform',      true,       ...
'p_values',      1,       ...
'iterations',    3);

syms x y
a = 1;
b = 2*pi;
% u(x,y) = sin (b*y); 
u(x,y) = exp (a*x) .* sin (b*y); 
main_make_exact_solution;

main_init;
gauss_n = gauss_n + 1
main_assemble_square_nmnn
main_static;
integrateNorms;
main_dump_iteration_results;
figure; lr.surf(u);          title('u'); colormap('turbo');

return

result_h     = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
result_dof   = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
result_uh_H1 = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
result_uh_L2 = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));

for iteration_p=1:numel(Convergence_rates.p_values)
  Problem.Polynomial_Degree = ones(1,2)*Convergence_rates.p_values(iteration_p);
  main_init;
  h_val_result = Problem.H_Max;
  dof_val_result = size(lr.cp, 2);
  for iteration_h=1:Convergence_rates.iterations

    main_init_iteration
    main_assemble_square_nmnn
    main_static;
    integrateNorms;
    main_dump_iteration_results

    % return

    result_h(    iteration_h, iteration_p) = h_val_result;
    result_dof(    iteration_h, iteration_p) = dof_val_result
    result_uh_H1(iteration_h, iteration_p) = sqrt(sum(velocity_error_H1_squared)/sum(u_H1_norm_squared))
    result_uh_L2(iteration_h, iteration_p) = sqrt(sum(velocity_error_L2_squared)/sum(u_L2_norm_squared))

    disp 'velocity (energy) convergence results'
    diff(log(result_uh_H1)) ./ diff(log(result_h))
    disp 'velocity convergence results'
    diff(log(result_uh_L2)) ./ diff(log(result_h))

% mrk_fn20_u_H1_norm_squared = find( function_velocity_error_H1_squared >= prctile(function_velocity_error_H1_squared, [Problem.Mark_Param]));
mrk_basis = find( function_velocity_error_H1_squared >= Problem.Mark_Param*max(function_velocity_error_H1_squared));

figure; lr.surf(u);          title('u'); colormap('turbo');
figure; lr.surf(lr.p(:,1)); axis off;
title('polynomial degree'); 
% colorbar; 
clim([2,8]); 
% colormap('lines');
colormap(mycolormap);
pbaspect([1 1 1]); hold off;

  end
end
main_dump_final_results

exprt = [result_h result_dof result_uh_H1 result_uh_L2];
fid = fopen('sinus_href_mark20_p1.csv','w');
fprintf(fid, 'h\t\tDOF\t\tenergy\t\tL2\n');
fprintf(fid, '%.16f\t %d\t %.16f\t %.16f\n',exprt');
fclose(fid);

% figure;
% lr.plot('enumeration');