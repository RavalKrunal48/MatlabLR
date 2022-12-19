clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Square geometry with known exact solution (no penetration boundary conditions)
%    u   =  sin(pi*x).*sin(pi*y);
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
'Static'            ,  false,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Latex'             ,  true,       ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10],    ...
'Ref_Strategy'      ,  'href',    ...
'Mark_Param'        ,   0.3 );

Convergence_rates = struct( ...
'uniform',      true,       ...
'p_values',      1,       ...
'iterations',    9);


mycolormap = [0 0.4470 0.7410
0.8500 0.3250 0.0980
0.9290 0.6940 0.1250
0.4940 0.1840 0.5560
0.4660 0.6740 0.1880
0.3010 0.7450 0.9330
0.6350 0.0780 0.1840];

syms x y
m = 2;
u(x,y) = sin(m*pi*x)*sin(m*pi*y); 
% u(x,y) = atan(100*(sqrt((x-1.25).^2 + (y+0.25).^2) - pi/3));
main_make_exact_solution;
% return
BC = cell(0);
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'value', 0)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'value', 0)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'value', 0)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'value', 0)];

if exist('Convergence_rates')
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
      main_assemble;
      main_static;
      integrateNorms;
      main_dump_iteration_results

      result_h(    iteration_h, iteration_p) = h_val_result;
      result_dof(    iteration_h, iteration_p) = dof_val_result
      result_uh_H1(iteration_h, iteration_p) = sqrt(sum(velocity_error_H1_squared)/sum(u_H1_norm_squared))
      result_uh_L2(iteration_h, iteration_p) = sqrt(sum(velocity_error_L2_squared)/sum(u_L2_norm_squared))

      disp 'velocity (energy) convergence results'
      diff(log(result_uh_H1)) ./ diff(log(result_h))
      disp 'velocity convergence results'
      diff(log(result_uh_L2)) ./ diff(log(result_h))

% mrk_fn20_u_H1_norm_squared = find( function_velocity_error_H1_squared >= prctile(function_velocity_error_H1_squared, [Problem.Mark_Param]));
mrk_basis = find( function_velocity_error_H1_squared >= (1-Problem.Mark_Param)*max(function_velocity_error_H1_squared));

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
else
  main_init;
  lr.localRaiseOrder(1:2:10, 'basis')
  lr.refine(1:10:100, 'basis')
  lr.localRaiseOrder(112, 'basis')
  lr.localRaiseOrder(lr.support{lr.getElementContaining(.499, .423)}, 'basis')
  lr.localRaiseOrder(lr.support{lr.getElementContaining(.499, .423)}, 'basis')
  gauss_n = gauss_n + 2
  main_assemble;

  if Problem.Static
    main_static;
    integrateNorms;
  end

  main_dump_iteration_results;
end
main_dump_final_results

exprt = [result_h result_dof result_uh_H1 result_uh_L2];
fid = fopen('sinus_href_mark20_p1.csv','w');
fprintf(fid, 'h\t\tDOF\t\tenergy\t\tL2\n');
fprintf(fid, '%.16f\t %d\t %.16f\t %.16f\n',exprt');
fclose(fid);

% figure;
% lr.plot('enumeration');