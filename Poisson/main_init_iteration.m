if iteration_h > 1

switch (Problem.Ref_Strategy)
    case 'href'
        disp('Knot insertion')
        lr.refine(mrk_basis, 'basis');
        h_val_result = h_val_result / 2;

    case 'kref'
        if mod(iteration_h,2)
            disp('Knot insertion')
            lr.refine(mrk_basis, 'basis');
            h_val_result = h_val_result / 2;
        else 
            disp('Degree Elevation');
            lr.localRaiseOrder(mrk_basis, 'basis');
            gauss_n = gauss_n+1;
            h_val_result = h_val_result;
        end  
end
  dof_val_result = size(lr.cp, 2)
end