function [T_orbit, Y_orbit] = getOrbit(data_probes, data, fv, inp_var)

    data.Y0 = [inp_var.x_0 inp_var.y_0 inp_var.z_0 inp_var.vx_0 inp_var.vy_0 inp_var.vz_0];
    data.T_end = inp_var.t_end;
    data.T_steps = inp_var.ode_steps;

    tic
    disp('Orbit evaluation')
    [T_orbit, Y_orbit] = getOrbitInnerFunction(data_probes, data);
    toc
    
    sphere_points = 400;
    Rs = 0.3;
    [X_sphere_S Y_sphere_S Z_sphere_S] = sphere(sphere_points);            % Reference sphere for Sun
    
    
    figure(10)
    hold on
    renderSTL2(fv);

    plot3(Y_orbit(1,1),Y_orbit(1,2),Y_orbit(1,3),'s','MarkerSize', 15)
    plot3(Y_orbit(:,1),Y_orbit(:,2),Y_orbit(:,3),'-','MarkerSize', 15)  

end

function [T, Y] = getOrbitInnerFunction(data_probes, data)

    % Function for the determination of the orbit

    Y0 = data.Y0;
    N = 1;
    T_end = data.T_end;
    
    % Reduced number of points in the data structure
    
    data_probes_reduced = data_probes;
    
    tspan = linspace(0, T_end, data.T_steps);
    
    options = odeset('RelTol',1e-4,'AbsTol',1e-4);
    table(:,1) = data_probes_reduced.x_rand;
    table(:,2) = data_probes_reduced.y_rand;
    table(:,3) = data_probes_reduced.z_rand;
    table(:,4) = data_probes_reduced.g_x_point;
    table(:,5) = data_probes_reduced.g_y_point;
    table(:,6) = data_probes_reduced.g_z_point;
    [T, Y]= ode113(@newtode, tspan, Y0, options, data_probes_reduced, table);

end

function dydt = newtode(t, Y, data_probes_reduced, table)
   
    % equations of motion governing (coupled) dynamics of the satellite
    t
    dydt1 = Y(4);
    dydt2 = Y(5);
    dydt3 = Y(6);
    
    dydt4 = griddata(data_probes_reduced.x_rand, data_probes_reduced.y_rand, data_probes_reduced.z_rand, data_probes_reduced.g_x_point, ...
                     Y(1), Y(2), Y(3));
                 
    dydt5 = griddata(data_probes_reduced.x_rand, data_probes_reduced.y_rand, data_probes_reduced.z_rand, data_probes_reduced.g_y_point, ...
                     Y(1), Y(2), Y(3));             
    
    dydt6 = griddata(data_probes_reduced.x_rand, data_probes_reduced.y_rand, data_probes_reduced.z_rand, data_probes_reduced.g_z_point, ...
                     Y(1), Y(2), Y(3));             
                             

    dydt = [dydt1 dydt2 dydt3 dydt4 dydt5 dydt6] ;
    
    dydt = dydt';
    
end