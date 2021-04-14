function [output_structure] = GRAPE_v1

    % GRAPE (Gravity iRregular shAPEs) is a Matlab software to compute the
    % gravitational field of irregular bodies using STL files to describe
    % their shape.
    
   
    % stl_file_name: file name of the stl file (with .stl extension)
    % N_spheres_in: number of spheres to be placed inside the STL file
    % N_iter_single_sphere: number of iterations requested before placing
    % one single sphere
    % m_tot: total mass of the object [kg]
    % rho: density of the spheres [kg/m^3]
    % max_size_m: maximum size of the object [meters]
    
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: importing the STL file for reading the external topography of
    % the object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    
    % Read the input file
    [inp_var] = readInputFile;
    
    % Volume definitions for probes
    data.x_max = inp_var.probes_x_max;
    data.x_min = inp_var.probes_x_min;
    data.y_max = inp_var.probes_y_max;
    data.y_min = inp_var.probes_y_min;
    data.z_max = inp_var.probes_z_max;
    data.z_min = inp_var.probes_z_min;
    
    % Nuber of total spheres
    data.N_spheres_in = inp_var.N_tot_spheres;
    
    % Number of iterations to find the i-esime sphere
    data.N_check_2_max = inp_var.N_iter;
    
    % Particle density
    data.m_tot = inp_var.tot_mass;                                                   
    data.rho_fake = 0;
    
    fv_not_scaled = stlread(inp_var.stl_filename);                                     
   
    dim = size(fv_not_scaled.faces);
    
    figure (1)
    renderSTL(fv_not_scaled);
    title('Original STL file (not scaled)');
    axis equal
    set(gca,'FontSize',24)
    
    % Step 2: creation of the data structure. 
    data.N_points_probe = inp_var.N_probes;
    data.G = 6.674e-11;
    data.N_points_net = dim(1);
    
    
    % Step3: SCALING
    [data, fv] = STLPoints(data, fv_not_scaled, inp_var.max_size);
    data.fv = fv;
    
    figure (2)
    renderSTL(fv);
    title('scaled STL file');
    axis equal
    set(gca,'FontSize',24)
    
    % Step4: find the barycenter
    [data.x_bary, data.y_bary, data.z_bary] = findCentroid(data.x_net, data.y_net, data.z_net);
      
   
    % Step5: loop over N (net) random points within the object
    
    N_check_2 = 1;
    condition = 1;
    condition_2 = 1;
    vol_spheres_tot = 0;
    const_1 = 4/3*pi; 
    alpha_transparency = 0.2;
    
    figure(3)
       
    while condition <= data.N_spheres_in
        
        x_rand = (data.max_x-data.min_x).*rand(1,1) + data.min_x;
        y_rand = (data.max_y-data.min_y).*rand(1,1) + data.min_y;
        z_rand = (data.max_z-data.min_z).*rand(1,1) + data.min_z;
        
        IN = inpolyhedron(data.fv, [x_rand, y_rand, z_rand]);
        
        % If here the random point is inside the polygon
        if (IN == 1)
            
            if condition <= 1
               
                % Creation of the first sphere
                
                [radius, pos_min, check] = findClosestPoint(data.x_net, data.y_net, data.z_net, x_rand, y_rand, z_rand, 0);
                spheres_struct.r(condition) = radius; 
                spheres_struct.vol(condition) = 4/3*pi*spheres_struct.r(condition)^3;
                spheres_struct.x_c(condition) = x_rand; 
                spheres_struct.y_c(condition) = y_rand; 
                spheres_struct.z_c(condition) = z_rand; 
                spheres_struct.id(condition) = condition; 
                
           
                
                % Total volume of the spheres
                vol_spheres_tot = const_1 * spheres_struct.r(condition)^3 + vol_spheres_tot;
                
                [points_x_c, points_y_c, points_z_c] = createSingleSphere(spheres_struct.r(condition), ...
                                                                          spheres_struct.x_c(condition), ...
                                                                          spheres_struct.y_c(condition), ...
                                                                          spheres_struct.z_c(condition));
                    
                hl = surf(points_x_c, points_y_c, points_z_c);
                colormap bone
                shading interp
                alpha(hl, alpha_transparency)
                hold on
                fprintf('Now filling the STL file with sphere number %i\n', condition);
                condition = condition + 1;
                
            else
                
                % Centers of all the spheres within the volume (dynamic
                % structure: it increases with time)
                
                x_pos_center_spheres = spheres_struct.x_c;
                y_pos_center_spheres = spheres_struct.y_c;
                z_pos_center_spheres = spheres_struct.z_c;
                r_spheres = spheres_struct.r;
               
                [dists_spheres, pos_min, check] = findClosestPoint(x_pos_center_spheres, y_pos_center_spheres, z_pos_center_spheres, x_rand, y_rand, z_rand, r_spheres);
            
                % At this stage we ignore if the point is internal to an
                % existing sphere or not: check!
                
                if (check==0)
    
                    [dists_wall, pos_min, check] = findClosestPoint(data.x_net, data.y_net, data.z_net, x_rand, y_rand, z_rand, 0);
                    
                    effective_r = min(dists_spheres, dists_wall);
                    
                    temporary.r(N_check_2) = effective_r;
                    temporary.x_c(N_check_2) = x_rand;
                    temporary.y_c(N_check_2) = y_rand;
                    temporary.z_c(N_check_2) = z_rand;
                    
                    
                    if (N_check_2 > data.N_check_2_max)
                    
                        [dist_net, pos_max] = max(temporary.r);
                        
                        spheres_struct.r(condition) = temporary.r(pos_max); 
                        spheres_struct.vol(condition) = 4/3*pi*spheres_struct.r(condition)^3;
%                         spheres_struct.m(condition) = spheres_struct.vol(condition) * data.rho;
%                         spheres_struct.m(condition) = 1;
                        spheres_struct.x_c(condition) = temporary.x_c(pos_max); 
                        spheres_struct.y_c(condition) = temporary.y_c(pos_max); 
                        spheres_struct.z_c(condition) = temporary.z_c(pos_max); 
                        spheres_struct.id(condition) = condition; 
                    
                        % Total volume of the spheres
                        vol_spheres_tot = const_1 * spheres_struct.r(condition)^3 + vol_spheres_tot;
                    
                        [points_x_c, points_y_c, points_z_c] = ...
                                                           createSingleSphere(spheres_struct.r(condition), ...
                                                                              spheres_struct.x_c(condition), ...
                                                                              spheres_struct.y_c(condition), ...
                                                                              spheres_struct.z_c(condition));
                    
                        hl = surf(points_x_c, points_y_c, points_z_c);
                        colormap bone
                        shading interp
                        alpha(hl, alpha_transparency)
                        hold on
                        title('STL representation with spheres')
                        set(gca,'FontSize',24)
                        
                        fprintf('Now filling the STL file with sphere number %i\n', condition);
                        condition = condition + 1;
                        dists = [];
                        check = [];
                        pos_min = [];
                        temporary = [];
                        N_check_2 = 1;
                        
                    end
                    
                    N_check_2 = N_check_2+1;
                    
                end
                
            end 
            
        end
        
    end % End While
   
    renderSTL(fv);
    axis equal

    data.rho_fake = data.m_tot/vol_spheres_tot;
    data_probes = parallelComputationField(data, spheres_struct);
    
    toc
    
%     quiver3(data_probes.x_rand, data_probes.y_rand, data_probes.z_rand, data_probes.g_x_point, data_probes.g_y_point, data_probes.g_z_point, 1)

    % Streamlines
    N_traj = 20;
    
    % Initial points: left wall
    
    [x_m_1, y_m_1, z_m_1] = meshgrid(data_probes.x_min, linspace(data_probes.y_min, data_probes.y_max, N_traj), linspace(data_probes.z_min, data_probes.z_max, N_traj));
    
    traj_x_1 = x_m_1(:);
    traj_y_1 = y_m_1(:);
    traj_z_1 = z_m_1(:);
    
    
    % Initial points: right wall
  
    [x_m_2, y_m_2, z_m_2] = meshgrid(data_probes.x_max, linspace(data_probes.y_min, data_probes.y_max, N_traj), linspace(data_probes.z_min, data_probes.z_max, N_traj));
    
    traj_x_2 = x_m_2(:);
    traj_y_2 = y_m_2(:);
    traj_z_2 = z_m_2(:);
    
    % Initial points: front wall
    
    [x_m_3, y_m_3, z_m_3] = meshgrid(linspace(data_probes.x_min, data_probes.x_max, N_traj), data_probes.y_min, linspace(data_probes.z_min, data_probes.z_max, N_traj));
    
    traj_x_3 = x_m_3(:);
    traj_y_3 = y_m_3(:);
    traj_z_3 = z_m_3(:);
   
    % Initial points: rear wall
    
    [x_m_4, y_m_4, z_m_4] = meshgrid(linspace(data_probes.x_min, data_probes.x_max, N_traj), data_probes.y_max, linspace(data_probes.z_min, data_probes.z_max, N_traj));
    
    traj_x_4 = x_m_4(:);
    traj_y_4 = y_m_4(:);
    traj_z_4 = z_m_4(:);
    
     % Initial points: bottom wall
 
    [x_m_5, y_m_5, z_m_5] = meshgrid(linspace(data_probes.x_min, data_probes.x_max, N_traj), linspace(data_probes.y_min, data_probes.y_max, N_traj), data_probes.z_min);
    
    traj_x_5 = x_m_5(:);
    traj_y_5 = y_m_5(:);
    traj_z_5 = z_m_5(:);
     
    % Initial points: up wall
   
    [x_m_6, y_m_6, z_m_6] = meshgrid(linspace(data_probes.x_min, data_probes.x_max, N_traj), linspace(data_probes.y_min, data_probes.y_max, N_traj), data_probes.z_max);
    
    traj_x_6 = x_m_6(:);
    traj_y_6 = y_m_6(:);
    traj_z_6 = z_m_6(:);
    
   
    sx = [traj_x_1 traj_x_2 traj_x_3 traj_x_4 traj_x_5 traj_x_6];
    sy = [traj_y_1 traj_y_2 traj_y_3 traj_y_4 traj_y_5 traj_y_6];
    sz = [traj_z_1 traj_z_2 traj_z_3 traj_z_4 traj_z_5 traj_z_6];
    
    x_stream = reshape(data_probes.x_rand, [data.N_points_probe, data.N_points_probe, data.N_points_probe]);
    y_stream = reshape(data_probes.y_rand, [data.N_points_probe, data.N_points_probe, data.N_points_probe]);
    z_stream = reshape(data_probes.z_rand, [data.N_points_probe, data.N_points_probe, data.N_points_probe]);
    
    g_x_stream = reshape(data_probes.g_x_point, [data.N_points_probe, data.N_points_probe, data.N_points_probe]);
    g_y_stream = reshape(data_probes.g_y_point, [data.N_points_probe, data.N_points_probe, data.N_points_probe]);
    g_z_stream = reshape(data_probes.g_z_point, [data.N_points_probe, data.N_points_probe, data.N_points_probe]);
 
    g_tot = sqrt(g_x_stream.^2+g_y_stream.^2+g_z_stream.^2);
    
    figure (4)
    hhh = contourslice(x_stream, y_stream, z_stream, ...
                        g_tot, inp_var.x_slice_fig4, [], []);  
    [vol_true, miao] = renderSTL2(fv);
    axis equal
    title('Gravitational acceleration |g| on the x plane (in m s^{-2})');                      
    set(hhh,'LineWidth', 1.5) 
    hcb4 = colorbar;
    set(gca,'FontSize',24)
    view([90, 0])       
    
    figure (5)
    hhh = contourslice(x_stream, y_stream, z_stream, ...
                        g_tot, [], inp_var.y_slice_fig5, []);  
    [vol_true, miao] = renderSTL2(fv);
    axis equal
    title('Gravitational acceleration |g| on the y plane (in m s^{-2})')                      
    set(hhh,'LineWidth', 1.5) 
    colorbar
    set(gca,'FontSize',24)
    view([0, 0])                  
    
    figure (6)
    hhh = contourslice(x_stream, y_stream, z_stream, ...
                        g_tot, [], [], inp_var.z_slice_fig6);  
    [vol_true, miao] = renderSTL2(fv);
    axis equal
    title('Gravitational acceleration |g| on the z plane (in m s^{-2})')                      
    set(hhh,'LineWidth', 1.5) 
    colorbar
    set(gca,'FontSize',24)
    view([0, 90])  
    
    
    figure (7)                                
    axis equal
    hhh = streamslice(x_stream, y_stream, z_stream, ...
                     g_x_stream, g_y_stream, g_z_stream, ...
                     inp_var.x_slice_fig7, [], [] , 0.3);
    set(hhh, 'Color', [1 0.6 0] ) 
    set(hhh,'LineWidth', 1.5)
    title('Gravitational acceleration |g| - field lines on the x plane')  
    view([90, 0])                        
    [vol_true, miao] = renderSTL2(fv);
    set(gca,'FontSize',24)
    
    
    figure (8)                                
    axis equal
    hhh = streamslice(x_stream, y_stream, z_stream, ...
                     g_x_stream, g_y_stream, g_z_stream, ...
                     [], inp_var.y_slice_fig8, [] , 0.3);
    set(hhh, 'Color', [1 0.6 0] ) 
    set(hhh,'LineWidth', 1.5)            
    title('Gravitational acceleration |g| - field lines on the y plane')  
    [vol_true, miao] = renderSTL2(fv);
    view([0, 0])   
    set(gca,'FontSize',24)
    
    
    figure (9)                                
    axis equal
    hhh = streamslice(x_stream, y_stream, z_stream, ...
                     g_x_stream, g_y_stream, g_z_stream, ...
                     [], [], inp_var.z_slice_fig9 , 0.3);
    set(hhh, 'Color', [1 0.6 0] ) 
    set(hhh,'LineWidth', 1.5)
    title('Gravitational acceleration |g| - field lines on the z plane')  
    view([0, 90])                          
    [vol_true, miao] = renderSTL2(fv);
    set(gca,'FontSize',24)
    
    ratio = vol_spheres_tot/vol_true;
    
    output_structure.data_probes = data_probes;
    
    if (inp_var.boolean_orbit>0)
        data.N_points_probe = 10;
        data_probes = parallelComputationField(data, spheres_struct);
        [output_structure.T_orbit, output_structure.Y_orbit] = getOrbit(data_probes, data, fv, inp_var);
    end
    

end

function [data, fv_scaled] = STLPoints(data, fv_not_scaled, real_d_x_m)

    % Not scaled object
    tab = fv_not_scaled.vertices;
    
    x_all = tab(:,1);                       
    y_all = tab(:,2);
    z_all = tab(:,3);
    
    min_x = min(x_all);
    min_y = min(y_all);
    min_z = min(z_all);
    max_x = max(x_all);
    max_y = max(y_all);
    max_z = max(z_all);
    
    % Scaled image to real dimensions
    diff_pxl = abs(max_x-min_x);
    c_factor = real_d_x_m/diff_pxl;
%     c_factor = 1;
    
    fv_scaled.faces = fv_not_scaled.faces;
    fv_scaled.vertices = c_factor * fv_not_scaled.vertices;
    
    % Working with scaled object
    tab = [];
    x_all = [];
    y_all = [];
    z_all = [];
    min_x = [];
    max_x = [];
    min_y = [];
    max_y = [];
    min_z = [];
    max_z = [];
    
    tab = fv_scaled.vertices;
    
    % Total number of points in the STL file
    x_all = tab(:,1);                       
    y_all = tab(:,2);
    z_all = tab(:,3);
    
    % Number of points in the STL file
    N_x_tot = length(x_all);
    N_y_tot = N_x_tot;
    N_z_tot = N_x_tot;
    data.N_points_tot = N_x_tot;
    
    % Creation of a random permutation of the indexes of all the points
    N_rand_idxs_tot = randperm(N_x_tot);

    % Take the first N elements in the random vector
    N_rand_idxs = N_rand_idxs_tot(1:data.N_points_net);
    
    x_net = x_all(N_rand_idxs);
    y_net = y_all(N_rand_idxs);
    z_net = z_all(N_rand_idxs);
    
    % Min points
    min_x = min(x_net);
    min_y = min(y_net);
    min_z = min(z_net);
    
    % Max points
    max_x = max(x_net);
    max_y = max(y_net);
    max_z = max(z_net);
    
    data.x_net = x_net;
    data.y_net = y_net;
    data.z_net = z_net;
    
    data.min_x = min_x;
    data.min_y = min_y;
    data.min_z = min_z;
    data.max_x = max_x;
    data.max_y = max_y;
    data.max_z = max_z;
    
   
        
end

function [x_c, y_c, z_c] = findCentroid(x, y, z)

    x_c = mean(x);
    y_c = mean(y);
    z_c = mean(z);

end

function [dist_net, pos_min, check, dists] = findClosestPoint(x, y, z, x_c, y_c, z_c, r_c)

    dists = sqrt( (x-x_c).^2 + (y-y_c).^2 + (z-z_c).^2);
    diff = dists-r_c;
    check = any((diff)<0);
    [dist_net, pos_min] = min(abs(diff));

end

function [points_x_c, points_y_c, points_z_c] = createSingleSphere(r, x_c, y_c, z_c)

    % General sphere generator for coating
    [x, y, z] = sphere(20);

    points_x_c = x * r + x_c;
    points_y_c = y * r + y_c;
    points_z_c = z * r + z_c;
    
end

function angle_deg = rad2deg(angle_rad)

    angle_deg = 180 .* angle_rad /pi;

end

function [totalVolume, totalArea] = stlVolume(p,t)

    % Given a surface triangulation, compute the volume enclosed using
    % divergence theorem.
    % Assumption:Triangle nodes are ordered correctly, i.e.,computed normal is outwards
    % Input: p: (3xnPoints), t: (3xnTriangles)
    % Output: total volume enclosed, and total area of surface  
    % Author: K. Suresh; suresh@engr.wisc.edu

    % Compute the vectors d13 and d12
    d13 = [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
    d12 = [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
    cr = cross(d13,d12,1);                                                 % Cross-product (vectorized)
    area = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);                     % Area of each triangle
    totalArea = sum(area);
    crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
    zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
    nz = -cr(3,:)./crNorm;                                                 % z component of normal for each triangle
    volume = area.*zMean.*nz;                                              % Contribution of each triangle
    totalVolume = sum(volume);                                             % Divergence theorem
    
end

function data_probes = parallelComputationField(data, spheres_struct)

    % T-Rex
    data_probes.x_min = data.x_min;
    data_probes.x_max = data.x_max;
    data_probes.y_min = data.y_min;
    data_probes.y_max = data.y_max;
    data_probes.z_min = data.z_min;
    data_probes.z_max = data.z_max;
  

   [x, y, z] = meshgrid(linspace(data_probes.x_min, data_probes.x_max, data.N_points_probe), ...
                        linspace(data_probes.y_min, data_probes.y_max, data.N_points_probe), ...
                        linspace(data_probes.z_min, data_probes.z_max, data.N_points_probe));                 
   
   x_rand = x(:);
   y_rand = y(:);
   z_rand = z(:);
   
   x_rand = x_rand';
   y_rand = y_rand';
   z_rand = z_rand';
    
   
%     k = 1;
%     max_x = data.x_bary+k;
%     min_x = data.x_bary-k;
%     max_y = data.y_bary+k;
%     min_y = data.y_bary-k;
%     max_z = data.z_bary+k;
%     min_z = data.z_bary-k;
%     
%     x_rand = (max_x-min_x).*rand(data.N_points_probe,1) + min_x;
%     y_rand = (max_y-min_y).*rand(data.N_points_probe,1) + min_y;
%     z_rand = (max_z-min_z).*rand(data.N_points_probe,1) + min_z;
%     
%     x_rand = x_rand';
%     y_rand = y_rand';
%     z_rand = z_rand';



%     x_rand_0 = linspace(min_x, max_x, data.N_points_probe);
%     y_rand_0 = linspace(min_y, max_y, data.N_points_probe);
%     z_rand_0 = linspace(min_z, max_z, data.N_points_probe); 
%     
%     [x_rand, y_rand, z_rand] = meshgrid(x_rand_0, y_rand_0, z_rand_0);     

 
    tic
    
    g_x_point = zeros(1,data.N_points_probe);
    g_y_point = zeros(1,data.N_points_probe);
    g_z_point = zeros(1,data.N_points_probe);

    parfor i=1:length(x_rand)
%     parfor i=1:data.N_points_probe
               
        % We use this function just to check if particles are inside
        % existing spheres.
        
        [dists_spheres, pos_min, check] = findClosestPoint(spheres_struct.x_c, ...
                                                           spheres_struct.y_c, ...
                                                           spheres_struct.z_c, ... 
                                                           x_rand(i), y_rand(i), z_rand(i), ...
                                                           spheres_struct.r);
        
        %  ++++++++++++ Random probe is inside an existing sphere +++++++++                                               
                                                       
        if (check==1)
            
            % Contribute of all the external points without the
            % internal one
            
            if (data.N_spheres_in > 1)
           
                % Along x axis
                r_x = x_rand(i) - spheres_struct.x_c(1:end ~= pos_min);
                
                % Along y axis
                r_y = y_rand(i) - spheres_struct.y_c(1:end ~= pos_min);
            
                % Along z axis
                r_z = z_rand(i) - spheres_struct.z_c(1:end ~= pos_min);
            
                r_tot = sqrt(r_x.^2+r_y.^2+r_z.^2);
            
%                 g_x = -data.G*spheres_struct.m(1:end ~= pos_min)./abs(r_tot.^3).*r_x;
                g_x = -data.G*spheres_struct.vol(1:end ~= pos_min)./abs(r_tot.^3).*r_x;
                g_x_point(i) = data.rho_fake*sum(g_x);
            
                % Along y axis
%                 g_y = -data.G*spheres_struct.m(1:end ~= pos_min)./abs(r_tot.^3).*r_y;
                g_y = -data.G*spheres_struct.vol(1:end ~= pos_min)./abs(r_tot.^3).*r_y;
                g_y_point(i) = data.rho_fake*sum(g_y);
            
                % Along z axis
%                 g_z = -data.G*spheres_struct.m(1:end ~= pos_min)./abs(r_tot.^3).*r_z;
                g_z = -data.G*spheres_struct.vol(1:end ~= pos_min)./abs(r_tot.^3).*r_z;
                g_z_point(i) = data.rho_fake*sum(g_z);
            
            else
                
                % If you have just one sphere, the contribute from external
                % points is zero (remember you are inside the sphere)
                g_x_point(i) = 0;
                g_y_point(i) = 0;
                g_z_point(i) = 0;
                
            end
            
            %%% +++++++++++ Here the additional point +++++++++++++++++ %%%
            % This contribute is given by the sphere that contains the
            % point
            
            % Along x axis
            r_x_in = x_rand(i) - spheres_struct.x_c(pos_min);
            
            % Along y axis
            r_y_in = y_rand(i) - spheres_struct.y_c(pos_min);
            
            % Along z axis
            r_z_in = z_rand(i) - spheres_struct.z_c(pos_min);
            
            r_tot_in = sqrt(r_x_in^2+r_y_in^2+r_z_in^2);
            
            m_internal = data.rho_fake * 4/3*pi*r_tot_in^3;
            
            g_x_in = -data.G*m_internal/abs(r_tot_in^3)*r_x_in;
            g_x_point(i) = g_x_point(i) + g_x_in;
            
            g_y_in = -data.G*m_internal/abs(r_tot_in^3)*r_y_in;
            g_y_point(i) = g_y_point(i) + g_y_in;
            
            g_z_in = -data.G*m_internal/abs(r_tot_in^3)*r_z_in;
            g_z_point(i) = g_z_point(i) + g_z_in;
            
        
        else
            
        %  ++++++++ Random probe is NOT inside an existing sphere +++++++++    
            
            % The random point does not belong to any sphere
           
            % Along x axis
            r_x = x_rand(i) - spheres_struct.x_c;
            
            % Along y axis
            r_y = y_rand(i) - spheres_struct.y_c;
            
            % Along z axis
            r_z = z_rand(i) - spheres_struct.z_c;
            
            r_tot = sqrt(r_x.^2+r_y.^2+r_z.^2);
            
%             g_x = -data.G*spheres_struct.m./abs(r_tot.^3).*r_x;
            g_x = -data.G*spheres_struct.vol./abs(r_tot.^3).*r_x;
            g_x_point(i) = data.rho_fake*sum(g_x);
            
            % Along y axis
%             g_y = -data.G*spheres_struct.m./abs(r_tot.^3).*r_y;
            g_y = -data.G*spheres_struct.vol./abs(r_tot.^3).*r_y;
            g_y_point(i) = data.rho_fake*sum(g_y);
            
            % Along z axis
%             g_z = -data.G*spheres_struct.m./abs(r_tot.^3).*r_z;
            g_z = -data.G*spheres_struct.vol./abs(r_tot.^3).*r_z;
            g_z_point(i) = data.rho_fake*sum(g_z);
            
        end
        
    end
    
    toc
    
    data_probes.x_rand = x_rand;
    data_probes.y_rand = y_rand;
    data_probes.z_rand = z_rand;
    data_probes.g_x_point = g_x_point;
    data_probes.g_y_point = g_y_point;
    data_probes.g_z_point = g_z_point;
    
 
%     figure(3)
%     quiver3(x_rand, y_rand, z_rand, g_x_in_point, g_y_in_point, g_z_in_point, 10)
%     hold on
%     plot3(spheres_struct.x_c, spheres_struct.y_c ,spheres_struct.z_c, 'or')
    
end

function [inp_var] = readInputFile

    fid = fopen('input.inp');
   
    stl_filename = textscan(fid,'%s\n','CommentStyle','%');                % Name of the folder for the single run
    stl_filename = num2str(cell2mat(stl_filename{:}));
    inp_var.stl_filename = stl_filename;
    
    N_tot_spheres = textscan(fid,'%s\n','CommentStyle','%');               % Number of spheres to use to fill the STL file [#]
    N_tot_spheres = str2double(N_tot_spheres{:});
    inp_var.N_tot_spheres = N_tot_spheres;
    
    N_iter = textscan(fid,'%s\n','CommentStyle','%');                      % Number of iterations before placing a single sphere [#]
    N_iter = str2double(N_iter{:});
    inp_var.N_iter = N_iter;
    
    tot_mass = textscan(fid,'%s\n','CommentStyle','%');                    % Total mass [kg]
    tot_mass = str2double(tot_mass{:});
    inp_var.tot_mass = tot_mass;
    
    max_size = textscan(fid,'%s\n','CommentStyle','%');                    % Maximum size of the object [m]
    max_size = str2double(max_size{:});
    inp_var.max_size = max_size;
    
    N_probes = textscan(fid,'%s\n','CommentStyle','%');                    % N_probes [#]
    N_probes = str2double(N_probes{:});
    inp_var.N_probes = N_probes;
    
    probes_x_max = textscan(fid,'%s\n','CommentStyle','%');                % Max x coordinate for the probe
    probes_x_max = str2double(probes_x_max{:});
    inp_var.probes_x_max = probes_x_max;
    
    probes_x_min = textscan(fid,'%s\n','CommentStyle','%');                % Min x coordinate for the probe
    probes_x_min = str2double(probes_x_min{:});
    inp_var.probes_x_min = probes_x_min;
    
    probes_y_max = textscan(fid,'%s\n','CommentStyle','%');                % Max y coordinate for the probe
    probes_y_max = str2double(probes_y_max{:});
    inp_var.probes_y_max = probes_y_max;
    
    probes_y_min = textscan(fid,'%s\n','CommentStyle','%');                % Min y coordinate for the probe
    probes_y_min = str2double(probes_y_min{:});
    inp_var.probes_y_min = probes_y_min;
    
    probes_z_max = textscan(fid,'%s\n','CommentStyle','%');                % Max z coordinate for the probe
    probes_z_max = str2double(probes_z_max{:});
    inp_var.probes_z_max = probes_z_max;
    
    probes_z_min = textscan(fid,'%s\n','CommentStyle','%');                % Min z coordinate for the probe
    probes_z_min = str2double(probes_z_min{:});
    inp_var.probes_z_min = probes_z_min;
    
    x_slice_fig4 = textscan(fid,'%s\n','CommentStyle','%');                % x slice fig.4
    x_slice_fig4 = str2double(x_slice_fig4{:});
    inp_var.x_slice_fig4 = x_slice_fig4;
    
    y_slice_fig5 = textscan(fid,'%s\n','CommentStyle','%');                % y slice fig.5
    y_slice_fig5 = str2double(y_slice_fig5{:});
    inp_var.y_slice_fig5 = y_slice_fig5;
    
    z_slice_fig6 = textscan(fid,'%s\n','CommentStyle','%');                % z slice fig.6
    z_slice_fig6 = str2double(z_slice_fig6{:});
    inp_var.z_slice_fig6 = z_slice_fig6;
    
    x_slice_fig7 = textscan(fid,'%s\n','CommentStyle','%');                % x slice fig.7
    x_slice_fig7 = str2double(x_slice_fig7{:});
    inp_var.x_slice_fig7 = x_slice_fig7;
    
    y_slice_fig8 = textscan(fid,'%s\n','CommentStyle','%');                % y slice fig.8
    y_slice_fig8 = str2double(y_slice_fig8{:});
    inp_var.y_slice_fig8 = y_slice_fig8;
    
    z_slice_fig9 = textscan(fid,'%s\n','CommentStyle','%');                % z slice fig.9
    z_slice_fig9 = str2double(z_slice_fig9{:});
    inp_var.z_slice_fig9 = z_slice_fig9;
    
    % Orbital determination
    
    boolean_orbit = textscan(fid,'%s\n','CommentStyle','%');               % Boolean variable to calculate the orbit (0 = no; 1 = yes)
    boolean_orbit = str2double(boolean_orbit{:});
    inp_var.boolean_orbit = boolean_orbit;
    
    x_0 = textscan(fid,'%s\n','CommentStyle','%');                         % Initial condition on the position - x component (x0) [m]
    x_0 = str2double(x_0{:});
    inp_var.x_0 = x_0;
    
    y_0 = textscan(fid,'%s\n','CommentStyle','%');                         % Initial condition on the position - y component (y0) [m]
    y_0 = str2double(y_0{:});
    inp_var.y_0 = y_0;
    
    z_0 = textscan(fid,'%s\n','CommentStyle','%');                         % Initial condition on the position - z component (z0) [m]
    z_0 = str2double(z_0{:});
    inp_var.z_0 = z_0;
    
    vx_0 = textscan(fid,'%s\n','CommentStyle','%');                        % Initial condition on the position - velocity component (vx0) [m/s]
    vx_0 = str2double(vx_0{:});
    inp_var.vx_0 = vx_0;
    
    vy_0 = textscan(fid,'%s\n','CommentStyle','%');                        % Initial condition on the position - velocity component (vy0) [m/s]
    vy_0 = str2double(vy_0{:});
    inp_var.vy_0 = vy_0;
    
    vz_0 = textscan(fid,'%s\n','CommentStyle','%');                        % Initial condition on the position - velocity component (vz0) [m/s]
    vz_0 = str2double(vz_0{:});
    inp_var.vz_0 = vz_0;
    
    t_end = textscan(fid,'%s\n','CommentStyle','%');                       % Final time for the integration [s]
    t_end = str2double(t_end{:});
    inp_var.t_end = t_end;
    
    ode_steps = textscan(fid,'%s\n','CommentStyle','%');                   % Number of steps for the ODE
    ode_steps = str2double(ode_steps{:});
    inp_var.ode_steps = ode_steps;
    
    fclose(fid);
    
end
