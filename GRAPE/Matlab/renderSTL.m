function totalVolume = renderSTL(fv)

%     figure(1)
    hold on;
    patch(fv,'FaceColor', [0.8 0.8 1.0], ...
         'EdgeColor', 'none', ...
         'FaceLighting', 'gouraud', ...
         'AmbientStrength', 0.85);

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');
    alpha(0.2)                                                             % Flag: ho aggiunto questo solo per vedere le sfere dentro
    axis('image');
%     view([-135 35]);
    view([180 0]);
    
    [totalVolume, totalArea] = stlVolume(fv.vertices', fv.faces');
    
end