function plot_vector_field(mesh_V, mesh_T, warped_V, fig_title) 
    tangent_vec = sphere_log_map(warped_V,mesh_V);                      
    
    quiver3(mesh_V(:,1),mesh_V(:,2),mesh_V(:,3), ...
        tangent_vec(:,1), tangent_vec(:,2), tangent_vec(:,3),2, ...
        'LineWidth',1);  % 2 here is the scale parameter, change it to 0 to disable auto-scaling
    hold on
    trisurf(mesh_T, warped_V(:,1), warped_V(:,2), warped_V(:,3), 'FaceColor',[0.8,0.8,0.8], ...
        'EdgeAlpha', 0)
    hold off
    title(fig_title)
    axis off
    axis equal
end