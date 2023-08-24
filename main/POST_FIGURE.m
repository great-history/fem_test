function POST_FIGURE(NE,NODE,U)
    colormap('default');
    patch('Faces',NE,'Vertices',NODE,'FaceVertexCData',U,'FaceColor','interp','EdgeAlpha',1);
    grid off
    colorbar;
end