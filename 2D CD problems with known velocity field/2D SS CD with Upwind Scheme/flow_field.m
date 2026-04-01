% 1. Load the data 
phi = load('flow_field.txt');

% 2. Generate the plot
figure;
contourf(phi, 100, 'LineColor', 'none');

% 3. Apply standard CFD styling
colormap('jet');        
c = colorbar;            
ylabel(c, 'Scalar \phi'); 

title('2D SS Convection-Diffusion without Source with Upwind Scheme', 'FontSize', 14);
xlabel('Grid Node X ');
ylabel('Grid Node Y ');

% 4. Correct the physical orientation
set(gca, 'YDir', 'reverse');
axis equal tight;