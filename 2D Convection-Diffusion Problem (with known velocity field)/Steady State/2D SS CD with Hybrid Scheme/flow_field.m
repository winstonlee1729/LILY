% 1. Load the data 
phi = load('flow_field.txt');

% 2. Generate the plot
figure;
contourf(phi, 100, 'LineColor', 'none');

% 3. Apply standard CFD styling
colormap('jet');         
c = colorbar;            
ylabel(c, 'Scalar \phi'); 

title('2D SS Convection-Diffusion without Source with Hybrid Scheme', 'FontSize', 14);
xlabel('Grid Node X ');
ylabel('Grid Node Y ');

% 4. Correct the physical orientation
% In your Fortran code, Row 1 is the North boundary. 
% We reverse the Y-axis so Row 1 appears at the top of the screen.
set(gca, 'YDir', 'reverse');

% Lock the aspect ratio so the 1x1 domain looks like a perfect square
axis equal tight;