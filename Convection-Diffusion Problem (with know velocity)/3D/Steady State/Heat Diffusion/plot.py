import numpy as np
import plotly.graph_objects as go

def plot_full_3d_contour(filename='temperature_field.txt'):
    data = np.loadtxt(filename)

    nx, ny, nz = 100, 100, 100
    x_3d = data[:, 0].reshape((nz, ny, nx))
    y_3d = data[:, 1].reshape((nz, ny, nx))
    z_3d = data[:, 2].reshape((nz, ny, nx))
    phi_3d = data[:, 3].reshape((nz, ny, nx))

    step = 2 
    x_down = x_3d[::step, ::step, ::step].flatten()
    y_down = y_3d[::step, ::step, ::step].flatten()
    z_down = z_3d[::step, ::step, ::step].flatten()
    phi_down = phi_3d[::step, ::step, ::step].flatten()
    
    fig = go.Figure(data=go.Volume(
        x=x_down,
        y=y_down,
        z=z_down,
        value=phi_down,
        isomin=np.min(phi_down),
        isomax=np.max(phi_down),
        opacity=0.1,       
        surface_count=100,  
        colorscale='jet',
    ))
    
    fig.update_layout(
        title='3D Convection-Diffusion Temperature Field',
        scene=dict(
            xaxis_title='x (m)',
            yaxis_title='y (m)',
            zaxis_title='z (m)'
        )
    )

    fig.show()

if __name__ == '__main__':
    plot_full_3d_contour()