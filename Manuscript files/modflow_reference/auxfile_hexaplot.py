"""
This library contains several functions designed to help with the illustration of hexagonal grids

Functions:
    plot_hexagaons              : plots a specified data vector over a 2-D hexagon grid.
    create_alpha_mask           : creates an alpha shape (a concave hull), which is required for plotting contours; without it, the contour function extrapolates outside of the model area.
    plot_scattered_contour      : plots contour lines over an irregular grid, such as a hexagonal one.
    plot_hexagons_3d            : plots a 2-dimensional hexagon grid with specified z-dimensions

"""



def plot_hexagons (data, hexagon_grid_cores, hexagon_radius, hexagon_orientation = 0, colormap = 'steel', color = None, vmin = None, vmax = None, vincr = None, xlabel = None, ylabel = None, clabel = None, hide_colorbar = False,  **kwargs):
    
    """
    Call to plot a specified vector (positions relative to node IDs) in a hexagonal grid
    
    @params:
        data                - Required  : vector of values for hexagonal plot, positions corresponding to cell IDs (counting from zero)
        hexagon_grid_cores  - Required  : tessellated polygons over area of interest
        hexagon_radius      - Required  : radius of hexagons used for tessellation
        hexagon_orientation - Optional  : orientation of hexagon in clock-wise degrees [0 = flat top]
        colormap            - Optional  : specify a colormap as string
        vmin                - Optional  : externally specified min value for colorbar
        vmax                - Optional  : externally specified max value for colorbar
        vincr               - Optional  : specified value increment for colorbar
        xlabel              - Optional  : string for xlabel
        ylabel              - Optional  : string for ylabel
        clabel              - Optional  : string for colorbar label
        **kwargs            - Optional  : keyword arguments for matplotlib.patches.RegularPolygon
    """  
    import matplotlib
    import numpy as np
    import math
    
    #--------------------------------------------------------------------------
    # Prepare data for plotting
    #--------------------------------------------------------------------------
    
    # If not specified, define range of values
    if vmin == None or vmax == None:
        vmin = np.min(data)
        vmax = np.max(data)
    vrange = vmax-vmin
    if vincr == None:
        vincr = vrange/100
    
    # Snap value range to integers
    vmin = int(vmin/vincr)*vincr            # minimum value for colorbar
    vmax = (int(vmax/vincr)+1)*vincr        # maximum value for colorbar
    
    if color is None:
    
        # Retrieve colormap
        if colormap == 'steel':
            # Create colormap 'steel'
            from matplotlib.colors import LinearSegmentedColormap
            cmap_steel = [(0.007843137,0.305882353,0.443137255), (0.301960784,0.592156863,0.784313725),(0.623529412,0.776470588,0.882352941)]
            cm = LinearSegmentedColormap.from_list('steel', cmap_steel, N=100)
            cmaps = cm
        else:
            cmaps = colormap

    # Correct orientation
    orientation = math.radians(-hexagon_orientation+30)
    
    # Hexagon radius only goes to normal of sides
    edgepoint_distance = hexagon_radius/np.cos(np.deg2rad(30))
    
    # Retrieve colormap information
    if color is None:
        cmap = matplotlib.cm.get_cmap(cmaps)
            
    #--------------------------------------------------------------------------
    # Start plotting
    #--------------------------------------------------------------------------   
    
    # Create empty figure
    ax1 = matplotlib.pyplot.gca()
    
    # Plot hexagons
    for hex in range(len(hexagon_grid_cores[:,0])):
        
        # Retrieve color value
        if color is None:
            rgba = cmap((data[hex]-vmin)/(vrange))
            rgba = matplotlib.colors.rgb2hex(rgba)
        else:
            rgba = color
        
        # Add the patch
        ax1.add_patch(
        matplotlib.patches.RegularPolygon(
            (hexagon_grid_cores[hex,0], hexagon_grid_cores[hex,1]), # x and y
            6,                                                      # edges
            edgepoint_distance,
            orientation=orientation,
            facecolor = rgba,
            **kwargs)
        )
    
    # Determine meaningful colorbar steps
    if color is None:
        colorbar_increment = vincr
        colorbar_min = int(vmin/colorbar_increment)*colorbar_increment            # minimum value for colorbar
        colorbar_max = (int(vmax/colorbar_increment)+1)*colorbar_increment        # maximum value for colorbar
        colorbar_increment_numbers = int((colorbar_max-colorbar_min)/colorbar_increment+1)
        colorbar_steps = []
        for num in range(colorbar_increment_numbers):
            colorbar_steps = colorbar_steps + [colorbar_min+num*colorbar_increment]
    
    # Recompute the ax.dataLim
    ax1.relim()
    # Update ax.viewLim using the new dataLim
    ax1.autoscale_view()
    
    # Create colorbar
    if hide_colorbar == False and color is None:
        norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
        sm = matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = matplotlib.pyplot.colorbar(sm)
    
    # Label plot
    if xlabel != None:
        matplotlib.pyplot.xlabel(xlabel)
    if ylabel != None:
        matplotlib.pyplot.ylabel(ylabel)
    if clabel != None and not hide_colorbar and color is None:
        cbar.set_label(clabel, rotation=270, labelpad=20)
                
def create_alpha_mask(points, distance_limit, resolution_x = 1000, resolution_y = 1000, visualization = True):
        
        """
        Creates interpolation grid, then masks over the alpha shape spanned up by points and defined by distance_limit.
        
        @params:
            points              - Required  : points spanning up alpha shape
            distance_limit      - Required  : distance threshold for removing Delaunay simplices  
            resolution_x        - Optional  : resolution for grid in x, default is 1000
            resolution_y        - Optional  : resolution for grid in y, default is 1000
            visualization       - Optional  : boolean for visualizing result, default is False
            
        Returns:
            grid_mask           : An array containing 1 for cells inside, and 0 for cells outside
        """  
        
        import numpy as np
        from scipy.spatial import Delaunay
        from matplotlib.collections import LineCollection
        import matplotlib.path as mplPath
        
        #----------------------------------------------------------------------
        # Create Grid
        #----------------------------------------------------------------------        
    
        # Create meshgrid
        xi = np.transpose(np.linspace(min(points[:,0]), max(points[:,0]), resolution_x))
        yi = np.transpose(np.linspace(min(points[:,1]), max(points[:,1]), resolution_y))
        X, Y = np.meshgrid(xi, yi)
        
        # Reshape into vector
        gridpoints_x = np.reshape(X, resolution_x*resolution_y)
        gridpoints_y = np.reshape(Y, resolution_x*resolution_y)
        
        # Combine into gridpoints array
        gridpoints = np.transpose(np.asarray((gridpoints_x, gridpoints_y)))
    
        #----------------------------------------------------------------------
        # Create Alpha Shape
        #----------------------------------------------------------------------
        
        # Start Delaunay triangulation
        tri = Delaunay(points)
        
        # Auxiliary function for plotting, if required
        if visualization == True:
            import matplotlib.pyplot as plt
            edges = set()
            edge_points = []
            def add_edge(i, j):
                """Add a line between the i-th and j-th points, if not in the list already"""
                if (i, j) in edges or (j, i) in edges:
                    # already added
                    return
                edges.add( (i, j) )
                edge_points.append(points[ [i, j] ])
    
        # Remove simplices outside of distance_limit
        simplex_flag = np.zeros(len(tri.simplices[:,0]))    # Flags bad simplices
        counter = 0
        for ia, ib, ic in tri.vertices:
            # ia, ib, ic = indices of corner points of the triangle
            if np.sqrt((points[ia,0]-points[ib,0])**2+(points[ia,1]-points[ib,1])**2) < distance_limit and \
                np.sqrt((points[ia,0]-points[ic,0])**2+(points[ia,1]-points[ic,1])**2) < distance_limit and \
                np.sqrt((points[ib,0]-points[ic,0])**2+(points[ib,1]-points[ic,1])**2) < distance_limit:
                # do nothing
                simplex_flag[counter] = 0
            else:
                # simplex has at least one side larger than threshold, flag it
                simplex_flag[counter] = 1
            counter += 1
        tri.simplices = tri.simplices[simplex_flag == 0,:]  # Remove bad simplices
        tri.vertices = tri.vertices[simplex_flag == 0,:]    # Remove bad simplices
        
        # Visualize, if requested
        if visualization == True:
            # Mark all remaining simplices
            for ia, ib, ic in tri.vertices:
                add_edge(ia, ib)
                add_edge(ib, ic)
                add_edge(ic, ia)
            # Draw them
            lines = LineCollection(edge_points)
            plt.figure()
            plt.gca().add_collection(lines)
            plt.plot(points[:,0], points[:,1], 'o')
            
        #----------------------------------------------------------------------
        # Mask over Alpha Shape
        #----------------------------------------------------------------------
            
        # Prepare point flag
        flag_gridpoints = np.zeros(len(gridpoints[:,0]), dtype = np.int)
        
        # Evaluate gridpoints
        for sim in range(len(tri.simplices[:,0])):
            
            # Print progress bar
            cv = sim
            mv = len(tri.simplices[:,0])-1
            print('\r%s |%s| %s%% %s' % ('Masking: ', '\033[33m'+'█' * int(50 * cv // mv) + '-' * (50 - int(50 * cv // mv))+'\033[0m', ("{0:." + str(1) + "f}").format(100 * (cv / float(mv))), ' Complete'), end = '\r')
            
            # Create simplex path
            bbPath = mplPath.Path(np.array([points[tri.simplices[sim,0],:],
                                        points[tri.simplices[sim,1],:],
                                        points[tri.simplices[sim,2],:],
                                        points[tri.simplices[sim,0],:]]))
            
            
            
            # Flag points that are inside this simplex
            for gridpts in range(len(gridpoints[:,0])):
                if flag_gridpoints[gridpts] == 0: # only process points not already allocated
                    if bbPath.contains_point((gridpoints[gridpts,0],gridpoints[gridpts,1])) == True:
                        flag_gridpoints[gridpts] = 1
        
        # Plot, if required
        if visualization == True:
            plt.scatter(gridpoints[flag_gridpoints == 1,0], gridpoints[flag_gridpoints == 1,1],color = 'g')
            plt.scatter(gridpoints[flag_gridpoints == 0,0], gridpoints[flag_gridpoints == 0,1],color = 'r')
        
        # Reshape flag_gridpoints into a 2D array
        global grid_mask
        grid_mask = np.reshape(flag_gridpoints,(resolution_y,resolution_x))
        
        # Return result
        return grid_mask
    

      
def plot_scattered_contour(x, y, data, resolution_x=1000, resolution_y=1000, 
    grid_mask = None, vmin = None, vmax = None, vincr = None, suppress_clabel = False,
    **kwargs):
    
    """
    Call to plot contour of scattered data
    
    @params:
        x                   - Required  : x-coordinate
        y                   - Required  : y-coordinate
        data                - Required  : data for the contours
        resolution_x        - Optional  : resolution of auxiliary grid in x
        resolution_y        - Optional  : resolution of auxiliary grid in y
        grid_mask           - Optional  : mask array of dimension [resolution_y,resolution_x]
        vmin                - Optional  : min value for contour
        vmax                - Optional  : max value for contour
        vincr               - Optional  : increment for contour
        suppress_clabel     - Optional  : Flag wether contours should be labeld, False by default
        **kwargs            - Optional  : keyword arguments for matplotlib.patches.RegularPolygon
    """  
    
    import numpy as np
    import matplotlib
    import scipy
    
    #--------------------------------------------------------------------------
    # Integrity checks
    #--------------------------------------------------------------------------
    
    # Check if grid_mask matches meshgrid dimensions
    if len(grid_mask) != 1:
        if len(grid_mask[:,0]) != resolution_y or len(grid_mask[0,:]) != resolution_x:
            raise Exception('Grid mask dimensions must match resolution in x and y!')
            
    # Check if one of the cells has dried; this algorithm can't handle that yet
    if vmin < -1000:
        print('\033[31m'+'WARNING:'+'\033[0m'+' Dried cells detected. Contour not printed.')
        return
    
    # Extract vmin and vmax, if not specified
    if vmin == None or vmax == None:
        vmin = np.min(data)
        vmax = np.max(data)
    # Set vincr, if not specified
    if vincr == None:
        vincr = (vmax-vmin)/10
    
    # Snap value range to integers
    vmin = int(vmin/vincr)*vincr            # minimum value for colorbar
    vmax = (int(vmax/vincr)+1)*vincr        # maximum value for colorbar
        
    #--------------------------------------------------------------------------
    # Prepare data for plotting
    #--------------------------------------------------------------------------
    
    # Convert source material into required format
    source = np.transpose(np.asarray([x,y]))
    
    # Create and convert target material
    xi = np.transpose(np.linspace(min(x), max(x), resolution_x))
    yi = np.transpose(np.linspace(min(y), max(y), resolution_y))
    X, Y = np.meshgrid(xi, yi)
    target = np.transpose(np.asarray([X,Y]))
    
    # Interpolate and transpose
    Z = scipy.interpolate.griddata(source, data, target)
    Z = np.transpose(Z)
    
    # Mask values, if grid_mask was specified
    if len(grid_mask) != 1:
        Z[grid_mask == 0] = float('NaN')

    # Define function for masking
    levels = np.arange(vmin,vmax,vincr)
    
    #--------------------------------------------------------------------------
    # Plot that shit
    #--------------------------------------------------------------------------
    
    CS = matplotlib.pyplot.contour(xi,yi,Z,levels=levels,**kwargs)
    if not suppress_clabel:
        matplotlib.pyplot.clabel(CS, inline=1, inline_spacing = 0)
        
    return

def plot_hexagons_3d(grid, zdim, hexagon_radius, hexagon_orientation = 0, xlabel = 'x', ylabel = 'y', zlabel = 'z', clabel = 'depth', depth_colormap = 'steel', alpha = 1, **kwargs):
    
    """
    Call to tessellate a given polygon with hexagons
    
    @params:
        grid                - Required  : x-y-coordinates of center of hexagons, array of form [nx2]
        zdim                - Required  : bottom and top elevation of hexagon cells, array of form [nx2]
        hexagon_radius      - Required  : radius of hexagons used for tessellation
        hexagon_orientation - Required  : orientation of hexagon in clock-wise degrees [0 = flat top]
        xlabel              - Optional  : label for x-axis
        ylabel              - Optional  : label for y-axis
        zlabel              - Optional  : label for z-axis
        clabel              - Optional  : label for colorbar
        depth_colormap      - Optional  : string of colormap, if requested
        alpha               - Optional  : alpha value for transparency of polygons, default is 1
        **kwargs            - Optional  : keyword arguments for Poly3DCollection
    """
    
    # PLOT 3D
    
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    import math
    
    if depth_colormap == 'steel':
        # Create colormap 'steel'
        from matplotlib.colors import LinearSegmentedColormap
        cmap_steel = [(0.007843137,0.305882353,0.443137255), (0.301960784,0.592156863,0.784313725),(0.623529412,0.776470588,0.882352941)]
        cm = LinearSegmentedColormap.from_list('steel', cmap_steel, N=100)
        cmaps = cm
    else:
        cmaps = depth_colormap
    
    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Hexagon radius only goes to normal of sides
    edgepoint_distance = hexagon_radius/np.cos(np.deg2rad(30))
    
    # Determine depth range, if colorbar is requested
    vmin = np.min(zdim[:,1]-zdim[:,0])
    vmax = np.max(zdim[:,1]-zdim[:,0])
    c_range = vmax-vmin
    
    # Plot hexagons
    for hex in range(len(grid[:,0])):
        
        # Reset coordinate variables
        x = []
        y = []
        
        # Read top and bottom elevation
        zbot = zdim[hex,0]
        ztop = zdim[hex,1]
        
        # Pre-allocate memory for coordinate matrix
        Z = np.zeros((12,3))
        
        # Determine cell color, if requested
        if depth_colormap != 'None':
            import matplotlib
            # Retrieve colormap information
            cmap = matplotlib.cm.get_cmap(cmaps)
            rgba = cmap((ztop-zbot-vmin)/c_range) #cmap((zbot-vmin)/(vmax-vmin))
            rgba = list(rgba)
            rgba[3] = alpha
#            rgba = matplotlib.colors.rgb2hex(rgba)
        
        # Plot grid
        counter = 0
        for angle in range(0-hexagon_orientation, 420-hexagon_orientation, 60):
            
            # Coordinates of edge point
            x = np.append(x,grid[hex,0]+math.cos(math.radians(angle)) * edgepoint_distance)
            y = np.append(y,grid[hex,1]+math.sin(math.radians(angle)) * edgepoint_distance)
        
            # Write into coordinate matrix
            if counter < 6:
                Z[counter,0] = grid[hex,0]+math.cos(math.radians(angle)) * edgepoint_distance
                Z[counter,1] = grid[hex,1]+math.sin(math.radians(angle)) * edgepoint_distance
                Z[counter,2] = ztop
                Z[6+counter,0] = grid[hex,0]+math.cos(math.radians(angle)) * edgepoint_distance
                Z[6+counter,1] = grid[hex,1]+math.sin(math.radians(angle)) * edgepoint_distance
                Z[6+counter,2] = zbot
            
            counter += 1
    
        # Vertices of hexagon sides
        verts = [[Z[0],Z[1],Z[7],Z[6]],
                 [Z[1],Z[2],Z[8],Z[7]], 
                 [Z[2],Z[3],Z[9],Z[8]], 
                 [Z[3],Z[4],Z[10],Z[9]], 
                 [Z[4],Z[5],Z[11],Z[10]],
                 [Z[5],Z[0],Z[6],Z[11]]]
        
        
        if depth_colormap != 'None':
            # Plot hexagon side
            face = Poly3DCollection(verts, 
                                     **kwargs)
            face.set_facecolor(rgba)
            ax.add_collection3d(face)
        else:
            face = Poly3DCollection(verts, 
                                     **kwargs)
            face.set_facecolor(rgba)
            ax.add_collection3d(face)
            
        # Vertices of hexagon top
        verts = [[Z[0],Z[1],Z[2],Z[3],Z[4],Z[5]]]
        # Plot hexagon top
        if depth_colormap != 'None':
            # Plot hexagon side
            face = Poly3DCollection(verts, 
                                     **kwargs)
            face.set_facecolor(rgba)
            ax.add_collection3d(face)
        else:
            face = Poly3DCollection(verts, 
                                     **kwargs)
            face.set_facecolor(rgba)
            ax.add_collection3d(face)
        
        # Vertices of hexagon bot
        verts = [[Z[6],Z[7],Z[8],Z[9],Z[10],Z[11]]]
        # Plot hexagon bot
        if depth_colormap != 'None':
            # Plot hexagon side
            face = Poly3DCollection(verts, 
                                     **kwargs)
            face.set_facecolor(rgba)
            ax.add_collection3d(face)
        else:
            face = Poly3DCollection(verts, 
                                     **kwargs)
            face.set_facecolor(rgba)
            ax.add_collection3d(face)
    
    # Determine meaningful colorbar steps, if colorbar was requested
    if depth_colormap != 'None':
        colorbar_increment = 0.1
        colorbar_min = int(vmin/colorbar_increment)*colorbar_increment            # minimum value for colorbar
        colorbar_max = (int(vmax/colorbar_increment)+1)*colorbar_increment        # maximum value for colorbar
        colorbar_increment_numbers = int((colorbar_max-colorbar_min)/colorbar_increment+1)
        colorbar_steps = []
        for num in range(colorbar_increment_numbers):
            colorbar_steps = colorbar_steps + [colorbar_min+num*colorbar_increment]
        # Create colorbar
        norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
        sm = matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = matplotlib.pyplot.colorbar(sm)
        cbar.set_label(clabel, rotation=270, labelpad=20)
    
    # Label axes
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    
    # Equal aspect scaling doesn't work yet, manual workaround
    # Designate array of edges
    xyzlims = np.zeros((3,2))
    xyzlims[0,0] = np.min(grid[:,0])
    xyzlims[0,1] = np.max(grid[:,0])
    xyzlims[1,0] = np.min(grid[:,1])
    xyzlims[1,1] = np.max(grid[:,1])
    xyzlims[2,0] = np.min(zdim)
    xyzlims[2,1] = np.max(zdim)
    
    # Determine maximal range
    maxrange = np.max([xyzlims[0,1]-xyzlims[0,0],xyzlims[1,1]-xyzlims[1,0],xyzlims[2,1]-xyzlims[2,0]])
    
    # Determine difference to maximal range
    xdif = maxrange - (xyzlims[0,1]-xyzlims[0,0])
    ydif = maxrange - (xyzlims[1,1]-xyzlims[1,0])
    zdif = maxrange - (xyzlims[2,1]-xyzlims[2,0])
    
    # Set axis limits -> equal aspect
    ax.set_xlim3d(xyzlims[0,0]-xdif/2,xyzlims[0,1]+xdif/2)
    ax.set_ylim3d(xyzlims[1,0]-ydif/2,xyzlims[1,1]+ydif/2)
    ax.set_zlim3d(xyzlims[2,0]-zdif/2,xyzlims[2,1]+zdif/2)
    
    # Show result
    plt.show()
    
def vulture_plot(incr = 1, elev = 40., fps = 50):
    
    """
    Creates a short animated .gif providing a flight around the 3-D model, requiring an open, compatible 3D figure
    
    @params:
        incr                - Optional  : degree increment for rotation frames; defines temporal resolution of .gif (default = 1)
        elev                - Optional  : elevation angle for camera (default = 40)
        fps                 - Optional  : frames per second for resulting .gif; defines speed of .gif display (default 50)
    """
    
    # Import libraries
    import imageio
    import os
    import matplotlib.pyplot as plt
    
    # Retrieve axis
    ax = plt.gca()
    
    # Rotate, save and compile vulture plot
    images = []
    for cv in range(0,360,incr):
        
        # Rotate image
        ax.view_init(elev=40., azim=cv)
        plt.show()
        
        # Save it as temporary file
        plt.savefig("dummy.png")
        
        # Append it to saved movie
        images.append(imageio.imread("dummy.png"))
        
        # Remove temporary file
        os.remove("dummy.png")

        # Print progress bar
        mv = 359 # max value
        print('\r%s |%s| %s%% %s' % ('Printing: ', '\033[33m'+'█' * int(50 * cv // mv) + '-' * (50 - int(50 * cv // mv))+'\033[0m', ("{0:." + str(1) + "f}").format(100 * (cv / float(mv))), ' Complete'), end = '\r')
    
    # Compile .gif
    imageio.mimsave('output_quick.gif', images,fps=fps)
    
def visualize_genealogy(genealogy,weights = None, rejuvenation = None,colormap = 'jet'):
    
    """
    Creates an inline figure visualizing the particle genealogy over one resampling step.
    
    @params:
        genealogy           - Required  : vector describing genealogy of resampled particles, referring to indices
        weights             - Optional  : weight of particles prior to resampling
        rejuvenation        - Optional  : vector of booleans describing whether particles were rejuvenated
        colormap            - Optional  : colormap string for visualization
    """
    
    import numpy as np
    from IPython import get_ipython
    import matplotlib
    import matplotlib.pyplot as plt
    
    # Determine number of particles
    n_particles = len(genealogy)
    
    # Assign optional variables, if not provided
    if weights is None == True:
        weights = np.ones(n_particles)
#    if rejuvenation is None == True:
#        rejuvenation = np.ones((n_particles),dtype = np.bool)
        
    # Switch to inline printing
    get_ipython().run_line_magic('matplotlib', 'inline')
    
    # Create dummy features for the legend
    full_line = plt.Line2D([], [], color='black',label='inherited')
    dashed_line = plt.Line2D([], [], linestyle = '--', color='black',label='rejuvenated')
    particle = plt.Line2D([], [], linestyle = 'None', marker ='.', color='black',label='particle')
    
    # Plot legend
    plt.legend(handles=[dashed_line,full_line,particle],bbox_to_anchor=(0., -0.05, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)

    # Determine colormap for particles
    cmap = matplotlib.cm.get_cmap(colormap)
    
    # Extract particle colors
    rgba = [None] * n_particles
    for n in range(n_particles): 
        rgba[n] = matplotlib.colors.rgb2hex(cmap(n/(n_particles-1)))
        
    # Create plot
    for n in range(n_particles): 
        plt.plot([genealogy[n],n],[1,2],'--',c=rgba[genealogy[n]])
        # Draw genealogy of current particle
#        if rejuvenation[n] == False:
#            plt.plot([genealogy[n],n],[1,2],c=rgba[genealogy[n]])
#        else:
#            plt.plot([genealogy[n],n],[1,2],c='w')
#            plt.plot([genealogy[n],n],[1,2],'--',c=rgba[genealogy[n]])
            
        # Scatter previous and current particle index
        if weights[n] == 0: # Particle weight is zero - print as greyscale
            plt.scatter(n,1,s = weights[n]/np.max(weights)*55+5,c='xkcd:medium grey') 
        else:
            plt.scatter(n,1,s = weights[n]/np.max(weights)*55+5,c=rgba[n])      
        plt.scatter(n,2,s=20,c=rgba[n])
    
    # Deactivate axes
    plt.axis('off')
    
    # Show, and revert to automatic printing
    plt.show()    
    get_ipython().run_line_magic('matplotlib', 'qt5')