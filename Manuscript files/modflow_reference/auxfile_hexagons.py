

def hexagon_noise(noise_radius, hexagon_radius, cg_xm, cg_ym, cg_ID, anisotropy =1, visualization = False):
    
    """
    Create a noise field over the hexagonal grid, normalized between 0 and 1
    
    @params:
        noise_radius        - Required  : correlation radius for noise filter, defines on which scale noise correlates; in [lenuni]
        hexagon_radius      - Required  : radius of hexagons used for tessellation; in [lenuni]
        cg_xm               - Required  : correspondence grid in matrix format, x coordinates
        cg_ym               - Required  : correspondence grid in matrix format, y coordinates
        cg_ID               - Required  : IDs of correspondence grid relating to hexagon cells
        anisotropy          - Required  : ratio of x correlation vs y correlation, default is 1
        visualization       - Optional  : boolean that defines whether the process is printed
    
    Returns
        noise_field         : vector corresponding to hexagon cells, noise values scaled between 0 and 1
    """
    import numpy as np
    from scipy import signal
    
    if isinstance(noise_radius, int) == False:
        # noise_radius is NOT a scalar: select a random value in between
        noise_radius = np.random.uniform(noise_radius[0],noise_radius[1],1)
    
    # Extract grid dimensionsa
    xrange = len(cg_xm[0,:])
    yrange = len(cg_xm[:,0])
    
    # Create random values of excess dimensions, convolution will crop it
    noise = np.random.rand(xrange*10,yrange*10)
    
    # Determine filter size
    size_x = noise_radius/hexagon_radius*np.sqrt(3)/anisotropy
    size_y = noise_radius/hexagon_radius*anisotropy
    
    # Create meshgrid, create filter
    x, y = np.mgrid[-size_x:size_x+1, -size_y:size_y+1]
    g = np.exp(-0.333*(x**2/float(size_x)+y**2/float(size_y)))
    
    # Normalize filter
    filter = g/g.sum()
     
    # Convolve noise and filter
    convolved_noise = signal.convolve(noise,filter,mode='valid')[:yrange,:xrange]
    
    # Normalize range
    convolved_noise = (convolved_noise - convolved_noise.min())/(convolved_noise.max() - convolved_noise.min())
    
    # Reshape to vector
    convolved_noise_vec = np.reshape(convolved_noise, xrange*yrange)
    
    if visualization == True:
        
        import matplotlib.pyplot as plt
        
        # Plot Noise
        plt.figure()
        plt.imshow(noise)
        plt.title('raw noise')
        plt.show()
        
        # Plot Filter
        plt.figure()
        plt.imshow(filter)
        plt.show()
        plt.title('filter')
        
        # Plot Noise over grid
        plt.figure()
        plt.scatter(np.reshape(cg_xm,(xrange*yrange)),
                    np.reshape(cg_ym,(xrange*yrange)),
                    c = 'xkcd:light grey')
        plt.scatter(np.reshape(cg_xm,(xrange*yrange))[cg_ID],
                    np.reshape(cg_ym,(xrange*yrange))[cg_ID],
                    c = convolved_noise_vec[cg_ID])
        plt.show()
        plt.title('noise over correspondence grid')
        
    # Reduce noise over correspondence grid to cell vector
#    global noise_field
    noise_field = convolved_noise_vec[cg_ID]
    
    noise_field -= np.min(noise_field)
    noise_field /= np.max(noise_field)
    
    # Return result
    return noise_field