def rotate(polygon,angle):
    """
    Rotates polygon clock-wise around angle
    
    @params:
        polygon     - Required  : coordinates of closed two-dimensional n-gon over which we want to tessellate, as array of dimensions [n,2]
        angle       - Required  : rotation angle around origin (0,0) in degrees
    """
    
    import math
    import numpy as np
    
    # Prepare variables
    angle = -math.radians(angle)
    polygon_rotated = []
    
    # Rotate each corner point
    for pt in polygon :
        polygon_rotated.append(( pt[0]*math.cos(angle)-pt[1]*math.sin(angle) , pt[0]*math.sin(angle)+pt[1]*math.cos(angle)) )
    
    polygon_rotated = np.asarray(polygon_rotated)
    
    # Return result
    return polygon_rotated

def polygon_offset(polygon,offset,visualize = False):
    
    """
    Call to offset (shrink or expand) a polygon by a given distance
    @params:
        polygon     - Required  : list of x and y coordinates of closed polygon points
        offset      - Required  : offset in units; positive is outwards
        visualize   - Optional  : visualizes output, if desired
        
    Returns
        offset_polygon          : array of dimensions(polygon) with new vertices
    """
    
    import numpy as np
    import math
    
    global offset_polygon
    
    def perp(a) :
        b = np.empty_like(a)
        b[0] = -a[1]
        b[1] = a[0]
        return b
    
    def seg_intersect(a1,a2,b1,b2) :
        da = a2-a1
        db = b2-b1
        dp = a1-b1
        dap = perp(da)
        denom = np.dot( dap, db)
        num = np.dot( dap, dp )
        return (num / denom.astype(float))*db + b1
    
    # Prepare edges for direction determination
    edges = np.zeros(len(polygon[:,0])-1)
    for vert in range(len(polygon[:,0])-1):
         edges[vert] = (polygon[vert+1,0]-polygon[vert,0])*\
             (polygon[vert+1,1]+polygon[vert,1])
    
    # Determine direction
    if np.sum(edges) < 0:
        clockwise = False
    else:
        clockwise = True
        
    offset_edgepts = np.zeros((len(polygon[:,0])-1,2,2))    
    
    # Find and offset edges
    for edge in range(len(polygon[:,0])-1):
        
        # Normals are only based on relative distance, pivot shifts it to proper position
        pivot = np.asarray([(polygon[edge,0]+polygon[edge+1,0])/2,
                 (polygon[edge,1]+polygon[edge+1,1])/2])
        
        # Determine edge normal vector by points
        normal_pt_1 = [-(polygon[edge+1,1]-polygon[edge,1]),\
                       +(polygon[edge+1,0]-polygon[edge,0])]\
                       + pivot
        normal_pt_2 = [+(polygon[edge+1,1]-polygon[edge,1]),\
                       -(polygon[edge+1,0]-polygon[edge,0])]\
                       + pivot
        
        # Get the intersetction with the edge
        center = seg_intersect(np.asarray(polygon[edge,:]),np.asarray(polygon[edge+1,:]),np.asarray(normal_pt_1),np.asarray(normal_pt_2))
        
        if visualize == True:
            import matplotlib.pyplot as plt
            plt.scatter(center[0],center[1],c='r')
        
        # If points are clockwise, then normal_pt_2 is OUTSIDE, else INSIDE
        # Determine normal unit vector
        if clockwise == False:
            unit_vec = np.asarray([normal_pt_2[0]-center[0],normal_pt_2[1]-center[1]])/math.hypot(normal_pt_2[0]-center[0],normal_pt_2[1]-center[1])
        
        else:
            unit_vec = np.asarray([normal_pt_1[0]-center[0],normal_pt_1[1]-center[1]])/math.hypot(normal_pt_1[0]-center[0],normal_pt_1[1]-center[1])
                
        # Scale unit vector by offset
        unit_vec = unit_vec * offset
        
        # These are the offset edgepoints of the polygon
        offset_edgepts[edge,0,:] = polygon[edge,:] + unit_vec
        offset_edgepts[edge,1,:] = polygon[edge+1,:] + unit_vec
     
    # Pre-allocate offset polygon
    offset_polygon = polygon*0
        
    # Now find the new edge vertices
    for edge in range(len(offset_edgepts[:,0])):
        
        if edge != len(offset_edgepts[:,0])-1:
        
            offset_polygon[edge,:] = seg_intersect(
                    np.asarray(offset_edgepts[edge,0,:]),
                    np.asarray(offset_edgepts[edge,1,:]),
                    np.asarray(offset_edgepts[edge+1,0,:]),
                    np.asarray(offset_edgepts[edge+1,1,:]))
        else:
            offset_polygon[edge,:] = seg_intersect(
                    np.asarray(offset_edgepts[edge,0,:]),
                    np.asarray(offset_edgepts[edge,1,:]),
                    np.asarray(offset_edgepts[0,0,:]),
                    np.asarray(offset_edgepts[0,1,:]))
            offset_polygon[-1,:] = offset_polygon[0,:].copy()

    # Visualize output, if desired
    if visualize == True:
        plt.plot(polygon[:,0],polygon[:,1])
        plt.plot(offset_polygon[:,0],offset_polygon[:,1])
        plt.show()
    
    return offset_polygon

def hexagonal_tessellation(polygon, hexagon_seed, hexagon_radius, polygon_cutout = None, hexagon_orientation = 0, visualization = False):
    """
    Call to tessellate a given polygon with hexagons
    Updated for MODFLOW6
    
    @params:
        polygon             - Required  : coordinates of closed two-dimensional n-gon over which we want to tessellate, as array of dimensions [n,2]
        hexagon_seed        - Required  : point at which hexagon tessellation is seeded
        hexagon_radius      - Required  : radius of hexagons used for tessellation
        polygon_cutout      - Optional  : list of polygons within which no cells should be generated
        hexagon_orientation - Optional  : orientation of hexagon in clock-wise degrees [0 = flat top]
        visualization       - Optional  : boolean that defines whether a figure with the hexagons is printed
    
    Returns
        hexagon_gid_cores   : array of dimensions [ncel,2] specifying the x and y coordinates of each cell
        cg_xm               : correspondence grid, x in mesh format
        cg_ym               : correspondence grid, y in mesh format
        cg_xv               : correspondence grid, x in vector format
        cg_yv               : correspondence grid, y in vector format
        cg_ID               : correspondence grid, vector format IDs corresponding to cell IDs
        hex_area            : vector containing the surface area of each cell
        njag                : scalar number of all connections within the grid
        iac                 : vector indicating the number of connections plus 1 for each cell
        ja                  : list of cell number (n) followed by its connecting cell numbers (m) for each of the m cells connected to cell n

        
    """
    
    global hexagon_grid_cores
    
    # Import mplPath, required for assessing position in relation to polygon
    import math
    import numpy as np
    import matplotlib.path as mplPath
    import matplotlib.pyplot as plt
    
    # Define polygon rotation function
    def rotate(polygon,angle):
        """
        Rotates polygon clock-wise around angle
        
        @params:
            polygon             - Required  : coordinates of closed two-dimensional n-gon over which we want to tessellate, as array of dimensions [n,2]
            angle               - Required  : rotation angle around origin (0,0) in degrees
        """
        
        # Prepare variables
        angle = -math.radians(angle)
        polygon_rotated = []
        
        # Rotate each corner point
        for pt in polygon :
            polygon_rotated.append(( pt[0]*math.cos(angle)-pt[1]*math.sin(angle) , pt[0]*math.sin(angle)+pt[1]*math.cos(angle)) )
        
        polygon_rotated = np.asarray(polygon_rotated)
        
        # Return result
        return polygon_rotated
    
    # Return polygon for tessellation
    #   the specified angle is negative, as we are turning the polygon counter-clockwise instead to achieve the desired hexagon orientation
    polygon_rotated = rotate(polygon = polygon, angle = - hexagon_orientation)
    
    # Calculate path based on polygon outline
    polygonpath = mplPath.Path(polygon_rotated)
    
    if polygon_cutout is not None:
        # Initialize an empty list
        polygonpath_cutout = []
        # Add all polygon paths for the cutouts
        for poly in polygon_cutout:
            polygonpath_cutout.append(mplPath.Path(rotate(polygon = poly, angle = - hexagon_orientation)))
    
    # Calculate x and y increments
    increment_y = hexagon_radius
    increment_x = hexagon_radius*math.tan(math.radians(60))

    # Rotate hexagon_seed as well
    hexagon_seed_rotated = np.zeros(2)
    hexagon_seed_rotated[0] = hexagon_seed[0]*math.cos(math.radians(hexagon_orientation))-hexagon_seed[1]*math.sin(math.radians(hexagon_orientation))
    hexagon_seed_rotated[1] = hexagon_seed[0]*math.sin(math.radians(hexagon_orientation))+hexagon_seed[1]*math.cos(math.radians(hexagon_orientation))

    # Determine range over which points will be seeded, just considering points in rectangular orientation
    tessellation_range_x = [np.floor((np.min(polygon_rotated[:,0])-(hexagon_seed_rotated[0]))/increment_x/2), \
                                       np.ceil((np.max(polygon_rotated[:,0])-hexagon_seed_rotated[0])/increment_x/2)]
    tessellation_range_y = [np.floor((np.min(polygon_rotated[:,1])-(hexagon_seed_rotated[1]))/increment_y/2), \
                                       np.ceil((np.max(polygon_rotated[:,1])-hexagon_seed_rotated[1])/increment_y/2)]
    
    # create array for proposal points
    hexagon_points_x = []
    hexagon_points_y = []
    
    # add first set of in-phase polygons
    for idx_row in range(tessellation_range_y[0].astype(int),tessellation_range_y[1].astype(int)+1):
        for idx_col in range(tessellation_range_x[0].astype(int),tessellation_range_x[1].astype(int)+1):
            
            #==================================================================
            # In-Phase point
            #==================================================================
            
            # Proposed point
            proposal_point_x = hexagon_seed_rotated[0]+idx_col*increment_x*2
            proposal_point_y = hexagon_seed_rotated[1]+idx_row*increment_y*2
            
            if polygon_cutout is None:
            
                # If it's in the polygon, accept it
                if polygonpath.contains_point((proposal_point_x,proposal_point_y)) == True:
                    # Append coordinates
                    hexagon_points_x = np.append(hexagon_points_x,proposal_point_x)
                    hexagon_points_y = np.append(hexagon_points_y,proposal_point_y)
            
            else: # Also consider cutouts
                
                # If it's in the polygon, accept it
                if polygonpath.contains_point((proposal_point_x,proposal_point_y)) == True:
                    
                    # Go through all cutouts, check if the point lies within
                    passed_check = True
                    for polypath in polygonpath_cutout:
                        if polypath.contains_point((proposal_point_x,proposal_point_y)) == True: 
                            # Point is within one of the cutouts
                            passed_check = False
                    
                    # Only append the points if it does not lie within one of the cutouts
                    if passed_check:
                        # Append coordinates
                        hexagon_points_x = np.append(hexagon_points_x,proposal_point_x)
                        hexagon_points_y = np.append(hexagon_points_y,proposal_point_y)
                
            #==================================================================
            # Out-Of-Phase point (in positive x and positive y)
            #==================================================================
            
            # Proposed point
            proposal_point_x = hexagon_seed_rotated[0]+idx_col*increment_x*2+increment_x
            proposal_point_y = hexagon_seed_rotated[1]+idx_row*increment_y*2+increment_y
            
            if polygon_cutout is None:
            
                # If it's in the polygon, accept it
                if polygonpath.contains_point((proposal_point_x,proposal_point_y)) == True:
                    # Append coordinates
                    hexagon_points_x = np.append(hexagon_points_x,proposal_point_x)
                    hexagon_points_y = np.append(hexagon_points_y,proposal_point_y)
                
            else: # Also consider cutouts
                
                # If it's in the polygon, accept it
                if polygonpath.contains_point((proposal_point_x,proposal_point_y)) == True:
                    
                    # Go through all cutouts, check if the point lies within
                    passed_check = True
                    for polypath in polygonpath_cutout:
                        if polypath.contains_point((proposal_point_x,proposal_point_y)) == True: 
                            # Point is within one of the cutouts
                            passed_check = False
                    
                    # Only append the points if it does not lie within one of the cutouts
                    if passed_check:
                        # Append coordinates
                        hexagon_points_x = np.append(hexagon_points_x,proposal_point_x)
                        hexagon_points_y = np.append(hexagon_points_y,proposal_point_y)

    
    # Compile output variable
    hexagon_grid_cores = np.column_stack((hexagon_points_x, hexagon_points_y))
    
    #==========================================================================
    # Create Correspondence grid
    #==========================================================================
    
    # Sort grid cores in 0-rotation setting, required for correspondence grid
    import operator
    hexagon_grid_cores = sorted(np.around(hexagon_grid_cores,decimals = 6), key=operator.itemgetter(0))
    hexagon_grid_cores = sorted(np.around(hexagon_grid_cores,decimals = 6), key=operator.itemgetter(1), reverse=True)
    hexagon_grid_cores = np.asarray(hexagon_grid_cores, dtype = np.float64)
    
    # Grid spacing
    dy = hexagon_radius
    dx = hexagon_radius*math.tan(math.radians(60))
    # np.sqrt(4*hexagon_radius**2)
    
    # Grid boundaries
    xmin = np.min(hexagon_grid_cores[:,0])
    xmax = np.max(hexagon_grid_cores[:,0])
    ymin = np.min(hexagon_grid_cores[:,1])
    ymax = np.max(hexagon_grid_cores[:,1])
    
    # Grid range
    xrange = np.int(np.ceil((xmax-xmin)/dx))
    yrange = np.int(np.ceil((ymax-ymin)/dy))
    
    # Meshgrid dummies   
    X = np.arange(start = 0,stop = xrange)*dx + xmin
    Y = np.arange(start = 0,stop = yrange)*dy + ymin

    # Correspondence grid: mesh format
    cg_xm, cg_ym = np.meshgrid(X,Y)
    
    # Correspondence grid: vector format
    cg_xv = np.reshape(cg_xm, xrange*yrange)
    cg_yv = np.reshape(cg_ym, xrange*yrange)
    
    # Create rotated version of correspondence grid
    cg_rotated_back = rotate(polygon = np.column_stack((cg_xv,cg_yv)), angle = hexagon_orientation)
    cg_xv_r         = cg_rotated_back[:,0].copy()
    cg_yv_r         = cg_rotated_back[:,1].copy()
    cg_xm_r         = np.reshape(cg_xv_r,cg_xm.shape)
    cg_ym_r         = np.reshape(cg_yv_r,cg_ym.shape)
    
    #==========================================================================
    # Correspondence grid allocation finished, rotate system back
    #==========================================================================
    
    # Tessellation complete, rotate back
    #   hexagon_grid_cores isn't a polygon, but the function works nonetheless
    hexagon_grid_cores = rotate(polygon = hexagon_grid_cores, angle = hexagon_orientation)
    
    # Round the hexagon grid coordinates to snap them to regular positions
    hexagon_grid_cores = np.around(hexagon_grid_cores,decimals = 12)
    
    # Sort the rotated hexagon grid cores again; cg_ID remains true because the cg field doesn't care about rotation
    #   Sort polygons by rows, then columns
    import operator
    hexagon_grid_cores = sorted(np.around(hexagon_grid_cores,decimals = 6), key=operator.itemgetter(0))
    hexagon_grid_cores = sorted(np.around(hexagon_grid_cores,decimals = 6), key=operator.itemgetter(1), reverse=True)
    hexagon_grid_cores = np.asarray(hexagon_grid_cores, dtype = np.float64)
    
    #==========================================================================
    # Find which grid cores (in true rotation) relate to cg
    #==========================================================================
    
    # Concatenated dummy, create and rotate
    x_and_y = np.transpose(np.concatenate(([cg_xv],[cg_yv])))
    x_and_y = rotate(polygon = x_and_y, angle = hexagon_orientation)
    
    # Grid snapper
    def grid_snapper(node, nodes):
        nodes = np.asarray(nodes)
        dist_2 = np.sum((nodes - node)**2, axis=1)
        return np.argmin(dist_2)
    
    # Pre-allocate memory for cg_ID
    cg_ID = np.zeros(len(hexagon_grid_cores[:,0]),dtype = np.int)
    
    # Designate cg_ID by snapping hexagon_grid_cores to cg
    for i in range(len(hexagon_grid_cores[:,0])):
        cg_ID[i] = grid_snapper(node = hexagon_grid_cores[i,:],nodes = x_and_y)
        
    
    #==========================================================================
    # Determine flopy grid variables
    #==========================================================================
    
    import math
    import scipy.spatial
    
    # Determine flow area, assuming hexagons of fixed z-depth
    hexagon_surface_area = (math.tan(math.radians(30)) * hexagon_radius) * hexagon_radius * 6
    
    # MODFLOW USG cell connection variables
    hex_area = []   # surface area of hexagon
    iac = []        # number of connections (+1 for self) per cell
    ja = []         # cell IDs for those connections (negative for self)
    
    # Create KDTree for quick distance searching
    tree = scipy.spatial.KDTree(data        = hexagon_grid_cores,
                                leafsize    = 1000)
    
    ncel = len(hexagon_grid_cores[:,0])
    
    njag    = 0
    
    # Check for neighbours by going through list of polygons
    for mainhex in range(ncel):
    
        print('\rProgress:|%s| %s%% \r' % ('\033[36m'+'█' * int(50 * (mainhex+1) // ncel) + '-' * (50 - int(50 * (mainhex+1) // ncel))+'\033[0m', ("{0:." + str(1) + "f}").format(100 * ((mainhex+1) / float(ncel)))), end = '\r')
        
#        # Allocate space for distance variable
#        hex_connected = np.ones(len(hexagon_grid_cores[:,0]))*np.NAN
        
        njag += 1

        ja_dummy = [mainhex]

        # Query the tree for distance
        d,i = tree.query(hexagon_grid_cores[mainhex,:], 
                     k = 7) # We cannot have more than seven neighbours (including self)
        for subhex in range(len(d)):
            # The seven nearest cells aren't necessarily neighbouring
            if d[subhex] < hexagon_radius*2.1 and i[subhex] != mainhex:  # Also returns own index; 2.1 to ensure float precision doesn't mess up
                ja_dummy.append(i[subhex])
                njag += 1
                
        ja_dummy = [ja_dummy[0]] + list(np.sort(np.asarray(ja_dummy[1:],dtype=int)))
           
        ja.append(ja_dummy.copy())
#        ja.append(tuple(ja_dummy.copy()))
        
#        ja += ja_dummy.copy()
                
#                hex_connected[i[subhex]] = 0    # Mark with 0, all else remains NaN

        # Update iac
        iac.append(len(ja_dummy)) # This appends the number of connections (including self)
        
#        # Assign cell ID to connected hexagons
#        # COMMENT: Cell ID may not be 0, MODFLOW can't handle that
#        for subhex in range(len(hex_connected)):
#            if hex_connected[subhex] == 0:  # Save marked cells as connections
#                # Save cell ID
#                if subhex == mainhex:
#                    hex_connected[subhex] = -subhex-1   # save self ID as negative, is reverted by FloPy
#                else:
#                    hex_connected[subhex] = subhex+1    # Plus and minus one are necessary because FloPy counts from cell 1, not 0
#                    
#        
#        # Collapse and sort vector
#        hex_connected = np.asarray(sorted(hex_connected[np.isnan(hex_connected) == False]),dtype = np.int32)
#
#        # Update ja
#        ja = np.append(ja,hex_connected)
    
#    ja = ja.astype(int)
#    iac = iac.astype(int)
#    njag = len(ja)
    
    # Designate hex_area
    hex_area = np.ones(ncel) * hexagon_surface_area
        
    #==========================================================================
    # Finished cg_ID allocation
    #==========================================================================
    
    # Print if visualization is requested
    if visualization == True:
        
        color_orange = [253/255,102/255,0]
        color_blue = [0,120/255,255/255]
        
        fig, ax = plt.subplots()
        
        # Plot it for visualization
        plt.plot(polygon[:,0],polygon[:,1],color = color_orange)
        
        # Hexagon radius only goes to normal of sides
        edgepoint_distance = hexagon_radius/np.cos(np.deg2rad(30))
        
        # Plot hexagons
        for hex in range(len(hexagon_grid_cores[:,0])):
            x = []
            y = []
            for angle in range(0-hexagon_orientation, 420-hexagon_orientation, 60):
                x = np.append(x,hexagon_grid_cores[hex,0]+math.cos(math.radians(angle)) * edgepoint_distance)
                y = np.append(y,hexagon_grid_cores[hex,1]+math.sin(math.radians(angle)) * edgepoint_distance)
            plt.plot(x, y, color = color_blue)
            
        # Plot hexagon seed
        x = []
        y = []
        for angle in range(0-hexagon_orientation, 420-hexagon_orientation, 60):
            x = np.append(x,hexagon_seed[0]+math.cos(math.radians(angle)) * edgepoint_distance)
            y = np.append(y,hexagon_seed[1]+math.sin(math.radians(angle)) * edgepoint_distance)
        plt.plot(x, y, linewidth=2, color = color_orange)
           
        scatterlabel=list(range(len(hexagon_grid_cores[:,0])))

        
        for i, txt in enumerate(scatterlabel):
            ax.annotate(txt, (hexagon_grid_cores[i,:]), ha='center', va='center', color = color_blue)
        
        plt.axis('scaled')
        
        
    return hexagon_grid_cores, cg_xm, cg_ym, cg_xv, cg_yv, cg_ID, cg_xm_r, cg_ym_r, cg_xv_r, cg_yv_r, hex_area, njag, iac, ja

#%%

def depth_specific_usgfiles(hexagon_grid_cores,hexagon_radius,hexagon_depth,hex_area, njag, iac, ja, guarantee_connection = None):
    
    """
    Call to build MODFLOW USG connection files for hexagons, based on previous hexagon tessellation
    
    @params:
        hexagon_grid_cores  - Required  : tessellated polygons over area of interest
        hexagon_radius      - Required  : radius of hexagons used for tessellation
        hexagon_depth       - Required  : depth of hexagons, required for calculation of flow-perpendicular area; of structure [nx2], designated as [[zbot,top],...]
        hex_area            - Required  : vector containing the surface area of each cell
        njag                - Required  : scalar number of all connections within the grid
        iac                 - Required  : vector indicating the number of connections plus 1 for each cell
        ja                  - Required  : list of cell number (n) followed by its connecting cell numbers (m) for each of the m cells connected to cell n
        guarantee_connection- Optional  : can be a float with a z height, guarantees that all neighbouring cells are connected at least by this

    Returns:
        cl12        : vector containing distance to grid interface for each cell
        fahl        : vector containing the interface area between two cells for each connection
    """  
    
    import numpy as np
    import math
    
    ncel = len(hexagon_grid_cores[:,0])
    
    fahl = []
    cl12 = []
    
    hexagon_side = math.tan(math.radians(30))*2* hexagon_radius
    
    # Go through all connections
    for n in range(ncel):
        
        dummy_fahl = []
        dummy_cl12 = []
        
        for m,con in enumerate(ja[n]):
            
            if m == 0:
                dummy_fahl.append(0.)
                dummy_cl12.append(0.)
            else:
                interface_height = max(0, \
                                   min(hexagon_depth[n,1], \
                                       hexagon_depth[con,1]) - \
                                   max(hexagon_depth[n,0], \
                                       hexagon_depth[con,0]))
                                   
                # Guarantee a minimum connection, if desired
                if interface_height == 0 and guarantee_connection is not None:
                    interface_height = guarantee_connection
                                       
                dummy_fahl.append(interface_height*hexagon_side)
                dummy_cl12.append(hexagon_radius)
                
        fahl.append(dummy_fahl.copy())
        cl12.append(dummy_cl12.copy())
        
#        if ja[idx] < 0:
#            
#            # Extract last negative index
#            ln_idx      = -ja[idx]
#            
#            fahl[idx]   = 0
#            cl12[idx]   = 0
#            
#        else:
#            
##            interface_height = max(0, \
##                                   min(hexagon_depth[ln_idx-1,1], \
##                                       hexagon_depth[ja[idx]-1,1]) - \
##                                   max(hexagon_depth[ln_idx-1,0], \
##                                       hexagon_depth[ja[idx]-1,0]))               
#            interface_height = max(0, \
#                                   min(hexagon_depth[ln_idx-1,1], \
#                                       hexagon_depth[ja[idx]-1,1]) - \
#                                   max(hexagon_depth[ln_idx-1,0], \
#                                       hexagon_depth[ja[idx]-1,0])) 
#                                   
#                   
#            # Guarantee a minimum connection, if desired
#            if interface_height == 0 and guarantee_connection is not None:
#                interface_height = guarantee_connection
#                                   
#            fahl[idx]   = interface_height*hexagon_side
#            cl12[idx]   = hexagon_radius
    
    # Return results
    return cl12,fahl


#%%
    
#def extract_ztop(hexagon_grid_cores,filename):
#    
#    """
#    This function extracts the cell top elevations from a DEM.
#    
#    @params:
#        hexagon_grid_cores  - Required  : tessellated polygons over area of interest
#        filename            - Required  : list of cell number (n) followed by its connecting cell numbers (m) for each of the m cells connected to cell n
#
#    Returns:
#        ztop        : vector containing distance to grid interface for each cell
#    """  
#    
#    import numpy as np
#    import pickle
##    import PIL
#    
#    
##    img = PIL.Image.open("ztop_smoothed.tiff").convert("L")
##    imgarr = np.array(img)
#    
#    dictionary = pickle.load(open(filename,'rb'))
#    
#    ztop = np.zeros(len(hexagon_grid_cores[:,0]))
#    
#    dem = dictionary['dem']
#    X = dictionary['X'][0,:]
#    Y = dictionary['Y'][:,0]
#    
#    for cel in range(len(hexagon_grid_cores[:,0])):
#        
#        dist_x = abs(hexagon_grid_cores[cel,0]-X)
#        idx_x = np.where(dist_x == np.min(dist_x))[0]
#        if len(idx_x) > 0: # If two points have the same distance, use the first
#            idx_x = idx_x[0]
#            
#        dist_y = abs(hexagon_grid_cores[cel,1]-Y)
#        idx_y = np.where(dist_y == np.min(dist_y))[0]
#        if len(idx_y) > 0: # If two points have the same distance, use the first
#            idx_y = idx_y[0]
#            
#        ztop[cel] = dem[idx_y,idx_x]
#        
##    imgarr /= 255
##    imgarr *= (598.3332-513.6999)
##    imgarr += 513.6999
#    
#    return ztop

def extract_ztop(hexagon_grid_cores,filename,hexagon_radius,array=None):
    
    """
    This function extracts the cell top elevations from a DEM.
    
    @params:
        hexagon_grid_cores  - Required  : tessellated polygons over area of interest
        filename            - Required  : list of cell number (n) followed by its connecting cell numbers (m) for each of the m cells connected to cell n
        hexagon_radius      - Required  : radius of the hexagonal cells, in m
        array               - Optional  : array with arbitrary entries

    Returns:
        ztop        : vector containing distance to grid interface for each cell
    """  
    
    import numpy as np
    import pickle
    
    dictionary = pickle.load(open(filename,'rb'))
    
    ztop = np.zeros(len(hexagon_grid_cores[:,0]))
    
    if array is None:
        dem = dictionary['dem']
    else:
        dem = array.copy()
#        dem = np.asarray(dem,dtype=float)
#        dem = (dem-np.min(dem))/(np.max(dem)-np.min(dem))
#        dem = dem*(np.max(dictionary['dem'])-np.min(dictionary['dem']))+np.min(dictionary['dem'])
        
    X = dictionary['X'][0,:]
    Y = dictionary['Y'][:,0]
    
    dx = np.diff(X)[0]
    
#    # Should be ceil, but we subtract one entry for center, so floor
#    cell_buffer = int(np.ceil(hexagon_radius/dx))
    
    indices = []
    
    for cel in range(len(hexagon_grid_cores[:,0])):
        
        dist_x = abs(hexagon_grid_cores[cel,0]-X)
        idx_x = np.where(dist_x == np.min(dist_x))[0]
        if len(idx_x) > 0: # If two points have the same distance, use the first
            idx_x = idx_x[0]
            
        dist_y = abs(hexagon_grid_cores[cel,1]-Y)
        idx_y = np.where(dist_y == np.min(dist_y))[0]
        if len(idx_y) > 0: # If two points have the same distance, use the first
            idx_y = idx_y[0]
            
        # idx_x and idx_y give us the center of the cell in the DEM array
#        elevations  = dem[idx_y-cell_buffer:idx_y+cell_buffer,idx_x-cell_buffer:idx_x+cell_buffer]
#        ztop[cel]   = np.max(elevations)
        
        ztop[cel]   = dem[idx_y,idx_x]
        
        indices.append([idx_y,idx_x])
    
    return ztop, indices


#%%


def refine_hexagons(hexagon_grid_cores,hexagon_radius, hexagon_depth, hex_area,
                    iac ,ja ,cl12, fahl, refine_cells, refine_levels,ztop,zbot):
    
    import math
    import numpy as np
    import copy
    
    """
    Refine cells counts, does not work for directly neighbouring cells.
    
    
    """
    
    dict_refine_cells = {}
    
    for hx in range(len(refine_cells)):
        
#        refine_cells_inner = {}
#        
#        refine_cells_inner[refine_cells[hx]]  = 
#        refine_cells_inner['neighbours']  = []
        
        cel = refine_cells[hx]+1 # Flopy here requires cells to count from 1, for whatever reason
        
        ja_idx = np.where(ja == -cel)[0]
        cons    = int(iac[cel]-1)
        
        if cons < 6:
            raise Exception('Refinement stopped - for now, it is only implemented for a full neighbourhood.')
        
        connections = ja[np.arange(ja_idx+1,ja_idx+cons+1,1,dtype=int)] # These are the cell indices +1
        
        # Now we have to sever the connections between the old hexagons
        # First, we find where the connections are defined
        con_dix = []
        for c1 in np.arange(ja_idx+1,ja_idx+cons+1,1,dtype=int):
            
            # Each c1 is where the other cell connects to hx
            
            # Now find index of second cell
            idx = np.where(ja == -ja[c1])[0][0]
            counter = 1
            continue_flag = True
            while continue_flag:
                # Write new value until stop
                if ja[idx+counter] > 0:
                    # Only save it if it connects to hx
                    if ja[idx+counter] == cel:
                        con_dix.append([c1,idx+counter])
                    counter += 1
                else:
                    # New cell in ja? Stop looking.
                    continue_flag = False
        
        # Learn how the neighbours are connected
        # Connection list knows, for each neighbour, to which two others it connects
        connection_list = []

        
        for c in connections:
            dummy = []
            idx = np.where(ja == -c)[0]
            counter = 1
            continue_flag = True
            while continue_flag:
                # Write new value until stop
                if ja[idx+counter] > 0:
                    # Only save it if it connects to hx
                    if ja[idx+counter] in connections:
                        dummy.append(ja[idx+counter][0])
                    counter += 1
                else:
                    # New cell in ja? Stop looking.
                    continue_flag = False
            connection_list.append(dummy.copy()) 
        
        # Re-write connections to indexadditions
        connection_list_relabelled = copy.deepcopy(connection_list)
        
        for i,entr1 in enumerate(connection_list_relabelled):
            for j,entr2 in enumerate(entr1):
                connection_list_relabelled[i][j] = np.where(connections == entr2)[0][0]
        
        # Go through all desired levels
        maxcel = -np.min(ja)
        newcells_idx = np.arange(6*refine_levels[hx])+1
        newcells_idx = -newcells_idx - maxcel # These are the new indices, negative for self-reference
        
        jastrings = []
        cl12strings = []
        fahlstrings = []
        
        cl12center = np.max(cl12)/(refine_levels[hx]*2+1)

        internal_height = ztop[refine_cells[hx]]-zbot[refine_cells[hx]]
        
        hex_area[refine_cells[hx]] = 6/2*cl12center**2
        
        fahllateral = (hexagon_radius*refine_levels[hx]*2/(refine_levels[hx]*2+1) - hexagon_radius*(refine_levels[hx]-1)*2/(refine_levels[hx]*2+1))*internal_height
        
        # Now alter the connections piece by piece
        for lvl in range(refine_levels[hx]):
            
            if lvl == 0: # We are in the first level - connect outwards
                
                outercons = connections.copy()
                
                cl12lateral = hexagon_radius*refine_levels[hx]*2/(refine_levels[hx]*2+1)
                
                fahlinwards = internal_height * math.tan(math.radians(30))  * 2 \
                    * hexagon_radius*refine_levels[hx]*2/(refine_levels[hx]*2+1)
                
                
                dummy_area = \
                    0.5*(hexagon_radius*1)**2 - \
                    0.5*(hexagon_radius*(refine_levels[hx]-1)*2/(refine_levels[hx]*2+1))**2
                
                hex_area = np.concatenate((hex_area,np.ones(6)*dummy_area))
                
                
                for nc in range(6): # Go through all six new cells
                    
                    iac = np.concatenate((iac,np.ones(1)*5)) # Append new connections
                    
                    ix =lvl*6+nc
                    
                    # Pick the old connection to center and rewire it to the new cell
                    ja[con_dix[nc][1]] = -newcells_idx[ix]
                    cl12[con_dix[nc][0]] = cl12center # change cl12 of center one
                    
                    # Write a dummy list of connections for jastring 
                    dummy = [newcells_idx[ix]]
                    dummy.append(outercons[nc]) # Index of outer cell
                    for entry in connection_list_relabelled[nc]:
                        dummy.append(-newcells_idx[ix + entry-nc]) # Append neighbour on same level
                    # Connect inner cell
                    if refine_levels[hx] > 1: # If there is a second refinement level, connect to that
                        dummy.append(-newcells_idx[nc+6])
                    else: # Else, connect to core
                        dummy.append(cel)
                        # Also change core reference
                        fahl[con_dix[nc][0]] = fahlinwards # Also has to be changed
                        ja[con_dix[nc][0]] = -newcells_idx[ix] 
                    jastrings.append(copy.deepcopy(dummy))
                        

                    # Write a dummy list of connections for cl12string 
                    dummy = [0]
                    dummy.append(cl12center) # Cl12 to outer cell
                    for entry in connection_list_relabelled[nc]:
                        dummy.append(cl12lateral) # Append neighbour on same level
                    dummy.append(cl12center) # Cl12 to inner cell
                    cl12strings.append(copy.deepcopy(dummy))
                    
                    
                    # Write a dummy list of fahl entries
                    dummy = [0] # Self is always zero
                    dummy.append(fahl[con_dix[nc][1]]) # Cl12 to outer cell
                    for entry in connection_list_relabelled[nc]:
                        dummy.append(fahllateral) # Connect to all neighbours
                    if refine_levels[hx] > 1: # If there is a second refinement level, connect to that
                        dummy.append(fahlinwards)
                    else: # Else, connect to core
                        dummy.append(fahlinwards)
                        # Also change core reference
                        fahl[con_dix[nc][0]] = fahlinwards
                    fahlstrings.append(copy.deepcopy(dummy))
                    

                    outercons[nc] = -newcells_idx[ix]

            else: # We are in an outer layer
                
                cl12lateral = hexagon_radius*(refine_levels[hx]-lvl)*2/(refine_levels[hx]*2+1)
                
                fahlinwards = internal_height * math.tan(math.radians(30))  * 2 \
                    * hexagon_radius*(refine_levels[hx]-lvl)*2/(refine_levels[hx]*2+1)
                    
                dummy_area = \
                    0.5*(hexagon_radius*(refine_levels[hx]-lvl+1)*2/(refine_levels[hx]*2+1))**2 - \
                    0.5*(hexagon_radius*(refine_levels[hx]-lvl)*2/(refine_levels[hx]*2+1))**2
                
                hex_area = np.concatenate((hex_area,np.ones(6)*dummy_area))
                
                for nc in range(6): # Go through all six new cells
                    
                    iac = np.concatenate((iac,np.ones(1)*5)) # Append new connections
                    
                    ix =lvl*6+nc
                    
                    # Write a dummy list of connections for jastring 
                    dummy = [newcells_idx[ix]]
                    dummy.append(outercons[nc]) # Index of outer cell
                    for entry in connection_list_relabelled[nc]:
                        dummy.append(-newcells_idx[ix + entry-nc]) # Append neighbour on same level lvl*6 + entry-nc
                    # Connect inner cell
                    if refine_levels[hx] > lvl+1: # If there is a further refinement level, connect to that
                        dummy.append(-newcells_idx[ix+6])
                    else: # Else, connect to core
                        dummy.append(cel)
                        # Also change core reference
                        fahl[con_dix[nc][0]] = fahlinwards # Also has to be changed
                        ja[con_dix[nc][0]] = -newcells_idx[ix] 
                    jastrings.append(copy.deepcopy(dummy))
                    
                    # Write a dummy list of connections for cl12string 
                    dummy = [0]
                    dummy.append(cl12center) # Cl12 to outer cell
                    for entry in connection_list_relabelled[nc]:
                        dummy.append(cl12lateral) # Append neighbour on same level
                    dummy.append(cl12center) # Cl12 to inner cell
                    cl12strings.append(copy.deepcopy(dummy))
                    
                    # Write a dummy list of fahl entries
                    dummy = [0] # Self is always zero
                    dummy.append(fahl[con_dix[nc][1]]) # Cl12 to outer cell
                    for entry in connection_list_relabelled[nc]:
                        dummy.append(fahllateral) # Connect to all neighbours
                    if refine_levels[hx] > 1: # If there is a second refinement level, connect to that
                        dummy.append(fahlinwards)
                    else: # Else, connect to core
                        dummy.append(fahlinwards)
                        # Also change core reference
                        fahl[con_dix[nc][0]] = fahlinwards
                    fahlstrings.append(copy.deepcopy(dummy))
                    
                    outercons[nc] = -newcells_idx[ix]
                
        # Expand vectors
        for entry in jastrings:
            ja = np.concatenate((ja,np.asarray(entry)))
        for entry in cl12strings:
            cl12 = np.concatenate((cl12,np.asarray(entry)))  
        for entry in fahlstrings:
            fahl = np.concatenate((fahl,np.asarray(entry)))     
            
    # Number of connections in an unstructured grid
    njag = len(ja)
                
    return (hex_area,njag,iac,ja,cl12,fahl)
                

def hexagon_usgfiles (hexagon_grid_cores,hexagon_radius, hexagon_depth):
    
    """
    Call to build MODFLOW USG connection files for hexagons, based on previous hexagon tessellation
    
    @params:
        hexagon_grid_cores  - Required  : tessellated polygons over area of interest
        hexagon_radius      - Required  : radius of hexagons used for tessellation
        hexagon_depth       - Required  : depth of hexagons, required for calculation of flow-perpendicular area; of structure [nx2], designated as [[zbot,top],...]
        
    Returns:
        hex_area    : vector containing the surface area of each cell
        njag        : scalar number of all connections within the grid
        iac         : vector indicating the number of connections plus 1 for each cell
        ja          : list of cell number (n) followed by its connecting cell numbers (m) for each of the m cells connected to cell n
        cl12        : vector containing distance to grid interface for each cell
        fahl        : vector containing the interface area between two cells for each connection
    """  
    
    global hex_area, iac, ja, cl12, fahl
    
    import math
    import numpy as np
    import scipy.spatial
    
    # Determine flow area, assuming hexagons of fixed z-depth
    hexagon_flow_area = math.tan(math.radians(30)) * 2 * hexagon_radius * hexagon_depth
    hexagon_surface_area = (math.tan(math.radians(30)) * hexagon_radius) * hexagon_radius * 6
    
    # MODFLOW USG cell connection variables
    hex_area = []   # surface area of hexagon
    iac = []        # number of connections (+1 for self) per cell
    ja = []         # cell IDs for those connections (negative for self)
    cl12 = []       # distance of cell center to respective flow face
    fahl = []       # flow area for each connection between hexagons (0 for self)
    
    # Create KDTree for quick distance searching
    tree = scipy.spatial.KDTree(data        = hexagon_grid_cores,
                                leafsize    = 1000)
    
    ncel = len(hexagon_grid_cores[:,0])
    
    # Check for neighbours by going through list of polygons
    for mainhex in range(ncel):
    
        print('\rProgress:|%s| %s%% \r' % ('\033[36m'+'█' * int(50 * (mainhex+1) // ncel) + '-' * (50 - int(50 * (mainhex+1) // ncel))+'\033[0m', ("{0:." + str(1) + "f}").format(100 * ((mainhex+1) / float(ncel)))), end = '\r')
        
        # Allocate space for distance variable
        hex_connected = np.ones(len(hexagon_grid_cores[:,0]))*np.NAN
        # Query the tree for distance
        d,i = tree.query(hexagon_grid_cores[mainhex,:], 
                     k = 7) # We cannot have more than seven neighbours (including self)
        for subhex in range(len(d)):
            # The seven nearest cells aren't necessarily neighbouring
            if d[subhex] < hexagon_radius*2.1:  # Also returns own index; 2.1 to ensure float precision doesn't mess up
                hex_connected[i[subhex]] = 0    # Mark with 0, all else remains NaN
        
#        # Go through each hexagon, calculate distance to current hex
#        for subhex in range(len(hexagon_grid_cores[:,0])):
#            
#            hex_distance = np.sqrt((hexagon_grid_cores[mainhex,0]-hexagon_grid_cores[subhex,0])**2 + \
#                (hexagon_grid_cores[mainhex,1]-hexagon_grid_cores[subhex,1])**2)
#            
#            # Find all connected hexagons
#            if hex_distance < hexagon_radius*2.1:   # Also returns own index; 2.1 to ensure float precision doesn't mess up
#                hex_connected[subhex] = 0           # Mark with 0, all else remains NaN

        # Update iac
        iac = np.append(iac,sum(i == 0 for i in hex_connected)) # This appends the number of connections (including self)
        
        # Assign cell ID to connected hexagons
        # COMMENT: Cell ID may not be 0, MODFLOW can't handle that
        for subhex in range(len(hex_connected)):
            if hex_connected[subhex] == 0:  # Save marked cells as connections
                # Save cell ID
                if subhex == mainhex:
                    hex_connected[subhex] = -subhex-1   # save self ID as negative, is reverted by FloPy
                else:
                    hex_connected[subhex] = subhex+1    # Plus and minus one are necessary because FloPy counts from cell 1, not 0
                    
        
        # Collapse and sort vector
        hex_connected = np.asarray(sorted(hex_connected[np.isnan(hex_connected) == False]),dtype = np.int32)
        
        # Determine interface area, use information of connections
        for subhex in range(len(hex_connected)):
            if subhex == 0: # self
                fahl_addition = 0   # area to itself is zero
            else:
                # find interface between two cells
                index = hex_connected[subhex] - 1   # revert cell index to row
                # Here we check if there is a flow area between:
                #   We compare the top elevation for hex and subhex, choosing the lower one
                #   We compare the bot elevation for hex and subjex, choosing the larger one
                #   The difference gives us the height; if they do not overlap, set to zero
                interface_height = max(0, \
                                       min(hexagon_depth[mainhex,1], \
                                           hexagon_depth[index,1]) - \
                                       max(hexagon_depth[mainhex,0], \
                                           hexagon_depth[index,0]))
                fahl_addition = np.append(fahl_addition,interface_height*\
                                          math.tan(math.radians(30)) \
                                          * 2 * hexagon_radius)
        
        # Update ja
        ja = np.append(ja,hex_connected)
        
        # Update fahl
#        fahl_addition = hex_connected.copy().astype(np.float64)
#        fahl_addition[fahl_addition >= 0] = hexagon_flow_area
#        fahl_addition[fahl_addition < 0] = 0
        fahl = np.append(fahl,fahl_addition)
    
    
    # Fahl cannot deal with cell 0
    fahl[0] = 0 # first value is self, TEMPORARY
    
    # Designate cl12
    cl12 = fahl.copy()
    cl12[cl12 != 0] = hexagon_radius
    
    # Designate hex_area
    hex_area = np.ones(len(hexagon_grid_cores[:,0])) * hexagon_surface_area
    
    # Number of connections in an unstructured grid
    njag = len(ja)
        
    # Enter new line of progress bar
    print()
    
    # Return results
    return (hex_area,njag,iac,ja,cl12,fahl)


    
def mask_value_over_polygon(grid,value,polygon):
    
    """
    Projects a specified value on all grid points that lie within a polygon
    
    @params:
        grid                - Required  : grid on which value will be projected, of dimensions [ncell,2]
        value               - Required  : value to be written inside the polygon
        polygon             - Required  : extent of polygon, of dimensions [ncell,2]
    
    Returns:
        poly_mask   : vector containing [value] for each grid cell within the polygon, 0 everywhere else
    """  
    
    import numpy as np
    import matplotlib.path as mplPath
    
    # Create polygon path
    bbPath = mplPath.Path(np.array(polygon))
    
    # Designate returned variable
    global poly_mask
    poly_mask = np.zeros(len(grid[:,0]))
    
    # Flag points that are inside this simplex
    for gridpts in range(len(poly_mask)):
        if bbPath.contains_point((grid[gridpts,0],grid[gridpts,1])) == True:
            poly_mask[gridpts] = value
            
    # Return result
    return poly_mask
    
def gaussian_smoother(grid,values,sd,verbose = False):
    
    """
    Projects a specified value on all grid points that lie within a polygon
    
    @params:
        grid        - Required  : grid on which value will be projected, of dimensions [ncell,2]
        values      - Required  : vector of values corresponding to grid's cells
        sd          - Required  : standard deviation of gaussian smoother
        verbose     - Optional  : boolean to determine whether a progress bar shall be printed
    
    Returns:
        values      : vector containing [value] for each grid cell within the polygon, 0 everywhere else
    """  
    
    def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', color = '\033[0m'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s|%s| %s%% %s' % (prefix, color+bar+'\033[0m', percent, suffix), end = '\r')
        # Print New Line on Complete
        if iteration == total: 
            print()
    
    import numpy as np
    import scipy
    
    # Copy original values into a safe variable
    values_prior    = values.copy()
    
    # Reset values variable
    values          = np.zeros(len(values))
    
    # Designate smoothing weight
    for cell in range(len(values)):
        
        if verbose == True:
            printProgressBar(cell, len(values)-1, prefix = 'Smoothing:', suffix = 'Complete', length = 50, color = '\033[35m')
        
        # Pre-allocate space
        weight      = np.zeros(len(values))
        
        # Determine distances
        for othercell in range(len(values)):
            
            distance = np.sqrt(\
                (grid[cell,0]-\
                 grid[othercell,0])       **2 +\
                (grid[cell,1]-\
                 grid[othercell,1])       **2)
                
            # Calculate and normalize weight
            weight[othercell] = scipy.stats.norm(0, sd).pdf(distance)
            
        # Normalize weights
        weight = weight/np.sum(weight)
        
        # Smooth the distribution
        for othercell in range(len(values)):
            values[cell] = values[cell] + values_prior[othercell]*weight[othercell]
            
    # Return result
    return values
        
def hexagon_density(grid,searchradius,searchgrid,normalize = True):
    
    """
    Projects a specified value on all grid points that lie within a polygon
    
    @params:
        grid            - Required  : grid on which value will be projected, of dimensions [ncell,2], designating centers of hexagons
        searchradius    - Required  : radius around which neighbours are sought
        searchgrid      - Required  : grid on which neighbours are sought, normally equal to grid
        normalize       - Optional  : boolean whether density is normalized, default is True
    
    Returns:
        hexdens                     : vector containing [density] for each grid cell 
    """
    
    import numpy as np
    
    global hexdens
    
    # Determine number of cells
    n_cells = len(grid[:,0])
    n_cells_sg = len(searchgrid[:,0])
    
    # Pre-allocate memory
    hexdens = np.zeros(n_cells)
    
    # Search for neighbours
    for cell in range(n_cells):
        
        # Go through each othercell in the grid, relative to cell
        for othercell in range(n_cells_sg):
            
            # Calculate distance
            distance = np.sqrt(\
                (grid[cell,0]-\
                 searchgrid[othercell,0])       **2 +\
                (grid[cell,1]-\
                 searchgrid[othercell,1])       **2)
            
            # Update density counter, if in range
            if distance <= searchradius:
                hexdens[cell] += 1
            
    # Normalize, if requested
    if normalize == True:
        
        # Set minimum to zero
        hexdens = hexdens - np.min(hexdens)
        
        # Check if density is uniform
        if np.sum(hexdens) == 0:
            raise Exception('Normalization failed: hexagon density is uniform. Consider a smaller searchradius.')
        
        # Normalize maximum to one
        hexdens = hexdens / np.max(hexdens)
        
    # Return result
    return hexdens

#%%
    
def generate_hexagon_vertex_list(hexagon_grid_cores,hexagon_radius,hexagon_orientation):
    
    """
    This function takes a hexagonal grid as input and returns a vertex representation
    
    
    """
    
    
    import math
    import numpy as np
    
    def rotate(polygon,angle):
        """
        Rotates polygon clock-wise around angle
        
        @params:
            polygon             - Required  : coordinates of closed two-dimensional n-gon over which we want to tessellate, as array of dimensions [n,2]
            angle               - Required  : rotation angle around origin (0,0) in degrees
        """
        
        # Reshape if one-dimensional
        if len(polygon.shape) == 1:
            polygon = polygon.reshape((1,polygon.shape[0]))
        
        # Prepare variables
        angle = -math.radians(angle)
        polygon_rotated = []
        
        # Rotate each corner point
        for pt in range(len(polygon[:,0])) :
            polygon_rotated.append(( polygon[pt,0]*math.cos(angle)-polygon[pt,1]*math.sin(angle) , polygon[pt,0]*math.sin(angle)+polygon[pt,1]*math.cos(angle)) )
        
        polygon_rotated = np.asarray(polygon_rotated)
        
        # Return result
        return polygon_rotated
    
    ncel = len(hexagon_grid_cores[:,0])
    
    # Define regular hexagon offset
    vertex_distance = hexagon_radius/math.cos(math.radians(30))
    
    # Define unrotated hexagon offset vertices
    hexagon_offset_vertices = np.zeros((6,2))
    hexagon_offset_vertices[:,1] = vertex_distance
    for idx,angle in enumerate(range(6)):
        hexagon_offset_vertices[idx,:] = rotate(hexagon_offset_vertices[idx,:],30+angle*60)
        
    # Rotate hexagon offset vertices by prescribed clockwise rotation (in degrees)
    hexagon_offset_vertices = rotate(hexagon_offset_vertices,hexagon_orientation)
    
    # Define grid by vertices
    vertex_list = []
    hexagon_vertices = []
    for cel in range(ncel):
                 
        print('\rGenerate:|%s| %s%% \r' % ('█' * int(50 * (cel+1) // ncel) + '-' * (50 - int(50 * (cel+1) // ncel)), ("{0:." + str(1) + "f}").format(100 * ((cel+1) / float(ncel)))), end = '\r')
        
        dummy_vertices = hexagon_grid_cores[cel,:] + hexagon_offset_vertices
        dummy_vertex_list = []
        
        if cel == 0:
            
            for hx in range(6):
                hexagon_vertices.append(dummy_vertices[hx,:].copy())
                dummy_vertex_list.append(hx)
            
        else:
            
            innovator = np.ones(6,dtype=bool)
            for hx in range(6):
#                
#                if not any((dummy_vertices[hx,:] == x).all() for x in hexagon_vertices):
                if all(abs(np.sum(dummy_vertices[hx,:]-np.asarray(hexagon_vertices),axis=1)) > 1E-3):
                    
                    dummy_vertex_list.append(len(hexagon_vertices))
                    hexagon_vertices.append(dummy_vertices[hx,:].copy())
                    
                else:
                    

                    dummy_vertex_list.append(
                    np.where(
                        abs(np.sum(dummy_vertices[hx,:]-np.asarray(hexagon_vertices),axis=1)) < 1E-3)[0][0])
                    
        vertex_list.append(dummy_vertex_list.copy())
        
    return vertex_list, hexagon_vertices
    