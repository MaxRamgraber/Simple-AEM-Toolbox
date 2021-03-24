# =============================================================================
"""
Max' utilities

This file contains a number of useful functions I have written over the years. 
I endeavour to keep it as much up-to-date as possible, but make no promises =P
"""
# =============================================================================
#%%
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '█', oae = False):
    """
    This function prints and updates a progress bar in the console. Mind that 
    printing other ouputs while this progress bar is not advised - these outputs
    will overwrite (and be overwritten by) this function.
    
    Parameters:
        
        iteration       - Required  : current iteration (Int)
        total           - Required  : total iterations (Int)
        prefix          - Optional  : prefix string (Str)
        suffix          - Optional  : suffix string (Str)
        decimals        - Optional  : positive number of decimals in percent complete (Int)
        length          - Optional  : character length of bar (Int)
        fill            - Optional  : bar fill character (Str)
        oae             - Optional  : overwrite after end (Bool)
    """
    
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total and not oae: 
        print()
        
#%%
def multivariate_plot(data,variable_names = None,variable_labels=None,axis_limits=None,**kwargs):
    
    """
    This function plots multivariate data as a collection of 2-D plots.
    
    Parameters:
        
        data            - Required  : a 2-D array of dimenions (number datapoints x dimensions parameter space)
                                      or a 3-D array of dimensions (number datapoints x dimensions parameter space x time)
        variable_names  - Optional  : list of strings specifying the names of the variables
        variable_labels - Optional  : list of strings specifying the axis labels for the variables
        axis_limits     - Optional  : list of 2-element-lists (for lower and upper specified axis limits) or None 
                                      (for unspecified axis limits)
    """
    
    import matplotlib.pyplot as plt
    
    # Find data dimensions
    N = data.shape[0]
    D = data.shape[1]
    if len(data.shape) == 3: # A time dimension is specified
        T = data.shape[2]
    else: # No time dimension is specified
        T = None
        
    # If variable names or labels are specified, check if dimensions match
    text = ''
    if variable_names is not None:
        if len(variable_names) != D:
            text += '\nNumber of variable name entries does not match data dimensions!'
    if variable_labels is not None:
        if len(variable_labels) != D:
            text += '\nNumber of variable axis label entries does not match data dimensions!'
    if axis_limits is not None:
        if len(axis_limits) != D:
            text += '\nNumber of variable axis limit entries does not match data dimensions!'
    if len(text) > 0:
        raise Exception(text)
    
    # Now plot everything
    if T is None: # If data array is two-dimensional (not a time series)
        for row in range(D):
            for col in range(D):
                plt.subplot(D,D,row*D+col+1)
                if row == col and variable_names is not None: # If this point lies along the diagonal
                    plt.text(0.5, 0.5, variable_names[row], horizontalalignment='center',verticalalignment='center', transform=plt.gca().transAxes)
                elif row != col: # Plot the two dimensions against each other
                    plt.scatter(data[:,col],data[:,row],**kwargs)
                if variable_labels is not None: # Variable labels were specified
                    if row == D-1:
                        plt.xlabel(variable_labels[col])
                    if col == 0:
                        plt.ylabel(variable_labels[row])
                if axis_limits is not None:
                    if axis_limits[row] is not None:
                        plt.ylim(axis_limits[row])
                    if axis_limits[col] is not None:
                        plt.xlim(axis_limits[col])
    else: # If data array is three-dimensional (a time series)
        for row in range(D):
            for col in range(D):
                plt.subplot(D,D,row*D+col+1)
                if row == col and variable_names is not None: # If this point lies along the diagonal
                    plt.text(0.5, 0.5, variable_names[row], horizontalalignment='center',verticalalignment='center', transform=plt.gca().transAxes)
                elif row != col: # Plot the two time series against each other
                    for n in range(N): # One entry at a time
                        plt.plot(data[n,col,:],data[n,row,:],**kwargs)
                if variable_labels is not None: # Variable labels were specified
                    if row == D-1:
                        plt.xlabel(variable_labels[col])
                    if col == 0:
                        plt.ylabel(variable_labels[row])
                if axis_limits is not None:
                    if axis_limits[row] is not None:
                        plt.ylim(axis_limits[row])
                    if axis_limits[col] is not None:
                        plt.xlim(axis_limits[col])
    return

#%%   
def multivariate_grid(data,parameters,variable_names = None,variable_labels=None,axis_limits=None,**kwargs):
    
    """
    This function plots multivariate data as a collection of 2-D plots.
    
    Parameters:
        
        data            - Required  : a 2-D array of dimenions (number datapoints x dimensions parameter space)
                                      or a 3-D array of dimensions (number datapoints x dimensions parameter space x time)
        parameters      - Required  : a list of (dimensions parameter space) entries, each entry being a (dimensions parameter space)-dimensional rectange of (levels of discretization)
        variable_names  - Optional  : list of strings specifying the names of the variables
        variable_labels - Optional  : list of strings specifying the axis labels for the variables
        axis_limits     - Optional  : list of 2-element-lists (for lower and upper specified axis limits) or None 
                                      (for unspecified axis limits)
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Find data dimensions 
    D = len(data.shape)
        
    # If variable names or labels are specified, check if dimensions match
    text = ''
    if variable_names is not None:
        if len(variable_names) != D:
            text += '\nNumber of variable name entries does not match data dimensions!'
    if variable_labels is not None:
        if len(variable_labels) != D:
            text += '\nNumber of variable axis label entries does not match data dimensions!'
    if axis_limits is not None:
        if len(axis_limits) != D:
            text += '\nNumber of variable axis limit entries does not match data dimensions!'
    if len(text) > 0:
        raise Exception(text)
    
    all_D = np.arange(D,dtype=int)
    
    # Now plot everything
    for row in range(D):
        for col in range(D):
            plt.subplot(D,D,row*D+col+1)
            if row == col and variable_names is not None: # If this point lies along the diagonal
                plt.text(0.5, 0.5, variable_names[row], horizontalalignment='center',verticalalignment='center', transform=plt.gca().transAxes)
            elif row != col: # Plot the two dimensions against each other
                
                X = np.mean(parameters[col],axis=tuple([x for x in all_D if (x != row) and (x != col)]))
                Y = np.mean(parameters[row],axis=tuple([x for x in all_D if (x != row) and (x != col)]))
                Z = np.mean(data,axis=tuple([x for x in all_D if (x != row) and (x != col)]))
                
                plt.contourf(X,Y,Z)
                
#                plt.imshow(Z)
                
            if variable_labels is not None: # Variable labels were specified
                if row == D-1:
                    plt.xlabel(variable_labels[col])
                if col == 0:
                    plt.ylabel(variable_labels[row])
#            if axis_limits is not None:
#                if axis_limits[row] is not None:
#                    plt.ylim(axis_limits[row])
#                if axis_limits[col] is not None:
#                    plt.xlim(axis_limits[col])
    return