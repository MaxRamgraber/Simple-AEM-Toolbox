#import numpy as np
#from mpmath import mpc,mpmathify,ellipfun,acos,ellipf
#import matplotlib.pyplot as plt


#plt.close('all')

class Model:
    
    def __init__(self,head_offset=0,aquifer_type='unconfined',domain_center=0+0j,
                 domain_radius=1,H = None,variables=[],priors=[],observations=[]):
        
        """
        This creates a model base object, to which we can add other elements.
        
        Parameters:
            
        head_offset     - [scalar]  : aquifer base elevation in [length units]
        aquifer_type    - [string]  : specifies the aquifer type; either 'confined', 'unconfined', or 'convertible'
        domain_center   - [complex] : x + iy coordinate of center of the circular, physical domain in [length units]; can also be specified as a vector of length 2
        domain_radius   - [scalar]  : radius of the circular domain in [length units]
        H               - [scalar]  : aquifer top elevation in [length units]; only used if the aquifer is 'confined' or 'convertible'
        
        If MCMC use is intended, we further require:
        
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['head_offset','H']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        observations    - [list]    : list of dictionaries, one for each hydraulic head observations; each dictionary must contain a 'location' and a 'head key', with a complex and real number, respectively
        """
        
        import numpy as np
        
        # Set potential scaling variables
        self.head_offset    = head_offset
        self.aquifer_type   = aquifer_type
        self.H              = H
        
        # Set domain scaling variables
        self.domain_center  = domain_center
        self.domain_radius  = domain_radius
        
        if not np.isscalar(self.domain_center):
            self.domain_center  = self.domain_center[0] + 1j*self.domain_center[1]
        
        # Check input for validity
        self.check_input()
        
        # Define a list for Analytic Elements
        self.elementlist    = []
        
        self.variables      = variables
        self.priors         = priors
        self.observations   = observations
        
        # This function scrapes the model and its elements for unknown variables,
        # then gives this instance three new variables:
        #   self.num_params         Number of unknown variables
        #   self.params             List of unknown variables
        #   self.param_names        List of names of unknown variables
        #   self.priors             List of prior dictionaries for unknow variables
        self.take_parameter_inventory()
        
        self.linear_solver  = False
        
        # Pre-allocate the function matrix and parameter vector for the linear solver
        self.matrix_solver  = []
        self.params_vector  = []
        
    def update(self):
        
        import copy
        
        # self.take_parameter_inventory()
        
        # Unpack the parameter vector to their respective instances -----------
        
        # Count through the parameters
        counter     = -1
        
        # Go through all unknown variables in the model class, if any
        for var in self.variables:
            
            # Replace the old variable
            counter += 1
            exec("self.%s = copy.copy(self.params[counter])" % var)
           
        # Count through the parameters
        counter     = -1
           
        # Go through all elements
        for e in self.elementlist:
            
            # Go through all the element's unknown variables
            for var in e.variables:
                
                # Replace the old variable
                counter += 1
                exec("e.%s = copy.copy(self.params[counter])" % var)
              
            # Update all other elements
            e.update()   
            
               
    def check_input(self):
            
        # Check if aquifer type is valid
        if self.aquifer_type != 'confined' and \
           self.aquifer_type != 'unconfined' and \
           self.aquifer_type != 'convertible':
            raise Exception("aquifer_type must be either 'confined', 'unconfined', or 'convertible'.")
            
        if (self.aquifer_type == 'confined' or self.aquifer_type == 'convertible') and \
            self.H is None:
            raise Exception("depth of confined layer 'H' must be specified if aquifer is confined or convertible.")
            
    def evaluate(self,z,mode='potential',derivatives='all',return_error_flag=False, suppress_warnings = False):
        
        import numpy as np
        import copy
        
        # Ensure that the evaluation points are complex
        z   = self.complexify(z)
        
        if return_error_flag:
            error_flag  = False
        
        self.update()
        
        # Inverse maps from disk to square,
        # Not inverse maps from square to disk
        
        # If there is at least one prescribed head element, prepare the linear solver
        if self.linear_solver:
            
            matrix,solution_vector = self.set_up_linear_system()
            
            # Find all elements which require the solver
            part_of_solver      = [(isinstance(e, ElementHeadBoundary) or isinstance(e, ElementNoFlowBoundary) or isinstance(e, ElementInhomogeneity)) for e in self.elementlist]
            part_of_solver      = [idx for idx,val in enumerate(part_of_solver) if val]
            not_part_of_solver  = [i for i in np.arange(len(self.elementlist)) if i not in part_of_solver]
            
            # Solve the system of linear equations
            param_vec   = np.linalg.solve(matrix,solution_vector)
            
            # Assign those parameters to each element
            counter = 0
            for idx in part_of_solver:
                
                # Extract the current element...
                e   = self.elementlist[idx]
                
                # ...and assign the correct strength
                e.strength = copy.copy(param_vec[counter:counter+e.segments])
            
                # Then update the entry counter
                counter += e.segments
        
        # =====================================================================
        # Now that the coefficients are set, evaluate the results
        # =====================================================================
        
        if mode == 'potential':
        
            # Coordinates in canonical space are the start values
            z_canonical     = copy.copy(z)
            
            z = np.zeros(z.shape,dtype=np.complex)
      
            # z *= 0 #-marked-
            for e in self.elementlist:
                z += e.evaluate(z_canonical)
                
        elif mode == 'gradient':
            
            # Coordinates in canonical space are the start values
            z_canonical     = copy.copy(z)
            
            z = np.zeros(z.shape,dtype=np.complex)
      
            # z *= 0 #-marked-
            for e in self.elementlist:
                z += e.evaluate_gradient(z_canonical,derivatives=derivatives)

        elif mode == 'head':
            
            # Coordinates in canonical space are the start values
            z_canonical     = copy.copy(z)
            
            z = np.zeros(z.shape,dtype=np.complex)
      
            # z *= 0 #-marked-
            for e in self.elementlist:
                z += e.evaluate(z_canonical)
             
            # First, get the base conductivity
            for e in self.elementlist:
                if isinstance(e, ElementMoebiusBase) or isinstance(e, ElementUniformBase):
                    temp_k = np.ones(z_canonical.shape)*e.k
            for e in self.elementlist:
                if isinstance(e, ElementInhomogeneity):
                    inside  = e.are_points_inside_polygon(z_canonical)
                    temp_k[inside] = e.k
                    
            # The hydraulic potential can never be negative; set it to zero 
            # (drying of an area) for any regions where it is negative, then
            # issue a warning
            if len(np.where(np.real(z) <= 0)[0]) > 0:
                if not suppress_warnings:
                    print('WARNING: negative or zero potential detected at some evaluation points. Consider lowering head_offset or prescribing prior limits.')
                z[np.where(np.real(z) <= 0)] = 1j*np.imag(z[np.where(np.real(z) <= 0)])
                if return_error_flag:
                    error_flag  = True
                
            if self.aquifer_type == 'confined':
                
                # Strack 1989, Eq. 8.12
                z   = (np.real(z) + 0.5*temp_k*self.H**2)/(temp_k*self.H) + \
                    1j*np.imag(z)
                    
            elif self.aquifer_type == 'unconfined':
                
                # Strack 1989, Eq. 8.13
                z   = np.sqrt(2*(np.real(z))/temp_k) + 1j*np.imag(z)
                
            elif self.aquifer_type == 'convertible':
                
                # Decide which equation to use for what points
                # confined:     Strack 1989, Eq. 8.12
                # unconfined:   Strack 1989, Eq. 8.13
                limit   = 0.5*temp_k/self.H**2
                index_conf      = np.where(np.real(z) >= limit)[0]
                index_unconf    = np.where(np.real(z) < limit)[0]
                
                # Handle the confined part
                z[index_conf] = \
                    (np.real(z[index_conf]) + 0.5*temp_k[index_conf]*self.H**2)/(temp_k[index_conf]*self.H) + \
                    1j*np.imag(z[index_conf])
                    
                # Handle the unconfined part
                z[index_unconf] = \
                    np.sqrt(2*(np.real(z[index_unconf]))/temp_k[index_unconf]) + 1j*np.imag(z[index_unconf])
                    
            # Offset the head
            z   += self.head_offset

        else:
            
            raise Exception("Mode must be either 'potential', 'gradient', or 'head'.")
            
        if return_error_flag:
            return z,error_flag
        else:
            return z
    
    def set_up_linear_system(self):
        
        """
        This function sets up the system of linear equations required to solve
        for the unknown coefficients of prescribed head boundaries, no-flow
        boundaries, and polygonal inhomogeneities.
        """
        
        import numpy as np
        import copy
        
        # Find all elements which require the solver
        # First, find all elements which are either Line Sinks, Doublets, or Inhomogeneities
        part_of_solver      = [(isinstance(e, ElementHeadBoundary) or isinstance(e, ElementNoFlowBoundary) or isinstance(e, ElementInhomogeneity)) for e in self.elementlist]
        # Only keep the elements which must be part of the linear system...
        part_of_solver      = [idx for idx,val in enumerate(part_of_solver) if val]
        # ...and prepare a second set of indices for its complement
        not_part_of_solver  = [i for i in np.arange(len(self.elementlist)) if i not in part_of_solver]
        
        # These elements invariably consist of segments - find out how many there are in total
        num_segments        = np.sum([self.elementlist[idx].segments for idx in part_of_solver])
        
        # =====================================================================
        # Now create the matrix
        # =====================================================================
        
        # Pre-allocate arrays for the linear solver
        matrix              = np.zeros((num_segments,num_segments))
        
        # The counter will keep track at what row we are
        row = 0
        
        # Go through all elements
        for i in part_of_solver:
            
            # Find the corresponding element
            e   = self.elementlist[i]
            
            # We need a second counter for the columns
            col = 0
            
            # e is the element we are currently looking at - the row -, now we 
            # must go through all other elements which are part of the solver
            # and check what they contribute to the control points of this element
            for i2 in part_of_solver:
                
                # Find the corresponding element
                e2  = self.elementlist[i2]
            
                # If the row element is a HeadLineSink, we must extract potentials
                if isinstance(e, ElementHeadBoundary):
                    
                    # Evaluate the contributions of this element to the control points
                    if e != e2:
                        block   = e2.evaluate(
                            z                   = e.zc,
                            detailed            = True,
                            override_parameters = True).T
                    else:
                        block   = e2.evaluate(
                            z                   = e.zc,
                            detailed            = True,
                            override_parameters = True,
                            evaluate_self       = True).T
                        
                
                elif isinstance(e, ElementNoFlowBoundary):
                    
                    # Evaluate the contributions of this element to the control points
                    block   = e2.evaluate_gradient(
                        z                   = e.zc,
                        detailed            = True,
                        derivatives         = 'phi',
                        override_parameters = True).T
                    
                    # Project the partial derivatives onto the normal vector
                    # The projection is a->b = <a,b>/||b||^2*b
                    # Let's try it with the inner product instead
                    # The normal vector is already normalized
    
                    # We should have as many normal vectors as we have control points
                    # Go through them all, and project each gradient onto the normal vector
                    for idx,nv in enumerate(e.segment_nvec):
                        
                        # Calculate the inner product between the returned partial
                        # derivatives and the segment's normal vector
                        block[idx,:] = np.inner(
                            np.column_stack(( 
                                np.real(block[idx,:]),
                                np.imag(block[idx,:]) )),
                            np.asarray([np.real(nv),np.imag(nv)]).T )[:,0]
                        
                elif isinstance(e, ElementInhomogeneity):
                    
                    # If this inhomogeneity evaluates itself
                    if i == i2:
            
                        # Retrieve own matrix contribution
                        block   = copy.copy(e2.block)
                        
                        # This contribution is incomplete, subtract A_star from
                        # its diagonal
                        
                        # Prepare a vector of outside conductivities; all are
                        # the background conductivity by default
                        for e3 in self.elementlist:
                            if isinstance(e3, ElementMoebiusBase) or isinstance(e3, ElementUniformBase):
                                A_star  = np.ones(e2.zc.shape)*e3.k/(e2.k - e3.k)
                        
                        # Get add matrix
                        addmat  = np.identity(block.shape[0])
                        np.fill_diagonal(addmat,A_star)
                        
                        # Subtract it from the retrieved block
                        block   -= addmat
                        
                    else:
            
                        # Evaluate the contributions of this element to the control points
                        block   = e2.evaluate(
                            z                   = e.zc,
                            detailed            = True,
                            override_parameters = True).T
                                    
                # Write this block into the matrix
                matrix[row:row+e.segments,col:col+e2.segments] = copy.copy(np.real(block))
            
                # Update the column counter
                col     += e2.segments
                
            # Update the row counter
            row     += e.segments
            
        # =====================================================================
        # Now create the solution_vector
        # =====================================================================
        
        # Pre-allocate spac efor the solution vector
        solution_vector     = np.zeros(num_segments)
        
        # The counter will keep track at what row we are
        counter = 0
        
        # Go through all elements
        for i in part_of_solver:
            
            # Find the corresponding element
            e   = self.elementlist[i]
            
            # If the element is a HeadLineSink, we must assign the difference
            # between the head target and the background contributions
            if isinstance(e, ElementHeadBoundary):
                
                # Step 1: Assign the head targets
                solution_vector[counter:counter+e.segments] = \
                    copy.copy(e.phi_target)
                # solution_vector[counter:counter+e.segments] = \
                #     copy.copy(e.head_target)
                
                # # Step 2: Background potential --------------------------------
                # solution_vector[counter:counter+e.segments] -= \
                #     np.real(self.evaluate(e.zc))
                    
                # Step 3: All elements ----------------------------------------
                for idx in not_part_of_solver:
                    solution_vector[counter:counter+e.segments] -= \
                        np.real(self.elementlist[idx].evaluate(e.zc))
                        
            # If the element is a no-flow boundary, we must assign the difference
            # between the head target and the background contributions
            if isinstance(e, ElementNoFlowBoundary):
                
                # # Step 1: Background gradient ---------------------------------
                # temp = self.evaluate_gradient(e.zc,derivatives='phi')
                
                # Step 2: Gradients from all elements -------------------------
                temp    = np.zeros(e.zc.shape,dtype=np.complex)
                for idx in not_part_of_solver:
                    temp += \
                        self.elementlist[idx].evaluate_gradient(e.zc,derivatives='phi')
                    
                # Step 3: Project gradients onto normal vector ----------------
                for ix,nv in enumerate(e.segment_nvec):
                    solution_vector[counter+ix] = \
                        -np.inner(
                            np.asarray([np.real(nv),np.imag(nv)])[:,0],
                            np.asarray([np.real(temp[ix]),np.imag(temp[ix])]) )
        
            # If the element is an Inhomogeneity, we must simply assign the potentials
            # induced by other elements
            if isinstance(e, ElementInhomogeneity):
                
                # # Step 1: Background potential --------------------------------
                # solution_vector[counter:counter+e.segments] -= \
                #     np.real(self.evaluate(e.zc))
                    
                # Step 2: All elements ----------------------------------------
                for idx in not_part_of_solver:
                    solution_vector[counter:counter+e.segments] -= \
                        np.real(self.elementlist[idx].evaluate(e.zc))
                        
            # Update the counter
            counter += e.segments
            
        self.matrix = matrix
        self.solvec = solution_vector
            
        return matrix, solution_vector
        
    
    def gradients(self,z):
        
        import numpy as np
        
        # Extract the gradients and return them
        grad    = np.zeros(z.shape,dtype=np.complex)
        
        for e in self.elementlist:
            
            grad    += e.evaluate_gradient(z)
        
        return grad
    
    def take_parameter_inventory(self):
        
        # Find the number of unknown variables
        self.num_params     = 0
        self.params     = []
        self.param_names    = []
        
        # First see if they are any unknown variables in the main model
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.num_params += 1
                exec("self.params += [self.%s]" % var)
                if 'name' in list(self.priors[idx].keys()):
                    self.param_names    += [self.priors[idx]['name']]   
                else:
                    self.param_names    += ['unknown']  
                
        # Check if the prior matches the number of parameters
        if len(self.priors) != self.num_params:
            raise Exception('Number of priors must match number of parameters. Number of parameters:'+str(self.num_params)+' / Number of priors:'+str(len(self.priors)))

    def logprior(self,params,priors,verbose=False):
        
        import numpy as np
        import scipy.stats
        import copy
        import math
        from toolbox_AEM import ElementMoebiusBase,ElementMoebiusOverlay
        
        def check_limits(params,var_dict):
            
            reject = None
            import numpy as np
            
            # Check if any limits are prescribed
            if 'limits' in list(var_dict.keys()): 
                
                if var_dict['limits'][0] is not None and type(var_dict['limits'][0]) != str:
                    
                    if np.isscalar(params):
                        if params <= var_dict['limits'][0]:
                            reject  = True
                            
                    else:
                        for entry in params:
                            if entry <= var_dict['limits'][0]:
                                reject  = True
                        
                if var_dict['limits'][1] is not None and type(var_dict['limits'][1]) != str:
                    
                    if np.isscalar(params):
                        if params >= var_dict['limits'][1]:
                            reject  = True
                            
                    else:
                        for entry in params:
                            if entry >= var_dict['limits'][1]:
                                reject  = True
                
                var_dict.pop('limits')
                
            return reject,var_dict
            
        # Find the base element
        MoebiusBase_index           = None
        for idx,e in enumerate(self.elementlist):
            if isinstance(e, ElementMoebiusBase) and 'r' in e.variables:
                MoebiusBase_index   = idx
                
        # Find any Möbius overlay element
        MoebiusOverlay_index        = None
        for idx,e in enumerate(self.elementlist):
            if isinstance(e, ElementMoebiusOverlay) and 'r' in e.variables:
                if MoebiusOverlay_index is None:
                    MoebiusOverlay_index    = [idx]
                else:
                    MoebiusOverlay_index    += [idx]
            
        
        logprior    = []
        reject      = False
        
        if MoebiusBase_index is not None:
            if self.elementlist[MoebiusBase_index].are_points_clockwise():
                # print('base clockwise')
                reject  = True
                
            # Check if the control points fulfill the minimum angular spacing
            r               = np.degrees(self.elementlist[MoebiusBase_index].r)
            angular_limit   = np.degrees(self.elementlist[MoebiusBase_index].angular_limit)
            if np.abs((r[0]-r[1] + 180) % 360 - 180) < angular_limit or \
               np.abs((r[1]-r[2] + 180) % 360 - 180) < angular_limit or \
               np.abs((r[2]-r[0] + 180) % 360 - 180) < angular_limit:
                # print('base angular limit violation')
                reject  = True
                
        if MoebiusOverlay_index is not None:
            for idx in MoebiusOverlay_index:
                if self.elementlist[idx].are_points_clockwise():
                    reject  = True
                    
                # Check if the control points fulfill the minimum angular spacing
                r               = np.degrees(self.elementlist[idx].r)
                angular_limit   = np.degrees(self.elementlist[idx].angular_limit)
                if np.abs((r[0]-r[1] + 180) % 360 - 180) < angular_limit or \
                   np.abs((r[1]-r[2] + 180) % 360 - 180) < angular_limit or \
                   np.abs((r[2]-r[0] + 180) % 360 - 180) < angular_limit:
                    reject  = True
                    
        # if not reject:
        #     print('WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOP')
                    

        for var in range(len(priors)):
            
            # If the logprior is to be rejected due to a violation of limits,
            # break the loop
            if reject:
                break   
        
            # Create a copy of this variable's prior dictionary
            var_dict  = copy.copy(priors[var])
            
            # Check if the user specified any converter
            converter_used = False
            if 'converter' in list(var_dict.keys()) and 'deconverter' in list(var_dict.keys()):
                
                # Activate the logarithmic boolean and save the base
                converter_used  = True
                
                # Extract the converter and the deconverter
                converter       = var_dict['converter']
                deconverter     = var_dict['deconverter']
                
                # Remove the keys from the dictionary
                var_dict.pop('converter')
                var_dict.pop('deconverter')
                
                # And convert the variable
                params[var] = converter(params[var])
                
            elif 'converter' in list(var_dict.keys()) or 'deconverter' in list(var_dict.keys()):
                raise Exception('Both a converter and a deconverter must be specified if either are used for a variable.')
                
            # Remove size if it was specified
            if 'size' in list(var_dict.keys()):
                del var_dict['size']
            
            # This prior is a univariate normal distribution
            if var_dict['distribution'] == 'norm' or var_dict['distribution'] == 'normal':
                
                # Remove the variable name from the dictionary
                var_dict.pop('distribution')
                if 'name' in list(var_dict.keys()): var_dict.pop('name')
                
                temp, var_dict = check_limits(params = params[var], var_dict = var_dict)
                if temp is not None:
                    reject  = True
                
                # Check if any limits are prescribed
                # if 'limits' in list(var_dict.keys()): 
                    
                #     if var_dict['limits'][0] is not None and type(var_dict['limits'][0]) != str:
                #         if params[var] < var_dict['limits'][0]:
                #             reject  = True
                            
                #     if var_dict['limits'][1] is not None and type(var_dict['limits'][1]) != str:
                #         if params[var] > var_dict['limits'][1]:
                #             reject  = True
                    
                #     var_dict.pop('limits')
        
                # Add to the logprior
                logprior += [np.sum(scipy.stats.norm.logpdf(x=params[var],**var_dict))]
                
            # This prior is a multivariate normal distribution
            elif var_dict['distribution'] == 'multivariate_normal' or var_dict['distribution'] == 'multivariate normal':
                
                # Remove the variable name from the dictionary
                var_dict.pop('distribution')
                if 'name' in list(var_dict.keys()): var_dict.pop('name')
                
                temp, var_dict = check_limits(params = params[var], var_dict = var_dict)
                if temp is not None:
                    reject  = True
                
                # # Check if any limits are prescribed
                # if 'limits' in list(var_dict.keys()): 
                    
                #     if var_dict['limits'][0] is not None and type(var_dict['limits'][0]) != str:
                #         if any(params[var] < var_dict['limits'][0]):
                #             reject  = True
                            
                #     if var_dict['limits'][1] is not None and type(var_dict['limits'][1]) != str:
                #         if any(params[var] > var_dict['limits'][1]):
                #             reject  = True
                    
                #     var_dict.pop('limits')
        
                # Add to the logprior
                logprior += [np.sum(scipy.stats.multivariate_normal.logpdf(x=params[var],**var_dict))]
        
            # This prior is a beta distribution
            elif var_dict['distribution'] == 'beta':
                
                # Remove the variable name from the dictionary
                var_dict.pop('distribution')
                if 'name' in list(var_dict.keys()): var_dict.pop('name')
                
                # A beta distribution has natural limits; if none are prescribed, add them
                if 'limits' not in list(var_dict.keys()): 
                    var_dict['limits'] = [0,1]
                
                temp, var_dict = check_limits(params = params[var], var_dict = var_dict)
                if temp is not None:
                    reject  = True
                
                # # Check if any limits are prescribed
                # if 'limits' in list(var_dict.keys()): 
                    
                #     if var_dict['limits'][0] is not None and type(var_dict['limits'][0]) != str:
                #         if params[var] < var_dict['limits'][0]:
                #             reject  = True
                            
                #     if var_dict['limits'][1] is not None and type(var_dict['limits'][1]) != str:
                #         if params[var] > var_dict['limits'][1]:
                #             reject  = True
                    
                #     var_dict.pop('limits')
        
                # Add to the logprior
                logprior += [np.sum(scipy.stats.beta.logpdf(x=params[var],**var_dict))]
        
            # This prior is an exponential distribution
            elif var_dict['distribution'] == 'expon' or var_dict['distribution'] == 'exponential':
                
                # SPECIAL EXCEPTION -------------------------------------------
                # The argument of the exponential distribution is considered as
                # its absolute value to permit evaluation of negative values.
                # Checking for limits happens before this, so negative limits
                # can be applied.
                # -------------------------------------------------------------
                
                # Remove the variable name from the dictionary
                var_dict.pop('distribution')
                if 'name' in list(var_dict.keys()): var_dict.pop('name')
                
                # An exponential distribution has natural limits; if none are prescribed, add them
                if 'limits' not in list(var_dict.keys()): 
                    var_dict['limits'] = [0,None]
                
                temp, var_dict = check_limits(params = params[var], var_dict = var_dict)
                if temp is not None:
                    reject  = True
                
                # # Check if any limits are prescribed
                # if 'limits' in list(var_dict.keys()): 
                    
                #     if var_dict['limits'][0] is not None and type(var_dict['limits'][0]) != str:
                #         if params[var] < var_dict['limits'][0]:
                #             reject  = True
                            
                #     if var_dict['limits'][1] is not None and type(var_dict['limits'][1]) != str:
                #         if params[var] > var_dict['limits'][1]:
                #             reject  = True
                    
                #     var_dict.pop('limits')
        
                # Add to the logprior
                logprior += [np.sum(scipy.stats.expon.logpdf(x=np.abs(params[var]),**var_dict))]
        
            # This prior is a von Mises distribution
            elif var_dict['distribution'] == 'vonmises' or var_dict['distribution'] == 'von mises' or var_dict['distribution'] == 'von Mises':
                
                # Remove the variable name from the dictionary
                var_dict.pop('distribution')
                if 'name' in list(var_dict.keys()): var_dict.pop('name')
                
                temp, var_dict = check_limits(params = params[var], var_dict = var_dict)
                if temp is not None:
                    reject  = True
                
                # # Check if any limits are prescribed
                # if 'limits' in list(var_dict.keys()): 
                    
                #     if var_dict['limits'][0] is not None and type(var_dict['limits'][0]) != str:
                #         if params[var] < var_dict['limits'][0]:
                #             reject  = True
                            
                #     if var_dict['limits'][1] is not None and type(var_dict['limits'][1]) != str:
                #         if params[var] > var_dict['limits'][1]:
                #             reject  = True
                    
                #     var_dict.pop('limits')
        
                # Add to the logprior
                logprior += [np.sum(scipy.stats.vonmises.logpdf(x=params[var],**var_dict))]
                
            else:
                
                raise Exception("Specified distribution name '" + str(var_dict['distribution']) + \
                                "' not understood. Valid distribution names are: 'norm', " + \
                                "'multivariate_normal', 'beta', 'expon', or 'vonmises'")
                    
            # If a converter are used, deconvert the variables
            if converter_used:
                params[var] = deconverter(params[var])
                
        # Check if variables have been prescribed as limits
        counter     = -1
        for idx,e in enumerate(self.elementlist):
            
            # If the logprior is to be rejected due to a violation of limits,
            # break the loop
            if reject:
                break
            
            # At what value was the counter at the start of this element
            counter_elementstart = counter
        
            # Check the priors list
            for prior in e.priors:
                
                # Increment the counter variable
                counter     += 1
                
                # Check if this prior entry has prescribed limits
                if 'limits' in prior:
                    
                    # Check if the lower limit is a string (a variable)
                    if type(prior['limits'][0]) == str:
                            
                        # Check where this variable is, store it
                        limit   = None
                        for idx,var in enumerate(e.variables):
                            if prior['limits'][0] == var:
                                limit   = params[counter_elementstart+1+idx]
                                
                        if limit is None: raise Exception("variable '"+prior['limits'][0]+"' marked as limit not part of the variables of element "+str(e))
                                
                        if params[counter] < limit:
                            reject      = True
                            
                    # Check if the upper limit is a string (a variable)
                    if type(prior['limits'][1]) == str:
                            
                        # Check where this variable is, store it
                        limit   = None
                        for idx,var in enumerate(e.variables):
                            if prior['limits'][1] == var:
                                limit   = params[counter_elementstart+1+idx]
                                
                        if limit is None: raise Exception("variable '"+prior['limits'][1]+"' marked as limit not part of the variables of element "+str(e))
                                
                        if params[counter] > limit:
                            reject      = True
                        
    
        # Return the logprior only if the sample isn't rejected
        if not reject:
            res     = logprior
        else:
            res     = None
            if verbose:
                print('Logprior calculation rejected because at least one variable violated prescribed limits.')
        
        return res
    
    def loglikelihood(self,observations,likelihood_dictionary,predictions = None):
        
        import numpy as np
        import scipy.stats
        import copy
        
        loglikelihood   = None
        
        obs_dict        = copy.deepcopy(likelihood_dictionary)
        
        # If no predictions have been provided, map forward
        if predictions is None:
            
            # Get the well positions
            z   = []
            for entry in observations:
                z += [copy.copy(entry['location'])]
            z   = np.asarray(z)
            
            predictions,error_flag = copy.copy(np.real(self.evaluate(
                z,
                mode='head',
                return_error_flag=True,
                suppress_warnings=True)))
            
            # If any of the predictions is NaN, raise an error flag
            if not error_flag and np.isnan(predictions).any():
                error_flag  = True
            
        else:
            error_flag = False
            
        predictions = np.asarray(predictions)
        
        # Create a vector of observations
        obs_vec     = []
        for entry in observations:
            obs_vec += [copy.copy(entry['head'])]
        obs_vec     = np.asarray(obs_vec)
            
        # Get the prediction residuals
        residuals   = obs_vec - predictions
        
        # This prior is a von Mises distribution
        if obs_dict['distribution'] == 'norm' or obs_dict['distribution'] == 'normal':
                
            # Remove superfluous keys
            obs_dict.pop('distribution')
            if 'name' in list(obs_dict.keys()):
                obs_dict.pop('name')
            if 'loc' in list(obs_dict.keys()):
                print("Warning: 'loc' specified for the loglikeligood pdf. This value is overwritten by the predictions.")
                obs_dict.pop('loc')
                
            # Add to the logprior
            loglikelihood = scipy.stats.norm.logpdf(x=obs_vec,loc=predictions,**obs_dict)
            
        # This prior is a von Mises distribution
        elif obs_dict['distribution'] == 'multivariate_normal' or obs_dict['distribution'] == 'multivariate normal':
            
            # Remove superfluous keys
            obs_dict.pop('distribution')
            if 'name' in list(obs_dict.keys()):
                obs_dict.pop('name')
            if 'mean' in list(obs_dict.keys()):
                print("Warning: 'mean' specified for the loglikeligood pdf. This value is overwritten by the predictions.")
                obs_dict.pop('mean')
                
            # Add to the logprior
            loglikelihood = np.sum(scipy.stats.multivariate_normal.logpdf(x=obs_vec,mean=predictions,**obs_dict))
        
        if error_flag:
            loglikelihood = None
            
        return loglikelihood, residuals
    
    def complexify(self,z):
        
        """
        This function takes the provided set of points and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,color='xkcd:grey'):
        
        """
        This function plots all elements in the elementlist
        """
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Plot the model domain
        plt.plot(
            np.cos(np.linspace(0,2*np.pi,361))*self.domain_radius + np.real(self.domain_center),
            np.sin(np.linspace(0,2*np.pi,361))*self.domain_radius + np.imag(self.domain_center),
            color=color)
        
        for e in self.elementlist:
            e.plot()
            
    def trace_gradient(self,p,direction='upgradient',stepsize=1,well_snap_distance = 1):
        
        """
        This function takes a regular grid defined by X and Y meshgrid arrays and 
        follows the gradient of a corresponding Z array starting from a point p.   
        
            p               : starting point for gradient tracing
            XY              : array of points between which the gradient tracing is interpolated
            Z               : hydraulic head at XY
            direction       : whether we trace 'upgradient' or 'downgradient'
            stepsize        : size of successive steps for gradient tracing
            visualize       : plots the results, if desired
            Del             : Delaunay triangulation, calculated if missing
            triang          : triangulation, calculated if missing
            thresh          : side length threshold above which vertices are rejected
        
        """
        
        if not direction == 'upgradient' and not direction == 'downgradient':
            raise Exception("direction must be either 'upgradient' or 'downgradient'.")
            
        import scipy.spatial
        import numpy as np
        import matplotlib.pyplot as plt
        import shapely.geometry
        
        ring    = np.column_stack((
            np.cos(np.linspace(0,2*np.pi,361)),
            np.sin(np.linspace(0,2*np.pi,361)) )) 
        ring    *= self.domain_radius
        ring    += np.asarray([np.real(self.domain_center),np.imag(self.domain_center)])
        
        # First, find all elements which could be stoppers
        stoppers    = []
        stoppers.append(shapely.geometry.LineString(ring))
        for e in self.elementlist:
            
            if isinstance(e, ElementHeadBoundary):
                # Head Boundaries are valid end points
                stoppers.append(shapely.geometry.LineString(e.line[:,:2]))
                
            if isinstance(e, ElementWell):
                # Wells are valid end points
                stoppers.append(shapely.geometry.Point(np.asarray([np.real(e.zc),np.imag(e.zc)])))
                
            if isinstance(e, ElementLineSink):
                # Line Sinks are valid end points
                stoppers.append(shapely.geometry.LineString(e.line[:,:2]))
                
            if isinstance(e, ElementNoFlowBoundary):
                # No-flow Boundaries are valid end points
                stoppers.append(shapely.geometry.LineString(e.line[:,:2]))
        
        def gradient(p1,p2,p3,z1,z2,z3):
            
            area    = abs((p1[0]*(p2[1]-p3[1])+p2[0]*(p3[1]-p1[1])+p3[0]*(p1[1]-p2[1]))/2)
            
            M      = np.asarray(
                [[p2[1]-p3[1],      p3[1]-p1[1],    p1[1]-p2[1]],
                 [p3[0]-p2[0],      p1[0]-p3[0],    p2[0]-p1[0]]])
        
            U       = np.asarray([z1,z2,z3]).reshape((3,1))
            
            # Solution based on http://pers.ge.imati.cnr.it/livesu/papers/MLP18/MLP18.pdf Equation 1
            return np.dot(M,U)[:,0]/(2*area)
        
        # Check if the start point is complex, if yes, turn it into a real vector
        if np.iscomplex(p).any():
            p   = np.asarray([np.real(p),np.imag(p)])
            
        # Depending on the direction, add a gradient
        if direction == 'upgradient':
            stepsize    = stepsize
        else:
            stepsize    = -stepsize
        
        # Set the repeater boolean to True
        repeater = True
        
        
        
        # Re-arrange the starting point into an array
        points = np.asarray(p).copy().reshape((1,2))
        # """
        # Get three points 
        testpoints  = np.asarray([
            points[-1,0] + 1j*points[-1,1],
            points[-1,0] + stepsize/100 + 1j*points[-1,1],
            points[-1,0] + 1j*points[-1,1] + 1j*stepsize/100])
        
        testpoints  = np.real(self.evaluate(testpoints,mode='head'))
        
        grad        = np.asarray([
            testpoints[1]-testpoints[0],
            testpoints[2]-testpoints[0]])/stepsize*100
        grad    = grad/np.linalg.norm(grad)
        # """
        
        # grad    = self.evaluate(
        #     z           = points,
        #     mode        = 'gradient',
        #     derivatives = 'phi')
        # # grad    = np.asarray([np.real(grad), np.imag(grad)])
        # grad    = grad/np.linalg.norm(grad)
        
        # And save the result to the points array
        points = np.row_stack((
            points.copy(),
            points + grad*stepsize))
        
        # Now start the while loop, trace until the end
        while repeater:
                
            # The last point in the array is the starting point
            p   = points[-1,:]
            
            # """
            testpoints  = np.asarray([
                points[-1,0] + 1j*points[-1,1],
                points[-1,0] + stepsize/100 + 1j*points[-1,1],
                points[-1,0] + 1j*points[-1,1] + 1j*stepsize/100])
            
            testpoints  = np.real(self.evaluate(testpoints,mode='head'))
            
            grad        = np.asarray([
                testpoints[1]-testpoints[0],
                testpoints[2]-testpoints[0]])/stepsize*100
        
            grad    = grad/np.linalg.norm(grad)
            # """
            
            # grad    = self.evaluate(
            #     z           = points[-1,:],
            #     mode        = 'gradient',
            #     derivatives = 'phi')
            # # grad    = np.asarray([np.real(grad), np.imag(grad)])
            # grad    = grad/np.linalg.norm(grad)
            
            # And append the next step to the list
            points = np.row_stack((
                points,
                points[-1,:] + grad*stepsize))
                
            
            line    = shapely.geometry.LineString(points[-2:,:])
            
            # Check for stopping elements
            for stop in stoppers:
                
                # If this stopper is a well, check for distance
                if stop.type == 'Point':
                    point = shapely.geometry.Point(points[-1,:])
                    if point.distance(stop) <= well_snap_distance:
                        points[-1,:] = np.asarray(point.xy)[:,0]
                        repeater = False
                
                # Else, we can check for intersection
                else:
                    if line.intersects(stop):
                        
                        if line.intersection(stop).type == 'Point':
                        
                            points[-1,:] = np.asarray(line.intersection(stop).xy)[:,0]
                            repeater = False
                            
                        else:
                            
                            print(type(line.intersection(stop)))
                            print((type(line.intersection(stop)) == 'Point'))
                            
                            points[-1,:] = np.asarray(line.intersection(stop)[0].xy)[:,0]
                            repeater = False

#            # Check for oscillation
#            p2p     = points[-3,:]-points[-2,:]
#            p1p     = points[-2,:]-points[-1,:]
#            if np.inner(p1p,p2p) < 0: 
#                # The trace direction has change by more than 90 degrees, i.e.
#                # turned back; stop iterating
#                points = points[:-1,:]
#                repeater = False
        
        return points
        

#%%

class ElementMoebiusBase:
    
    def __init__(self,model,r=None,a=None,b=None,c=None,d=None,head_min=0,head_max=1,
        k=1,variables=[],priors=[],proposals=[],angular_limit=1,use_SC=True):
        
        
        """
        This implements the Möbius base flow element, which can induce curving,
        converging, or diverging regional flow.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        r               - [vector]  : rotations of the three Moebius control points in counter-clockwise radians from East
        a               - [scalar]  : coefficient of the Moebius transformation; calculated from r if not specified
        b               - [scalar]  : coefficient of the Moebius transformation; calculated from r if not specified
        c               - [scalar]  : coefficient of the Moebius transformation; calculated from r if not specified
        d               - [scalar]  : coefficient of the Moebius transformation; calculated from r if not specified
        head_min        - [scalar]  : minimum hydraulic head (mapped to -1 in the unit square)
        head_max        - [scalar]  : maximum hydraulic head (mapped to +1 in the unit square)
        k               - [scalar]  : background hydraulic conductivity in canonical units (e.g., 1E-4 [length]/[time])
        use_SC          - [boolean] : a flag to denote whether the Möbius base uses the Schwarz-Christoffel transformation from the unit square to the unit disk; changes to this flag affect the flow field, particularly near the edges


        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['r','phi_min']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        angular_limit   - [scalar]  : a limit which prevents control points (A,B,C, or D) getting closer to each other than the specified value in radians; this acts as protection against improbable or unrealistic flow fields induced by the Möbius transformation
        """
        
        import numpy as np
        
        # Append the base to the elementlist
        self.model          = model
        model.elementlist.append(self)
        
        # Define an angular limit. This is designed to keep the Möbius control 
        # points from getting arbitrarily close to each other; defined in radians
        self.angular_limit  = angular_limit
        
        # Get the Schwarz-Christoffel flag
        self.use_SC         = use_SC
        
        # Set Moebius values
        self.r              = r
        self.a              = a
        self.b              = b
        self.c              = c
        self.d              = d
        
        # Set potential scaling variables
        self.head_min       = head_min
        self.head_max       = head_max
        
        # Assign the hydraulic conductivity of the base model
        self.k              = k
        
        # The model requires the base flow in terms of hydraulic potential (phi)
        # The function head_to_potential extracts the following variables:
        #   phi_min         hydraulic potential corresponding to head_min
        #   phi_max         hydraulic potential corresponding to head_max
        self.head_to_potential()
        
        # Check input for validity
        self.check_input()
        
        # Define the original control points in the Moebius base disk
        self.z0     = np.asarray(
            [np.complex(np.cos(-0.25*np.pi),np.sin(-0.25*np.pi)),
             np.complex(np.cos(0.25*np.pi),np.sin(0.25*np.pi)),
             np.complex(np.cos(0.75*np.pi),np.sin(0.75*np.pi))])
        
        # If only rotation is specified, get the Moebius coefficients
        if self.r is not None and (self.a is None and self.b is None and \
                                   self.c is None and self.d is None):
            # Find Moebius coefficients
            self.find_moebius_coefficients()
        
        self.variables      = variables
        self.priors         = priors
        self.proposals      = proposals
        
        self.Ke             = 1.854
        
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.variables  += [var]
                self.model.priors     += [self.priors[idx]]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']  
        
    def update(self):
        
        # If this model is updated, make sure to repeat any initialization
        # Find Moebius coefficients
        self.find_moebius_coefficients()
        
        # Flip h_min and h_max, if necessary
        if self.head_min > self.head_max:
            temp            = self.head_min
            self.head_min   = self.head_max
            self.head_min   = temp
            
        # Update potential
        self.head_to_potential()
            
    def check_input(self):
        
        import numpy as np
        
        # See if either control point rotations or a full set of Moebius
        # coefficients are specified
        if self.r is None and (self.a is None or self.b is None or \
                               self.c is None or self.d is None):
            raise Exception('Either control point rotations r or Moebius coefficients a, b, c, and d must be specified.')
        
        # Check if phi_min is smaller than phi_offset, switch if necessary
        if self.head_min > self.head_max:
            raise Exception('Minimum potential phi_min is larger than maximum potential phi_max.')
    
        # Check if the control points fulfill the minimum angular spacing
        r   = np.degrees(self.r)
        if np.abs((r[0]-r[1] + 180) % 360 - 180) < self.angular_limit or \
           np.abs((r[1]-r[2] + 180) % 360 - 180) < self.angular_limit or \
           np.abs((r[2]-r[0] + 180) % 360 - 180) < self.angular_limit:
            raise Exception('Control points '+str(self.r)+' are too close to each other. Define different control points or adjust the angular limit: '+str(self.angular_limit))
    
    def evaluate(self,z):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Coordinates in canonical space are the start values
        z_canonical     = copy.copy(z)
        
        # Scale the canonical disk to unity canonical disk
        z = (z - self.model.domain_center)/self.model.domain_radius
        
        # Map from canonical disk to Möbius base
        z = self.moebius(z,inverse=True)
        
        # Map from Möbius base to unit square
        if self.use_SC:
            z = self.disk_to_square(z)
        
        # Rescale the complex potential
        z   = (z+1)/2 * (self.phi_max-self.phi_min) + self.phi_min
            
        return z
    
    def evaluate_gradient(self,z,derivatives = 'all'):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Map from the canonical disk to Möbius base
        z_mb = (copy.copy(z) - self.model.domain_center)/self.model.domain_radius
        
        # dz_mb / dz_c
        grad_4  = 1/self.model.domain_radius
        
        # Map from Möbius base to unit disk
        z_ud = self.moebius(copy.copy(z_mb),inverse=True)
        
        # dz_ud / dz_mb
        grad_3  = (self.a*self.d-self.b*self.c)/(self.c*z_mb-self.a)**2
        
        if self.use_SC:
            grad_2  = 2/(self.Ke*np.sqrt(z_ud**4+1))
        
        grad_1  = (self.phi_max-self.phi_min)/2

        if self.use_SC:                    
            grad    = grad_1*grad_2*grad_3*grad_4
        else:
            grad    = grad_1*grad_3*grad_4
            
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
            
        return grad
    
    def complex_integral(self,func,a,b):
        
        """
        This implements the Gauss-Kronrod integration for complex-valued functions.
        We use this to evaluate the Legendre incomplete elliptic integral of the 
        first kind, since it is about ten times as fast as using mpmath's ellipf
        function. Since this integration is a major computational bottleneck of
        this function, we stick with this approach.
        
        The equations below are adapted from: 
        https://stackoverflow.com/questions/5965583/use-scipy-integrate-quad-to-integrate-complex-numbers
        """
        
        import scipy
        from scipy import array
        
        def quad_routine(func, a, b, x_list, w_list):
            c_1 = (b-a)/2.0
            c_2 = (b+a)/2.0
            eval_points = map(lambda x: c_1*x+c_2, x_list)
            func_evals = list(map(func, eval_points))    # Python 3: make a list here
            return c_1 * sum(array(func_evals) * array(w_list))
        
        def quad_gauss_7(func, a, b):
            x_gauss = [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, 0.405845151377397, 0.741531185599394, 0.949107912342759]
            w_gauss = array([0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277,0.129484966168870])
            return quad_routine(func,a,b,x_gauss, w_gauss)
        
        def quad_kronrod_15(func, a, b):
            x_kr = [-0.991455371120813,-0.949107912342759, -0.864864423359769, -0.741531185599394, -0.586087235467691,-0.405845151377397, -0.207784955007898, 0.0, 0.207784955007898,0.405845151377397, 0.586087235467691, 0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813]
            w_kr = [0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728, 0.204432940075298, 0.190350578064785, 0.169004726639267, 0.140653259715525,  0.104790010322250, 0.063092092629979, 0.022935322010529]
            return quad_routine(func,a,b,x_kr, w_kr)
        
        class Memorize:                     # Python 3: no need to inherit from object
            def __init__(self, func):
                self.func = func
                self.eval_points = {}
            def __call__(self, *args):
                if args not in self.eval_points:
                    self.eval_points[args] = self.func(*args)
                return self.eval_points[args]
        
        def quad(func,a,b):
            ''' Output is the 15 point estimate; and the estimated error '''
            func = Memorize(func) #  Memorize function to skip repeated function calls.
            g7 = quad_gauss_7(func,a,b)
            k15 = quad_kronrod_15(func,a,b)
            # I don't have much faith in this error estimate taken from wikipedia
            # without incorporating how it should scale with changing limits
            return [k15, (200*scipy.absolute(g7-k15))**1.5]
        
        return quad(func,a,b)
    
    def angle_to_unit_circle(self):
        
        import numpy as np
        
        # Angle must be provided in radians, counter-clockwise from 3 o'clock
        return np.cos(self.r)+1j*np.sin(self.r)
    
    def find_moebius_coefficients(self):
        
        import numpy as np
        
        # Find the images of the z0 control points
        w0  = self.angle_to_unit_circle()
        
        # Then calculate the four parameters for the corresponding Möbius map
        self.a = np.linalg.det(np.asarray(
            [[self.z0[0]*w0[0],     w0[0],          1],
             [self.z0[1]*w0[1],     w0[1],          1],
             [self.z0[2]*w0[2],     w0[2],          1]]))
        
        self.b = np.linalg.det(np.asarray(
            [[self.z0[0]*w0[0],     self.z0[0],     w0[0]],
             [self.z0[1]*w0[1],     self.z0[1],     w0[1]],
             [self.z0[2]*w0[2],     self.z0[2],     w0[2]]]))
        
        self.c = np.linalg.det(np.asarray(
            [[self.z0[0],           w0[0],          1],
             [self.z0[1],           w0[1],          1],
             [self.z0[2],           w0[2],          1]]))
        
        self.d = np.linalg.det(np.asarray(
            [[self.z0[0]*w0[0],     self.z0[0],     1],
             [self.z0[1]*w0[1],     self.z0[1],     1],
             [self.z0[2]*w0[2],     self.z0[2],     1]]))
        
        return
    
    def moebius(self,z,inverse=False):
        
        if not inverse:
            z = (self.a*z+self.b)/(self.c*z+self.d)
        else:
            z = (-self.d*z+self.b)/(self.c*z-self.a)
        
        return z
    
    def square_to_disk(self,z,k='default'):
        
        import numpy as np
        from mpmath import mpc,mpmathify,ellipfun
        
        if k == 'default': k = 1/mpmathify(np.sqrt(2))
        
        Ke = 1.854
        cn = ellipfun('cn')
    
        if type(z) is complex: 
            z = np.asarray([z])
        zf  = np.ndarray.flatten(z)
        w   = np.zeros(zf.shape)*1j
        
        pre_factor  = mpc(1,-1)/mpmathify(np.sqrt(2))
        mid_factor  = Ke*(mpc(1,1)/2)
    
        for idx,entry in enumerate(zf): # Go through all complex numbers
            
            # Calculate result
            temp = pre_factor*cn(
                u = mid_factor*entry-Ke,
                k = k)
            
            # Then place it into the array
            w[idx]  = np.complex(temp.real,temp.imag)
            
        # Now reshape the array back to its original shape
        z = w.reshape(z.shape).copy()
        
        return z
    
    def disk_to_square(self,z,k='default'):
        
        import numpy as np
    
        Ke = 1.854

        if type(z) is complex: 
            z = np.asarray([z])
        zf  = np.ndarray.flatten(z)
        w   = np.zeros(zf.shape)*1j

        # Using the Gauss-Kronrod integration is about 10 times faster than 
        # using the mpmath.ellipf function
        if k == 'default': k = 1/np.sqrt(2)
        m = k**2
        pre_factor  = (1-1j)/(-Ke)
        mid_factor  = (1+1j)/np.sqrt(2)
        
        temp = [pre_factor*self.complex_integral(
            func    = lambda t: (1-m*np.sin(t)**2)**(-0.5), 
            a       = 0, 
            b       = i)[0] + 1 - 1j for i in np.arccos(zf*mid_factor)]
    
        w = np.asarray(temp)
            
        # Now reshape the array back to its original shape
        z = w.reshape(z.shape).copy()
        
        return z
    
    def are_points_clockwise(self):
        
        import numpy as np
        
        verts = np.zeros((3,2))
        
        verts[0,:] = np.asarray([np.cos(self.r[0]),np.sin(self.r[0])])
        verts[1,:] = np.asarray([np.cos(self.r[1]),np.sin(self.r[1])])
        verts[2,:] = np.asarray([np.cos(self.r[2]),np.sin(self.r[2])])
        
        signed_area = 0
        for vtx in range(verts.shape[0]):
            x1 = verts[vtx,0]
            y1 = verts[vtx,1]
            if vtx == verts.shape[0]-1: # Last vertex
                x2 = verts[0,0]
                y2 = verts[0,1]
            else:
                x2 = verts[vtx+1,0]
                y2 = verts[vtx+1,1]
            signed_area += (x1 * y2 - x2 * y1)/2
            
        return (signed_area < 0)
       
    def head_to_potential(self):
        
        for idx,h in enumerate([self.head_min-self.model.head_offset,self.head_max-self.model.head_offset]):
        
            if self.model.aquifer_type == 'confined' or (self.model.aquifer_type == 'convertible' and h >= self.model.H):
                # Strack 1989, Eq. 8.6
                pot = self.k*self.model.H*h - 0.5*self.k*self.model.H**2
                    
            elif self.model.aquifer_type == 'unconfined' or (self.model.aquifer_type == 'convertible' and h < self.model.H):
                # Strack 1989, Eq. 8.7
                pot = 0.5*self.k*h**2
                
            if idx == 0:
                self.phi_min    = pot
            elif idx == 1:
                self.phi_max    = pot
                
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,label_offset = 1.1,fontsize=12,fontcolor='xkcd:dark grey',
             pointcolor='xkcd:dark grey',pointsize=50,horizontalalignment='center',
             verticalalignment='center',color_low = 'xkcd:cerulean',
             color_high = 'xkcd:orangish red',**kwargs):
        
        """
        This function plots the Möbius control/reference points on the unit disk.
        """
        
        import numpy as np
        import matplotlib.pyplot as plt
        import math
        
        # Get the coordinates of the control points
        z_A     = (1-1j)/np.abs(1-1j)
        z_A     = self.moebius(z_A,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_A     = np.asarray([np.real(z_A),np.imag(z_A)])
        
        z_B     = (1+1j)/np.abs(1+1j)
        z_B     = self.moebius(z_B,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_B     = np.asarray([np.real(z_B),np.imag(z_B)])
        
        z_C     = (-1+1j)/np.abs(-1+1j)
        z_C     = self.moebius(z_C,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_C     = np.asarray([np.real(z_C),np.imag(z_C)])
        
        z_D     = (-1-1j)/np.abs(-1-1j)
        z_D     = self.moebius(z_D,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_D     = np.asarray([np.real(z_D),np.imag(z_D)])
        
        dc      = self.model.domain_center
        if np.isscalar(dc):
            dc  = np.asarray([np.real(dc),np.imag(dc)])
        
        a_low   = np.linspace(
            math.atan2(
                z_C[1]-dc[1],
                z_C[0]-dc[0]),
            math.atan2(
                z_D[1]-dc[1],
                z_D[0]-dc[0]),
            360)
        if abs(a_low[0]-a_low[-1]) > np.pi:
            a_low = np.concatenate((
                np.linspace(np.min(a_low),-np.pi,360),
                np.linspace(np.pi,np.max(a_low),360) ))
        
        a_high  = np.linspace(
            math.atan2(
                z_A[1]-dc[1],
                z_A[0]-dc[0]),
            math.atan2(
                z_B[1]-dc[1],
                z_B[0]-dc[0]),
            360)
        if abs(a_high[0]-a_high[-1]) > np.pi:
            a_high = np.concatenate((
                np.linspace(np.min(a_high),-np.pi,360),
                np.linspace(np.pi,np.max(a_high),360) ))
            
        plt.plot(np.cos(a_low)*self.model.domain_radius + dc[0],
                 np.sin(a_low)*self.model.domain_radius + dc[1],
                 color = color_low,linewidth=2)
        plt.plot(np.cos(a_high)*self.model.domain_radius + dc[0],
                 np.sin(a_high)*self.model.domain_radius + dc[1],
                 color = color_high,linewidth=2)
        
        plt.scatter(z_A[0],z_A[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        plt.scatter(z_B[0],z_B[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        plt.scatter(z_C[0],z_C[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        plt.scatter(z_D[0],z_D[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        
        plt.text(z_A[0]*label_offset,z_A[1]*label_offset,'A',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        plt.text(z_B[0]*label_offset,z_B[1]*label_offset,'B',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        plt.text(z_C[0]*label_offset,z_C[1]*label_offset,'C',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        plt.text(z_D[0]*label_offset,z_D[1]*label_offset,'D',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        
#%%

class ElementUniformBase:
    
    def __init__(self,model,alpha=0,head_min=0,head_max=1,k=1,
                 variables=[],priors=[]):
        
        
        """
        This implements the uniform base flow element.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        alpha           - [scalar]  : direction of the uniform flow in counter-clockwise radians from East
        head_min        - [scalar]  : minimum hydraulic head (mapped to -1 in the unit square)
        head_max        - [scalar]  : maximum hydraulic head (mapped to +1 in the unit square)
        k               - [scalar]  : background hydraulic conductivity in canonical units (e.g., 1E-4 [length]/[time])

        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['r','phi_min']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        observations    - [list]    : list of dictionaries, one for each hydraulic head observations; each dictionary must contain a 'location' and a 'head key', with a complex and real number, respectively
        """
        
        import numpy as np
        
        # Append the base to the elementlist
        self.model          = model
        model.elementlist.append(self)
        
        # Set orientation value
        self.alpha          = alpha
        
        # Set potential scaling variables
        self.head_min       = head_min
        self.head_max       = head_max
        
        # Assign the hydraulic conductivity of the base model
        self.k              = k
        
        # The model requires the base flow in terms of hydraulic potential (phi)
        # The function head_to_potential extracts the following variables:
        #   phi_min         hydraulic potential corresponding to head_min
        #   phi_max         hydraulic potential corresponding to head_max
        self.head_to_potential()
        
        # Check input for validity
        self.check_input()
        
        self.variables      = variables
        self.priors         = priors
        
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.variables  += [var]
                self.model.priors     += [self.priors[idx]]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']  
        
    def update(self):
            
        # Update potential
        self.head_to_potential()
            
    def check_input(self):
        
        # Check if phi_min is smaller than phi_offset, switch if necessary
        if self.head_min > self.head_max:
            raise Exception('Minimum potential phi_min is larger than maximum potential phi_max.')
    
    def evaluate(self,z):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Coordinates in canonical space are the start values
        z_canonical     = copy.copy(z)
        
        # head_min and head_max lie on opposite points of the circular model domain
        Q               = (self.phi_max-self.phi_min)/(self.model.domain_radius*2)
        
        # Rotate the flow field
        z               = Q*z_canonical*np.exp(-1j*self.alpha)
        
        # And offset it to phi_min
        z               = z + (self.phi_max-self.phi_min)/2 + self.phi_min
            
        return z
    
    def evaluate_gradient(self,z,derivatives = 'all'):
        
        import numpy as np
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # head_min and head_max lie on opposite points of the circular model domain
        Q               = (self.phi_max-self.phi_min)/(self.model.domain_radius*2)
        
        # Extract the derivative
        grad            = Q*np.exp(-1j*self.alpha)
            
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
            
        return grad
       
    def head_to_potential(self):
        
        for idx,h in enumerate([self.head_min-self.model.head_offset,self.head_max-self.model.head_offset]):
        
            if self.model.aquifer_type == 'confined' or (self.model.aquifer_type == 'convertible' and h >= self.model.H):
                # Strack 1989, Eq. 8.6
                pot = self.k*self.model.H*h - 0.5*self.k*self.model.H**2
                    
            elif self.model.aquifer_type == 'unconfined' or (self.model.aquifer_type == 'convertible' and h < self.model.H):
                # Strack 1989, Eq. 8.7
                pot = 0.5*self.k*h**2
                
            if idx == 0:
                self.phi_min    = pot
            elif idx == 1:
                self.phi_max    = pot
                
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,color_low = 'xkcd:cerulean',color_high='xkcd:orangish red',
             s = 50, **kwargs):
        
        import numpy as np
        import matplotlib.pyplot as plt
        
        high    = np.asarray([
            np.cos(self.alpha)*self.model.domain_radius + np.real(self.model.domain_center),
            np.sin(self.alpha)*self.model.domain_radius + np.imag(self.model.domain_center)])
        low     = -high
        
        plt.scatter(low[0],low[1],s=50,color=color_low,zorder=11,**kwargs)
        plt.scatter(high[0],high[1],s=50,color=color_high,zorder=11,**kwargs)
        
        plt.arrow(low[0]*0.9,low[1]*0.9,low[0]*0.15,low[1]*0.15,color=color_low,
                  zorder=11,head_width = 50,width = 20)
        plt.arrow(high[0]*1.1,high[1]*1.1,-high[0]*0.15,-high[1]*0.15,color=color_high,
                  zorder=11,head_width = 50,width = 20)

#%%
        
class ElementHeadBoundary:
    
    def __init__(self, model, line, line_ht, segments = None, influence = None, 
                 connectivity = 1, connectivity_normdist = None,
                 variables = [], priors=[]):
        
        """
        This implements a prescribed head boundary.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        line            - [array]   : either a real N-by-2 matrix or complex vector of length N specifying the vertices of a line string tracing the element's path
        line_ht         - [vector]  : a real vector of length N specifying the corresponding prescribed hydraulic heads at the line string's vertices
        segments        - [scalar]  : this element has a subdivision function; if a finer resolution than the number of segments in 'line' is desired, specify a larger number here; the function will then subdivide 'line' and 'line_ht' so as to create segments of as equal length as possible
        influence       - [scalar]  : radius of zero influence of each line segment; set to twice the model domain_radius if unspecified
        connectivity    - [scalar]  : either a real scalar or vector of length M, specifying if the aquifer is fully connected (1) or unconnected (0) to the HeadBoundary
        connectivity_normdist - [vector]    : if connectivity is a vector, this specifies the normalized distances 0,...,d,...,1 along which the M connectivity nodes are placed; connectivity values are then linearly interpolated for each segment

        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['line_ht']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        """

        import numpy as np
        from scipy.interpolate import interp1d
        import copy
        
        # Connect this element to the solver
        self.model          = model
        model.elementlist.append(self)
        model.linear_solver = True
        
        # Prepare the stochastic variables
        self.variables      = variables
        self.priors         = priors
        
        # Initialize the head target and connectivity variables
        self.line_ht        = line_ht
        self.connectivity   = connectivity
        if np.isscalar(self.connectivity): # Connectivity provided is uniform
        
            self.connectivity_uniform   = True
            
        else: # Connectivity provided 
        
            self.connectivity_uniform   = False
            
            # Check if normalized distances were provided
            if connectivity_normdist is None:
                raise Exception('If connectivity is not uniform, a vector of equal length containing normalized distances (e.g., [0., 0.25, 0.6, 1.]) must be specified.')
            
            # Check if connectivity_normdist is valid
            if np.min(connectivity_normdist) < 0 or np.max(connectivity_normdist) > 1:
                raise Exception('connectivity_normdist values must be between 0 and 1. Current values: '+str(connectivity_normdist))
            
            # Check if connectivity_normdist is sorted
            if not (connectivity_normdist == np.sort(connectivity_normdist)).all():
                raise Exception('connectivity_normdist values must be provided in ascending order. Current values: '+str(connectivity_normdist))
            
            self.connectivity_normdist  = connectivity_normdist
            
        # ---------------------------------------------------------------------
        # Subdivide the provided no flow boundary into #segments pieces
        
        # Complexify the line, if it wasn't already complex
        line                = self.complexify(line)
        
        # The subdivision algorith requires the line coordinates as a real N-by-2 matrix
        line                = np.column_stack((
            np.real(line)[:,np.newaxis],
            np.imag(line)[:,np.newaxis]))
        
        # Make a copy of the line
        self.line_raw       = copy.copy(line)
        
        # Check if a subdivision has been specified
        if segments is None: # No subdivision required
            self.segments   = line.shape[0]-1
        else: # Otherwise, set target
            self.segments   = segments
        
        # A number of consistency checks
        if self.segments < self.line_raw.shape[0]-1:
            raise Exception('Number of segments '+str(self.segments)+" mustn't be smaller than number of line points "+str(line.shape[0])+'.')
        if len(line_ht) != line.shape[0]:
            raise Exception('Number of head prescriptions must equal number of vertices: '+str(len(line_ht))+' =/= '+str(line.shape[0]))
        
        
        if self.segments > self.line_raw.shape[0]:
            
            # Subdivide the line
            self.line       = self.subdivide_line(np.column_stack((line,self.line_ht)),self.segments)
            self.line_c     = copy.copy(self.line[:,0] + 1j*self.line[:,1])
            self.line_ht    = copy.copy(self.line[:,2])
            
        else:
            
            # Otherwise, reconstruct the line format
            self.line       = self.line_raw.copy()
            self.line_c     = self.line[:,0] + 1j*self.line[:,1]
            self.line_ht    = line_ht
            
        # ---------------------------------------------------------------------        
        
        # Assign the initial strength variables for each segment
        self.strength       = np.ones(self.segments)
        
        # Prepare the influence range for this line sink
        if influence is None:
            # If no influence range is specified, set it to twice the domain radius
            # to ensure that no point in the model domain will lie outside this range
            self.influence  = self.model.domain_radius*2
        else:
            self.influence  = influence
        
        # Prepare a few variables for this element
        self.L              = []    # Length of each line segment
        self.zc             = []    # Center of each line segment
        self.head_target    = []    # Head target at each line segment
        
        for seg in range(self.segments):
            
            self.L              += [np.abs(self.line_c[seg+1] - self.line_c[seg])]
            self.zc             += [(self.line_c[seg]+self.line_c[seg+1])/2]
            self.head_target    += [(self.line_ht[seg]+self.line_ht[seg+1])/2]
            
        # Convert list of segment centers to array
        self.zc             = np.asarray(self.zc)
        self.head_target    = np.asarray(self.head_target)
        
        # Now form a vector of cumulative distances
        self.cumdist        = []
        for seg in range(self.segments):
            if seg == 0:
                self.cumdist.append(np.abs(self.zc[0]-self.line_c[0]))
            else:
                self.cumdist.append(np.abs(self.zc[seg]-self.zc[seg-1]))
        self.cumdist        = np.cumsum(np.asarray(self.cumdist))
        self.cumdist        /= (self.cumdist[-1] + np.abs(self.zc[-1]-self.line_c[-1]))
        
        if not self.connectivity_uniform:
            
            # Interpolate the connectivity
            from scipy.interpolate import interp1d
            itp     = interp1d(self.connectivity_normdist,self.connectivity)
            self.connectivity_interpolated  = itp(self.cumdist)
        
        
        
        
        
        
        
            
        # Convert the head targets to potential targets
        self.set_potential_target()
        
        # Check if the prior matches the number of parameters
        if len(self.priors) != len(self.variables):
            raise Exception('Number of priors must match number of unknown variables. Number of priors: '+str(self.priors)+' / Number of unknown variables: '+str(len(self.variables)))
        
        # Go through all elements
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.priors     += [self.priors[idx]]
                self.model.variables  += [var]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']   
                       
    def update(self):
        
        import numpy as np
        
        # Update the potential targets
        self.set_potential_target()
        
        if not self.connectivity_uniform:
            
            # Interpolate the connectivity
            from scipy.interpolate import interp1d
            itp     = interp1d(self.connectivity_normdist,self.connectivity)
            self.connectivity_interpolated  = itp(self.cumdist)
        
        # self.zc             = self.xc + 1j*self.yc
        # self.L              = np.abs(self.z2 - self.z1)
        
        # influence_pt        = (self.z2-self.z1)*self.influence/self.L + self.z1
        # Z                   = (2*influence_pt-(self.z1+self.z2))/(self.z2-self.z1)
        # part1               = np.nan_to_num((Z+1)*np.log(Z+1))
        # part2               = np.nan_to_num((Z-1)*np.log(Z-1))
        # self.offset_outside = self.L / (4*np.pi) * (part1 - part2)

    def evaluate_gradient(self,z,detailed = False, derivatives = 'all', override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # 'detailed' returns the results as a matrix instead of a summed vector
        if detailed:
            grad    = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            grad    = np.zeros(z.shape, dtype = np.complex)
        
        for seg in range(self.segments):
        
            Z       = (2*z-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg])
            
            Z[np.where(np.abs(np.imag(Z)) < 1E-10)] = np.real(Z[np.where(np.abs(np.imag(Z)) < 1E-10)])
            
            # Now get the gradient d omega(z)/dZ
            if self.connectivity_uniform:
                # If the connectivity is uniform, i.e. does not vary along the boundary
                if not override_parameters:
                    temp    = self.strength[seg]*self.connectivity*self.L[seg]/4/np.pi*(np.log(Z+1) - np.log(Z-1))
                else:
                    temp    = self.L[seg]*self.connectivity/4/np.pi*(np.log(Z+1) - np.log(Z-1))
            else:
                # If the connectivity is not uniform, i.e. does vary along the boundary
                if not override_parameters:
                    temp    = self.strength[seg]*self.connectivity_interpolated[seg]*self.L[seg]/4/np.pi*(np.log(Z+1) - np.log(Z-1))
                else:
                    temp    = self.L[seg]*self.connectivity_interpolated[seg]/4/np.pi*(np.log(Z+1) - np.log(Z-1))
                    
            # To get d omega(z)/dz we can use the product rule
            #   d omega(z)/dz = d omega(z)/dZ * dZ/dz
            # hence:
            temp    = temp*2/(self.line_c[seg+1]-self.line_c[seg])
            
            if detailed:
                grad[seg,:]     = copy.copy(temp)
            else:
                grad    += temp
        
        
        
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
        
        
        return grad

    def evaluate(self,z,detailed = False, override_parameters = False,
                 evaluate_self = False):
        
        import copy
        import numpy as np
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        if detailed:
            res     = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            res     = np.zeros(z.shape, dtype = np.complex)
        
        for seg in range(self.segments):
            
            # Convert to local coordinates
            Z       = (2*z-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg])
            
            Z[np.where(np.abs(np.imag(Z)) < 1E-10)] = np.real(Z[np.where(np.abs(np.imag(Z)) < 1E-10)])
            
            if self.connectivity_uniform:
                # If the connectivity is uniform, i.e. does not vary along the boundary
                # Evaluate the complex potential offset by a distance in the 
                if not override_parameters:
                    temp    = self.strength[seg]*self.connectivity*self.L[seg]/4/np.pi * (\
                              (Z+1)*np.log(Z+1) - \
                              (Z-1)*np.log(Z-1) - \
                              (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                              (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
                elif override_parameters and evaluate_self:
                    # This is a unique addition. In order to make the connectivity 
                    # work as an element in conjunction with all others, we must 
                    # 'trick' it into thinking its connectivity is 1. This is only
                    # ever activated to evaluate the effects of a prescribed head
                    # boundary on itself (a diagonal block in setting up the matrix
                    # for the linear system of equations).
                    temp    = self.L[seg]/4/np.pi * (\
                              (Z+1)*np.log(Z+1) - \
                              (Z-1)*np.log(Z-1) - \
                              (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                              (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
                else:
                    temp    = self.L[seg]*self.connectivity/4/np.pi * (\
                              (Z+1)*np.log(Z+1) - \
                              (Z-1)*np.log(Z-1) - \
                              (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                              (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
            
            else:
                # If the connectivity is not uniform, i.e. does vary along the boundary
                
                # Evaluate the complex potential offset by a distance in the 
                if not override_parameters:
                    temp    = self.strength[seg]*self.connectivity_interpolated[seg]*self.L[seg]/4/np.pi * (\
                              (Z+1)*np.log(Z+1) - \
                              (Z-1)*np.log(Z-1) - \
                              (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                              (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
                elif override_parameters and evaluate_self:
                    # This is a unique addition. In order to make the connectivity 
                    # work as an element in conjunction with all others, we must 
                    # 'trick' it into thinking its connectivity is 1. This is only
                    # ever activated to evaluate the effects of a prescribed head
                    # boundary on itself (a diagonal block in setting up the matrix
                    # for the linear system of equations).
                    temp    = self.L[seg]/4/np.pi * (\
                              (Z+1)*np.log(Z+1) - \
                              (Z-1)*np.log(Z-1) - \
                              (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                              (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
                else:
                    temp    = self.L[seg]*self.connectivity_interpolated[seg]/4/np.pi * (\
                              (Z+1)*np.log(Z+1) - \
                              (Z-1)*np.log(Z-1) - \
                              (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                              (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
            
            
            # If evaluated directly at the endpoints, the result would be NaN
            # They should be zero, see Bakker 2009
            temp    = np.nan_to_num(temp)
            
            if detailed:
                res[seg,:]      = copy.copy(temp)
            else:
                res             += temp
            
        return res
    
    def subdivide_line(self,line,segments):
        
        import numpy as np
        import copy
        
        # If array is one-dimensional, reshape it appropriately
        if len(line.shape) == 1: line = line.reshape((line.shape[0],1))
        
        D   = line.shape[1]
        
        # Calculate the lengths of original segments
        length  = [np.linalg.norm(line[seg,:]-line[seg+1,:]) for seg in range(line.shape[0]-1)]
        
        # Normalize the length of the original segments
        length  /= np.sum(length)
        
        # Calculate the number of new segments we must create, the line already has
        # (#vertices-1) segments. We only require the difference
        new_segments = segments - line.shape[0] + 1
        
        # Calculate where those segments should go
        bins    = np.concatenate(( [0] , np.cumsum(length) ))
        
        # Add Extend the bin length a bit to prevent errors from arithmetic under- or overflow
        bins[0]     -= 1E-10
        bins[-1]    += 1E-10
        
        # Distribute vertices along the segments
        x       = np.linspace(0,1,new_segments+1)
        num_vertices = []
        for seg in range(line.shape[0]-1):
            if seg == 0:
                num_vertices += [len(np.where(x <= bins[1])[0])]
            else:
                num_vertices += [len(np.where(np.logical_and(
                        x > bins[seg],
                        x <= bins[seg+1]))[0])]
    
        # Subidivide the original segments
        new_vertices = []
        for seg in range(line.shape[0]-1):
            
            temp = None
            for d in range(D):
                if temp is None:
                    temp    = copy.copy(line[seg,d] + (line[seg+1,d]-line[seg,d]) * np.linspace(0,1,num_vertices[seg]+2)[1:-1])
                else:
                    temp    = np.column_stack((
                        temp,
                        copy.copy(line[seg,d] + (line[seg+1,d]-line[seg,d]) * np.linspace(0,1,num_vertices[seg]+2)[1:-1])))
                
            new_vertices += [copy.copy(temp)]
        
        # Create the seed for the new line
        new_line    = copy.copy(line[0,:].reshape((1,D)) )
        
        # Assemble the new line
        for seg in range(line.shape[0]-1):
            
            # Add the new segments, then the next original vertex
            new_line    = np.row_stack((
                new_line,
                new_vertices[seg],
                line[seg+1,:]))
        
        return new_line

    def set_potential_target(self):
        
        import copy
        import numpy as np
        
        # Get the hydraulic conductivities at the segment control points
        for e in self.model.elementlist:
            if isinstance(e, ElementMoebiusBase) or isinstance(e, ElementUniformBase):
                temp_k = np.ones(self.zc.shape)*e.k
        for e in self.model.elementlist:
            if isinstance(e, ElementInhomogeneity):
                inside  = e.are_points_inside_polygon(self.zc)
                temp_k[inside] = e.k
        
        # Create a list of hydraulic potential targets
        self.phi_target = copy.copy(self.head_target - self.model.head_offset)
        if self.model.aquifer_type == 'confined':
            # Strack 1989, Eq. 8.6
            self.phi_target = temp_k*self.model.H*self.phi_target - \
                0.5*temp_k*self.model.H**2
        elif self.model.aquifer_type == 'unconfined':
            # Strack 1989, Eq. 8.7
            self.phi_target = 0.5*temp_k*self.phi_target**2
        elif self.model.aquifer_type == 'convertible':
            # Find out which points are confined and which are unconfined
            index_conf      = np.where(self.phi_target >= self.model.H)[0]
            index_unconf    = np.where(self.phi_target < self.model.H)[0]
            # Account for the confined points
            # confined:     Strack 1989, Eq. 8.6
            self.phi_target[index_conf] = \
                temp_k[index_conf]*self.model.H*self.phi_target[index_conf] - \
                0.5*temp_k[index_conf]*self.model.H**2
            # unconfined:   Strack 1989, Eq. 8.7
            self.phi_target[index_unconf] = \
                0.5*temp_k[index_unconf]*self.phi_target[index_unconf]**2
    
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,col='xkcd:kermit green',zorder=12,linewidth=5,**kwargs):
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        plt.plot(np.real(self.line_c),np.imag(self.line_c),color=col,
                 zorder=zorder,linewidth=linewidth,**kwargs)

#%%

class ElementNoFlowBoundary:
    
    def __init__(self, model, line, segments = None,head_target = 0,
                 variables = [], priors=[]):
        
        """
        This implements a no-flow boundary.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        line            - [array]   : either a real N-by-2 matrix or complex vector of length N specifying the vertices of a line string tracing the element's path
        segments        - [scalar]  : this element has a subdivision function; if a finer resolution than the number of segments in 'line' is desired, specify a larger number here; the function will then subdivide 'line' and 'line_ht' so as to create segments of as equal length as possible

        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['line_ht']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        """

        import numpy as np
        from scipy.interpolate import interp1d
        import copy
        
        # Append this element to the specified model
        self.model          = model
        model.elementlist.append(self)
        model.linear_solver = True

        # ---------------------------------------------------------------------
        # Subdivide the provided no flow boundary into segments pieces
        
        # Complexify the line, if it wasn't already complex
        line                = self.complexify(line)
        
        # The subdivision algorith requires the line coordinates as a real N-by-2 matrix
        line                = np.column_stack((
            np.real(line)[:,np.newaxis],
            np.imag(line)[:,np.newaxis]))
        
        self.line_raw       = copy.copy(line)
        if segments is None:
            self.segments   = line.shape[0]-1
        else:
            self.segments   = segments
        
        if self.segments < self.line_raw.shape[0]-1:
            raise Exception('Prescribed number of line segments '+str(self.segments)+" mustn't be smaller than base number of segments "+str(line.shape[0]-1)+'.')
        
        if self.segments > self.line_raw.shape[0]-1:
            
            # Subdivide the line
            self.line   = self.subdivide_line(line,self.segments)
            self.line_c = self.line[:,0] + 1j*self.line[:,1]
            
        else:
            
            self.line       = self.line_raw.copy()
            self.line_c     = self.line[:,0] + 1j*self.line[:,1]
            
        # Also get the normal vector components to each segment
        self.line_nvec  = self.line[:,1] - 1j*self.line[:,0]
        self.line_nvec  = self.line_nvec/np.abs(self.line_nvec)

        # ---------------------------------------------------------------------
        
        # Get strength parameters for each vertex
        self.strength       = np.ones(self.segments)
        
        
        self.zc             = []
        self.segment_nvec   = []
        self.L              = []
        
        for seg in range(self.segments):
            
            self.zc             += [(self.line_c[seg]+self.line_c[seg+1])/2]
            
            # Calculate the normal vector to this segment
            self.segment_nvec   += [(self.line_c[seg]-self.line_c[seg+1])]
            self.segment_nvec[-1]= [np.imag(self.segment_nvec[-1])-1j*np.real(self.segment_nvec[-1])]
            
            self.L              += [np.abs(self.line_c[seg+1] - self.line_c[seg])]
            
        self.zc             = np.asarray(self.zc)
        
        # Extract target variables
        self.variables      = variables
        self.priors         = priors
        
        self.L              = np.asarray(self.L)
        
        
        # Check if the prior matches the number of parameters
        if len(self.priors) != len(self.variables):
            raise Exception('Number of priors must match number of unknown variables. Number of priors: '+str(self.priors)+' / Number of unknown variables: '+str(len(self.variables)))
        
        # Go through all elements
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.priors     += [self.priors[idx]]
                self.model.variables  += [var]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']   
                    
    def update(self):
        
        import numpy as np

    def evaluate_gradient(self,z,detailed = False,derivatives = 'all',override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # 'detailed' returns the results as a matrix instead of a summed vector
        if detailed:
            grad    = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            grad    = np.zeros(z.shape, dtype = np.complex)
        
        # Go through all line segments
        for seg in range(self.segments):
        
            # Convert z to the local variable Z
            Z       = (2*copy.copy(z)-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg]) #-marked- last influence
            
            # Add to the result file d Omega(Z) / dZ
            if not override_parameters:
                temp    = 1j*self.strength[seg]/(np.pi-np.pi*Z**2)
            else:
                temp    = 1j/(np.pi-np.pi*Z**2)
            
            # Multiply the result with dZ/dz to obtain dOmega(Z)/dz
            temp    = temp*2/(self.line_c[seg+1]-self.line_c[seg])
            
            if detailed:
                grad[seg,:]     = copy.copy(temp)
            else:
                grad            += copy.copy(temp)
        
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
        
        return grad
            
    def evaluate(self,z,detailed = False,override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # 'detailed' returns the results as a matrix instead of a summed vector
        if detailed:
            res     = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            res     = np.zeros(z.shape, dtype = np.complex)
        
        # Go through all line segments
        for seg in range(self.segments):
            
            # Convert z to the local variable Z
            Z       = (2*copy.copy(z)-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg]) #-marked- last influence

            # # If any Z values are at zero, offset them by a bit
            # indices = np.where(np.abs(Z) < 1E-10)
            
            # Z[indices]  = 1E-10
            
            
            indices     = np.where(np.logical_and(
                np.abs(Z-1) > 1E-10,
                np.abs(Z+1) > 1E-10) )[0]
            indices     = np.ones(Z.shape,dtype=bool)
            # indices_minus_1 = np.where(np.abs(Z+1) > 1E-10)[0]
            
            # Term 1
            term_1  = np.zeros(Z.shape,dtype=np.complex)
            term_1[indices]  = (Z[indices]+1)*np.log((Z[indices]-1)/(Z[indices]+1))
            
            term_2  = np.zeros(Z.shape,dtype=np.complex)
            term_2[indices] = (Z[indices]-1)*np.log((Z[indices]-1)/(Z[indices]+1))
            
            term_1  = np.nan_to_num(term_1)
            term_2  = np.nan_to_num(term_2)
            
            
            # indices_plus_1  = np.where(np.abs(Z-1) > 1E-10)[0]
            # indices_minus_1 = np.where(np.abs(Z+1) > 1E-10)[0]
            
            # # Term 1
            # term_1  = np.zeros(Z.shape,dtype=np.complex)
            # term_1[indices_plus_1]  = (Z[indices_plus_1]+1)*np.log((Z[indices_plus_1]-1)/(Z[indices_plus_1]+1))
            
            # term_2  = np.zeros(Z.shape,dtype=np.complex)
            # term_2[indices_minus_1] = (Z[indices_minus_1]-1)*np.log((Z[indices_minus_1]-1)/(Z[indices_minus_1]+1))

            # Blablablub
            if not override_parameters:
                temp        = self.strength[seg]/(4*np.pi*1j) * \
                    (term_1 - term_2)
            else:
                temp        = 1/(4*np.pi*1j) * \
                    (term_1 - term_2)

            # # Blablablub
            # if not override_parameters:
            #     temp        = self.strength[seg]/(4*np.pi*1j) * \
            #         ((Z+1)*np.log((Z-1)/(Z+1)) - (Z-1)*np.log((Z-1)/(Z+1)))
            # else:
            #     temp        = 1/(4*np.pi*1j) * \
            #         ((Z+1)*np.log((Z-1)/(Z+1)) - (Z-1)*np.log((Z-1)/(Z+1)))
            
            if detailed:
                res[seg,:]      = copy.copy(temp)
            else:
                res             += copy.copy(temp)
            
        return res
    
    def subdivide_line(self,line,segments):
        
        import numpy as np
        import copy
        
        # If array is one-dimensional, reshape it appropriately
        if len(line.shape) == 1: line = line.reshape((line.shape[0],1))
        
        D   = line.shape[1]
        
        # Calculate the lengths of original segments
        length  = [np.linalg.norm(line[seg,:]-line[seg+1,:]) for seg in range(line.shape[0]-1)]
        
        # Normalize the length of the original segments
        length  /= np.sum(length)
        
        # Calculate the number of new segments we must create, the line already has
        # (#vertices-1) segments. We only require the difference
        new_segments = segments - line.shape[0] + 1
        
        # Calculate where those segments should go
        bins    = np.concatenate(( [0] , np.cumsum(length) ))
        
        # Add Extend the bin length a bit to prevent errors from arithmetic under- or overflow
        bins[0]     -= 1E-10
        bins[-1]    += 1E-10
        
        # Distribute vertices along the segments
        x       = np.linspace(0,1,new_segments)
        num_vertices = []
        for seg in range(line.shape[0]-1):
            if seg == 0:
                num_vertices += [len(np.where(x <= bins[1])[0])]
            else:
                num_vertices += [len(np.where(np.logical_and(
                        x > bins[seg],
                        x <= bins[seg+1]))[0])]
    
        # Subidivide the original segments
        new_vertices = []
        for seg in range(line.shape[0]-1):
            
            temp = None
            for d in range(D):
                if temp is None:
                    temp    = copy.copy(line[seg,d] + (line[seg+1,d]-line[seg,d]) * np.linspace(0,1,num_vertices[seg]+2)[1:-1])
                else:
                    temp    = np.column_stack((
                        temp,
                        copy.copy(line[seg,d] + (line[seg+1,d]-line[seg,d]) * np.linspace(0,1,num_vertices[seg]+2)[1:-1])))
                
            new_vertices += [copy.copy(temp)]
        
        # Create the seed for the new line
        new_line    = copy.copy(line[0,:].reshape((1,D)) )
        
        # Assemble the new line
        for seg in range(line.shape[0]-1):
            
            # Add the new segments, then the next original vertex
            new_line    = np.row_stack((
                new_line,
                new_vertices[seg],
                line[seg+1,:]))
        
        return new_line
    
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,color='xkcd:dark grey',linewidth=5,zorder=10,**kwargs):
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        plt.plot(np.real(self.line_c),np.imag(self.line_c),color=color,
                 linewidth = linewidth, zorder = zorder, **kwargs)
    
#%%

class ElementInhomogeneity:
    
    def __init__(self, model, polygon, segments = None, k = 0.1,
                 variables = [], priors=[], snap_distance = 1E-10,
                 zero_cutoff = 1E-10, snap = True):

        """
        This implements a zonal hydraulic conductivity inhomogeneity.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        polygon         - [array]   : either a real N-by-2 matrix or complex vector of length N specifying the vertices of a polygon tracing the element's shape
        segments        - [scalar]  : this element has a subdivision function; if a finer resolution than the number of segments in 'polygon' is desired, specify a larger number here; the function will then subdivide 'line' and 'line_ht' so as to create segments of as equal length as possible
        k               - [scalar]  : hydraulic conductivity inside the inhomogeneity in canonical units (e.g., 1E-5 [length units]/[time units])

        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['line_ht']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        """

        import numpy as np
        import copy
        import matplotlib.path
        
        # Append this element to the specified model
        self.model          = model
        model.elementlist.append(self)
        model.linear_solver = True
        
        # Prepare the polygon variable
        self.polygon        = polygon
        
        self.polygon        = self.complexify(self.polygon)
        
        self.snap_distance  = snap_distance
        self.zero_cutoff    = zero_cutoff
        
        # Is the polygon closed? If not, close it temporarily
        if np.abs(self.polygon[0]-self.polygon[-1]) > self.snap_distance:
            self.polygon    = np.asarray(list(self.polygon)+[self.polygon[0]])
        
        # Also create an array with real coordinates
        self.polygon_XY     = np.column_stack((
            np.real(copy.copy(self.polygon))[:,np.newaxis],
            np.imag(copy.copy(self.polygon))[:,np.newaxis] ))

        # Is the polygon counter-clockwise? If not, correct it
        if self.are_vertices_clockwise(self.polygon_XY):
            self.polygon    = np.flip(self.polygon)
            self.polygon_XY = np.flipud(self.polygon_XY)
            
        # Do we wish to subdivide the polygon?
        # First, check if the user specified a desired segment count
        if segments is None:
            self.segments   = self.polygon.shape[0]-1
        else:
            self.segments   = segments
        
        if self.segments < self.polygon.shape[0]-1:
            raise Exception('Prescribed number of line segments '+str(self.segments)+" mustn't be smaller than the number of vertices "+str(polygon.shape[0]-1)+'.')
        
        # Subdivide the polygon, if desired
        if self.segments > self.polygon.shape[0]-1:
            self.polygon_XY = self.subdivide_line(self.polygon_XY,self.segments)
            self.polygon    = self.polygon_XY[:,0] + 1j*self.polygon_XY[:,1]
            
        # Un-close the polygon again
        self.polygon_XY     = self.polygon_XY[:-1,:]
        self.polygon        = self.polygon[:-1]
            
        # If vertex snapping is enabled, snap all outside vertices onto the domain edge
        if snap:
            self.snap_to_domain()

        # This is a hack: We shrink the polygon by a small amount. This ensures 
        # that no issues arise from evaluating points directly on the boundary, 
        # and allows us to consider inhomogeneities directly bounding each other; 
        # there might be other ways to solve this issue alternatively
        self.polygon_XY = self.shrink_polygon(
            polygon = self.polygon_XY,
            offset  = 1E-10)
        self.polygon    = self.polygon_XY[:,0] + 1j*self.polygon_XY[:,1]

        # The control points of the inhomogeneity are simply its vertices
        # This is required for the linear solver
        self.zc             = self.polygon
        
        # Raise an exception if this inhomogeneity intersects any of the previous
        # inhomogeneities
        for e in self.model.elementlist[:-1]:
            if isinstance(e, ElementInhomogeneity):
                if any(e.are_points_inside_polygon(self.zc)):
                    raise Exception('Inhomogeneities may not intersect each other.')
        
        # Create a path with the edges of the polygon
        # We can use this path to find out if evaluation points are inside or 
        # or outside the inhomogeneity.
        self.linepath   = matplotlib.path.Path(self.polygon_XY)
        
        # Get strength parameters for each vertex
        self.strength       = np.ones(self.segments)
        
        # Assign the hydraulic conductivity of the inhomogeneity
        self.k              = k
        
        # Extract target variables
        self.variables      = variables
        self.priors         = priors
        
        # Prepare the matrix block containing the effect of this element onto 
        # itself for future use in solving the linear system. The matrix requires
        # subtraction of the A_star variable from its diagonal entries for completion
        self.block          = self.matrix_contribution()
        
        # Check if the prior matches the number of parameters
        if len(self.priors) != len(self.variables):
            raise Exception('Number of priors must match number of unknown variables. Number of priors: '+str(self.priors)+' / Number of unknown variables: '+str(len(self.variables)))
        
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.priors     += [self.priors[idx]]
                self.model.variables  += [var]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']   
                    
    def update(self):
        
        import numpy as np
        
    def evaluate_gradient(self,z,detailed = False,derivatives = 'all',override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # 'detailed' returns the results as a matrix instead of a summed vector
        if detailed:
            grad    = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            grad    = np.zeros(z.shape, dtype = np.complex)
        
        # Go through all line segments
        for seg in range(self.segments):
            
            # Find the next (seg_plus) and previous (seg_minus) vertex
            if seg == self.segments-1:
                seg_plus    = 0
                seg_minus   = seg-1
            else:
                seg_plus    = seg+1
                seg_minus   = seg-1
        
            if override_parameters:
                
                # Calculate the gradient
                temp    = 1j/(2*np.pi)* (\
                    np.log((self.polygon[seg] - z)/(self.polygon[seg_minus]-z)) / \
                        (self.polygon[seg_minus]-self.polygon[seg]) - \
                    np.log((self.polygon[seg_plus] - z)/(self.polygon[seg]-z)) / \
                        (self.polygon[seg]-self.polygon[seg_plus]) )
                    
            else:
                
                # Calculate the gradient
                temp    = self.strength[seg]*1j/(2*np.pi)* (\
                    np.log((self.polygon[seg] - z)/(self.polygon[seg_minus]-z)) / \
                        (self.polygon[seg_minus]-self.polygon[seg]) - \
                    np.log((self.polygon[seg_plus] - z)/(self.polygon[seg]-z)) / \
                        (self.polygon[seg]-self.polygon[seg_plus]) )
                
            if detailed:
                grad[seg,:]     = copy.copy(temp)
            else:
                grad            += copy.copy(temp)
        
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
        
        return grad
            
    def evaluate(self,z,detailed = False,override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Define the segment influence functions
        def F(Z):
            return np.nan_to_num(-0.5*(Z-1)*np.log((Z-1)/(Z+1))) - 1
        def G(Z):
            return np.nan_to_num(+0.5*(Z+1)*np.log((Z-1)/(Z+1))) + 1
        
        # 'detailed' returns the results as a matrix instead of a summed vector
        if detailed:
            res     = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            res     = np.zeros(z.shape, dtype = np.complex)
        
        # Go through all line segments
        for seg in range(self.segments):
            
            temp    = np.zeros(z.shape,     dtype=np.complex)
            
            # Find the next (seg_plus) and previous (seg_minus) vertex
            if seg == self.segments-1:
                seg_plus    = 0
                seg_minus   = seg-1
            else:
                seg_plus    = seg+1
                seg_minus   = seg-1
            
            # Get a vector of distances
            dist    = np.abs(z - self.polygon[seg])
        
            if override_parameters:
                
                # First, use the standard solution for points which aren't snapped to vertices
                Z_before = \
                    (2*z - (self.polygon[seg_minus] + self.polygon[seg]))/(self.polygon[seg] - self.polygon[seg_minus])
                Z_after = \
                    (2*z - (self.polygon[seg] + self.polygon[seg_plus]))/(self.polygon[seg_plus] - self.polygon[seg])
                
                # # Prevent errors from underflow
                # Z_before[np.where(np.abs(np.imag(Z_before)) < 1E-10)]   = np.real(Z_before[np.where(np.abs(np.imag(Z_before)) < 1E-10)])
                # Z_after[np.where(np.abs(np.imag(Z_after)) < 1E-10)]     = np.real(Z_after[np.where(np.abs(np.imag(Z_after)) < 1E-10)])
                
                # Z_before[np.where(np.abs(np.real(Z_before)) < 1E-10)]   = 1j*np.imag(Z_before[np.where(np.abs(np.real(Z_before)) < 1E-10)])
                # Z_after[np.where(np.abs(np.real(Z_after)) < 1E-10)]     = 1j*np.imag(Z_after[np.where(np.abs(np.real(Z_after)) < 1E-10)])

                temp  = 1/(2*np.pi*1j)*(G(Z_before)+F(Z_after))
    
            else:
            
                # First, use the standard solution for points which aren't snapped to vertices
                Z_before = \
                    (2*z - (self.polygon[seg_minus] + self.polygon[seg]))/(self.polygon[seg] - self.polygon[seg_minus])
                Z_after = \
                    (2*z - (self.polygon[seg] + self.polygon[seg_plus]))/(self.polygon[seg_plus] - self.polygon[seg])
                
                # # Prevent errors from underflow
                # Z_before[np.where(np.abs(np.imag(Z_before)) < 1E-10)]   = np.real(Z_before[np.where(np.abs(np.imag(Z_before)) < 1E-10)])
                # Z_after[np.where(np.abs(np.imag(Z_after)) < 1E-10)]     = np.real(Z_after[np.where(np.abs(np.imag(Z_after)) < 1E-10)])
                
                # Z_before[np.where(np.abs(np.real(Z_before)) < 1E-10)]   = 1j*np.imag(Z_before[np.where(np.abs(np.real(Z_before)) < 1E-10)])
                # Z_after[np.where(np.abs(np.real(Z_after)) < 1E-10)]     = 1j*np.imag(Z_after[np.where(np.abs(np.real(Z_after)) < 1E-10)])
                
                temp  = self.strength[seg]/(2*np.pi*1j)*(G(Z_before)+F(Z_after))

            if detailed:
                res[seg,:]      = copy.copy(temp)
            else:
                res             += copy.copy(temp)
                
            self.res    = res
            
        return res
    
    def matrix_contribution(self):
        
        """
        This function writes a block into the matrix for the solution of the system
        of linear equations. It only evaluates the influence of the element itself
        onto itself. For its influence on other inhomogeneities, use the evaluate
        function with detailed = True.
        """
        
        import numpy as np
        import copy
        
        # The functions F and G sometimes return NaN, errors we catch through
        # np.nan_to_num. Suppress these error messages.
        import warnings
        warnings.filterwarnings('ignore')
        
        # Define the segment influence functions
        def F(Z):
            return np.nan_to_num(-0.5*(Z-1)*np.log((Z-1)/(Z+1))) - 1
        def G(Z):
            return np.nan_to_num(+0.5*(Z+1)*np.log((Z-1)/(Z+1))) + 1
        
        # We evaluate this block at its own vertices
        z   = self.polygon
            
        # Pre-allocate an empty matrix for the block
        block   = np.zeros((self.segments,self.segments))
        
        self.angles = []
        self.temp   = []
        
        # Go through all vertices in the polygon
        for seg in range(self.segments):
            
            # Set the previous, current, and next vertex of the polygon
            if seg == self.segments-1:
                seg_minus   = seg-1
                seg_center  = seg
                seg_plus    = 0
            else:
                seg_minus   = seg-1
                seg_center  = seg
                seg_plus    = seg+1
                
            self.temp.append([self.polygon[seg_plus]-self.polygon[seg_center],
                              self.polygon[seg_center]-self.polygon[seg_minus]])
            
            # To evaluate the effect of a vertex on itself, it is computationally
            # cleanest to evaluate it in terms of angles; these angles are 
            # calculated according to Strack 1989, Eq. 35.29 and 35.30
            newtemp     = np.angle(self.polygon[seg_plus]-self.polygon[seg_center]) - \
                          np.angle(self.polygon[seg_center]-self.polygon[seg_minus])
                          
            if newtemp < -np.pi:  newtemp += 2*np.pi
            if newtemp > +np.pi:  newtemp -= 2*np.pi
                          
            # self.angles.append(newtemp)
            # Sometimes, numerical imprecision causes the angle to fall outside
            # the range 0 and 2 pi; in that case, flip it back inside
            # if newtemp < 0:         newtemp += 2*np.pi
            # if newtemp > 2*np.pi:   newtemp -= 2*np.pi
            newtemp     -= np.pi
            

            
            self.angles.append(newtemp)
            
            # Write the diagonal entries of the matrix
            block[seg,seg] = 1/(2*np.pi)*newtemp
            
            # Here we would normally add the factor for the conductivity difference
            # to the diagonal entries; since we only prepare the matrix here, 
            # we skip it
            
            # # Determine the A_star variable (Strack 1989 35.4, 35.38)
            # A_star          = self.model.k/(self.k - self.model.k)
            # block[seg,seg]  -= A_star

            # Then handle all off-diagonal contributions
            for seg2 in range(self.segments):
                
                # Skip the diagonal
                if seg2 != seg:
                    
                    # Get the indices of the past, current, and next vertex
                    if seg2 == self.segments-1:
                        seg_minus   = seg2-1
                        seg_center  = seg2
                        seg_plus    = 0
                    else:
                        seg_minus   = seg2-1
                        seg_center  = seg2
                        seg_plus    = seg2+1
                    
                    # Calculate the local coordinates
                    Z_before = \
                        (2*z[seg] - (self.polygon[seg_minus] + self.polygon[seg_center]))/(self.polygon[seg_center] - self.polygon[seg_minus])
                    Z_after = \
                        (2*z[seg] - (self.polygon[seg_center] + self.polygon[seg_plus]))/(self.polygon[seg_plus] - self.polygon[seg_center])
                    
                    # And write the result into the correct matrix entries
                    block[seg,seg2] = copy.copy(np.real(1/(2*np.pi*1j)*(G(Z_before)+F(Z_after))))
                    
        return block
        
    def are_vertices_clockwise(self,line):
        
        """
        This function takes an string of 2-D vertices of a polygon, provided as a 
        N x 2 numpy array, and returns a boolean specifying whether the vertices 
        are provided in clock-wise or counter-clock-wise order.
        
        Parameters:
            
            line            - Required  : numpy array of polygon vertices (N by 2)
        """
        
        import numpy as np
        
        signed_area = 0
        for idx in range(line.shape[0]):
        
            x1 = line[idx,0]
            y1 = line[idx,1]
            if idx == line.shape[0]-1:
                x2 = line[0,0]
                y2 = line[0,1]
            else:
                x2 = line[idx+1,0]
                y2 = line[idx+1,1]
        
            signed_area += (x1 * y2 - x2 * y1)
        
        return (np.sign(signed_area) == -1.)


    def are_points_inside_polygon(self,z):
        
        import matplotlib.path
        import numpy as np
        
        indices = self.linepath.contains_points(
            np.column_stack((
                np.real(z),
                np.imag(z))))
        
        return indices
    
    def are_points_on_polygon(self,points,line = None,snap_distance = 1E-10):
        
        import numpy as np
        
        points  = np.column_stack((
            np.real(points)[:,np.newaxis],
            np.imag(points)[:,np.newaxis]))
        
        if line is None:
            # Get the polygon and close it
            line    = np.row_stack((self.polygon_XY,self.polygon_XY[0,:]))
        
        # Pre-allocate space for the indices
        indices     = np.zeros(points.shape[0],dtype=np.bool)
        
        # Go through all line segments
        for n in range(line.shape[0]-1):
            
            a   = line[n,:].reshape((1,2))
            b   = line[n+1,:].reshape((1,2))
            
            # Normalized tangent vectors
            d_ba = b - a
            d = np.divide(d_ba, (np.hypot(d_ba[:, 0], d_ba[:, 1]).reshape(-1, 1)))
        
            # Signed parallel distance components
            # Row-wise dot products of 2D vectors
            s = np.multiply(a - points, d).sum(axis=1)
            t = np.multiply(points - b, d).sum(axis=1)
        
            # Clamped parallel distance
            h = np.maximum.reduce([s, t, np.zeros(len(s))])
        
            # Perpendicular distance component
            # Row-wise cross products of 2D vectors  
            d_pa = points - a
            c = d_pa[:, 0] * d[:, 1] - d_pa[:, 1] * d[:, 0]
            
            # Calculate the distance
            d   = np.hypot(h, c)
            
            indices[np.where(d <= snap_distance)] = True
            
        return indices
    
    def subdivide_line(self,line,segments):
        
        import numpy as np
        import copy
        
        # If array is one-dimensional, reshape it appropriately
        if len(line.shape) == 1: line = line.reshape((line.shape[0],1))
        
        D   = line.shape[1]
        
        # Calculate the lengths of original segments
        length  = [np.linalg.norm(line[seg,:]-line[seg+1,:]) for seg in range(line.shape[0]-1)]
        
        # Normalize the length of the original segments
        length  /= np.sum(length)
        
        # Calculate the number of new segments we must create, the line already has
        # (#vertices-1) segments. We only require the difference
        new_segments = segments - line.shape[0] + 1
        
        # Calculate where those segments should go
        bins    = np.concatenate(( [0] , np.cumsum(length) ))
        
        # Add Extend the bin length a bit to prevent errors from arithmetic under- or overflow
        bins[0]     -= 1E-10
        bins[-1]    += 1E-10
        
        # Distribute vertices along the segments
        x       = np.linspace(0,1,new_segments)
        num_vertices = []
        for seg in range(line.shape[0]-1):
            if seg == 0:
                num_vertices += [len(np.where(x <= bins[1])[0])]
            else:
                num_vertices += [len(np.where(np.logical_and(
                        x > bins[seg],
                        x <= bins[seg+1]))[0])]
    
        # Subidivide the original segments
        new_vertices = []
        for seg in range(line.shape[0]-1):
            
            temp = None
            for d in range(D):
                if temp is None:
                    temp    = copy.copy(np.linspace(line[seg,d],line[seg+1,d],num_vertices[seg]+2)[1:-1])
                else:
                    temp    = np.column_stack((
                        temp,
                        copy.copy(np.linspace(line[seg,d],line[seg+1,d],num_vertices[seg]+2)[1:-1]) ))
            new_vertices += [copy.copy(temp)]
        
        # Create the seed for the new line
        new_line    = copy.copy(line[0,:].reshape((1,D)) )
        
        # Assemble the new line
        for seg in range(line.shape[0]-1):
            
            # Add the new segments, then the next original vertex
            new_line    = np.row_stack((
                new_line,
                new_vertices[seg],
                line[seg+1,:]))
        
        return new_line
    
    def shrink_polygon(self,polygon, offset = 1):
    
        """
        This function shrinks a user-provided polygon.
        
        Parameters:
            
            polygon         - Required  : a 2-D array of polygon vertices
            offset          - Required  : a scalar defining the distance by which we wish to shrink the polygon (default = 1)
        """
    
        import numpy as np
        import copy
        import math
        
        def angle(x1, y1, x2, y2):
            numer = (x1*x2 + y1*y2)
            denom = np.sqrt((x1**2 + y1**2) * (x2**2 + y2**2))
            print(numer)
            print(denom)
            print( math.acos(numer/denom) )
            return math.acos(numer/denom) 
        
        def cross_sign(x1, y1, x2, y2):
            return x1*y2 > x2*y1
        
        # If the polygon is closed, un-close it
        closed = False
        if np.linalg.norm(polygon[0,:]-polygon[-1,:]) < 1E-10:
            polygon = polygon[:-1,:]
            closed  = True
        
        # Make sure polygon is counter-clockwise
        if self.are_vertices_clockwise(np.row_stack((polygon,polygon[0,:]))):
            polygon     = np.flipud(polygon)
            
        polygon_shrinked = copy.copy(polygon)
        
        for idx in range(polygon.shape[0]):
            
            if idx == polygon.shape[0]-1:
                vtx_before  = idx-1
                vtx_center  = idx
                vtx_after   = 0
            else:
                vtx_before  = idx-1
                vtx_center  = idx
                vtx_after   = idx+1
                
            side_before = polygon[vtx_center,:] - polygon[vtx_before,:]
            side_after  = polygon[vtx_after,:] - polygon[vtx_center,:]
            
            side_before /= np.linalg.norm(side_before)
            side_after  /= np.linalg.norm(side_after)
            
            nvec_before = np.asarray([-side_before[1],  side_before[0]])
            nvec_after  = np.asarray([-side_after[1],   side_after[0]])
            
            vtx1_before = polygon[vtx_before,:] + nvec_before*offset
            vtx2_before = polygon[vtx_center,:] + nvec_before*offset
            
            vtx1_after  = polygon[vtx_center,:] + nvec_after*offset
            vtx2_after  = polygon[vtx_after,:]  + nvec_after*offset
            
            p       = vtx1_before
            r       = (vtx2_before-vtx1_before)
            
            q       = vtx1_after
            s       = (vtx2_after-vtx1_after)
            
            if np.cross(r,s) == 0:
                
                # Lines are collinear
                polygon_shrinked[idx,:] = vtx2_before
                
            else:
            
                # Lines are not collinear
                t       = np.cross(q - p,s)/(np.cross(r,s))
                
                # This is the intersection point
                polygon_shrinked[idx,:] = p + t*r
                
        if closed:
            polygon_shrinked = np.row_stack((
                polygon_shrinked,
                polygon_shrinked[0,:]))
            
        return polygon_shrinked
    
    def snap_to_domain(self):
        
        """
        This function takes the user-specified polygon and snaps any vertices 
        outside the model domain onto the domain's edge.
        """
        
        import numpy as np
        
        # Calculate the distances of all edge vertices from the center
        dist = np.abs(self.polygon-self.model.domain_center)
        
        # Any vertex with a distance larger than the domain radius lies outside
        indices = np.where(dist > self.model.domain_radius)[0]
        
        # Go through all outside vertices
        for idx in indices:
            
             #Center the vertex to the model
            temp = self.polygon[idx]-self.model.domain_center
            
            # And collapse its norm to unity
            temp /= dist[idx]
            
            # Then scale it to the boundary
            temp *= self.model.domain_radius
            
            # And translate it back to global variables
            temp += self.model.domain_center
            
            # Then save it to the polygon variables
            self.polygon[idx]       = temp
        
        # It is possible that some vertices have folded onto themselves
        repeat = 0
        while repeat != 2:
            
            # Calculate the interior angles
            self.angles     = np.zeros(self.polygon.shape[0])
            for idx in range(self.polygon.shape[0]):
        
                # Set the previous, current, and next vertex of the polygon
                if idx == self.polygon.shape[0]-1:
                    seg_minus   = idx-1
                    seg_center  = idx
                    seg_plus    = 0
                else:
                    seg_minus   = idx-1
                    seg_center  = idx
                    seg_plus    = idx+1
        
                # This is the routine to calculate the interior angle of a vertex
                newtemp     = np.angle(self.polygon[seg_plus]-self.polygon[seg_center]) - \
                              np.angle(self.polygon[seg_center]-self.polygon[seg_minus])
                
                newtemp     -= np.pi
                
                # Restrict it to the range between -pi and +pi
                while newtemp < -np.pi:  newtemp += 2*np.pi
                while newtemp > +np.pi:  newtemp -= 2*np.pi
                
                # Save the angles to the list
                self.angles[idx] = newtemp
                
            # All vertices whose interior angles are less than five degrees are
            # considered degenerate and are removed
            indices = np.ones(self.polygon.shape[0],dtype=bool)
            indices[np.where(np.abs(self.angles) < np.radians(5))[0]]   = False
            
            # Remove the degenerate vertices
            self.polygon = self.polygon[indices]
            
            if np.sum(indices) == len(indices):
                repeat += 1
        
        # Update the dependent variables
        self.polygon_XY     = np.column_stack((
            np.real(self.polygon)[:,np.newaxis],
            np.imag(self.polygon)[:,np.newaxis] ))
        self.segments       = self.polygon.shape[0]
        
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,facecolor='xkcd:silver',edgecolor='xkcd:dark grey',zorder=1,
             alpha=0.5,linewidth=2,**kwargs):
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        plt.fill(np.real(self.polygon),np.imag(self.polygon),edgecolor=edgecolor,
                 facecolor=facecolor,alpha=alpha,zorder=zorder,linewidth=linewidth,**kwargs)
        
#%%

class ElementAreaSink:
    
    def __init__(self, model, polygon, segments = None, strength = 1,
                 variables = [], priors=[], snap_distance = 1E-10,
                 snap = False, influence = None):

        """
        This implements an area sink.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        polygon         - [array]   : either a real N-by-2 matrix or complex vector of length N specifying the vertices of a polygon tracing the element's path
        segments        - [scalar]  : this element has a subdivision function; if a finer resolution than the number of segments in 'line' is desired, specify a larger number here; the function will then subdivide 'line' and 'line_ht' so as to create segments of as equal length as possible
        strength        - [scalar]  : injection or extraction rate of this element in [length unit]^3/[length unit]^2/[time unit]

        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['line_ht']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        """

        import numpy as np
        import copy
        import matplotlib.path
        import math
        
        # Append this element to the specified model
        self.model          = model
        model.elementlist.append(self)
        
        # This element adds water, so it also requires an influence range
        if influence is None:
            self.influence = self.model.domain_radius*2
        else:
            self.influence = influence
        
        # Complexify the polygon, if it isn't already complex
        polygon             = self.complexify(polygon)
        
        # Prepare the polygon variable
        self.polygon        = polygon
        
        # Is the polygon closed? If not, close it temporarily
        self.snap_distance  = snap_distance
        if np.abs(self.polygon[0]-self.polygon[-1]) > self.snap_distance:
            self.polygon    = np.asarray(list(self.polygon)+[self.polygon[0]])
        
        # Also create an array with real coordinates
        self.polygon_XY     = np.column_stack((
            np.real(copy.copy(self.polygon))[:,np.newaxis],
            np.imag(copy.copy(self.polygon))[:,np.newaxis] ))

        # Is the polygon counter-clockwise? If not, correct it
        if self.are_vertices_clockwise(self.polygon_XY):
            self.polygon    = np.flip(self.polygon)
            self.polygon_XY = np.flipud(self.polygon_XY)
            
        # Do we wish to subdivide the polygon?
        # First, check if the user specified a desired segment count
        if segments is None:
            self.segments   = self.polygon.shape[0]-1
        else:
            self.segments   = segments
        
        if self.segments < self.polygon.shape[0]-1:
            raise Exception('Prescribed number of line segments '+str(self.segments)+" mustn't be smaller than the number of vertices "+str(polygon.shape[0]-1)+'.')
        
        # Subdivide the polygon, if desired
        if self.segments > self.polygon.shape[0]-1:
            self.polygon_XY = self.subdivide_line(self.polygon_XY,self.segments)
            self.polygon    = self.polygon_XY[:,0] + 1j*self.polygon_XY[:,1]
            
        # This is a hack: We shrink the polygon by a small amount. This should ensure 
        # that no issues arise from evaluating points directly on the boundary; 
        # there might be other ways to solve this issue alternatively
        self.polygon_XY = self.shrink_polygon(
            polygon = self.polygon_XY,
            offset  = 1E-10)
        self.polygon    = self.polygon_XY[:,0] + 1j*self.polygon_XY[:,1]
            
        # Un-close the polygon again
        self.polygon_XY     = self.polygon_XY[:-1,:]
        self.polygon        = self.polygon[:-1]
            
        # If vertex snapping is enabled, snap all outside vertices onto the domain edge
        if snap:
            self.snap_to_domain()
            
        # =====================================================================
        # Now some area-sink-specific work
        # =====================================================================
        
        # Get the angles of all segments to the x axis
        # required for the local coordinates, Strack 1989, 37.19
        self.alpha = np.zeros(self.segments)
        for seg in range(self.segments):
            if seg == self.segments-1:
                nextseg = 0
            else:
                nextseg = seg+1
            
            # Get the side vector, then normalize it    
            temp            = self.polygon[nextseg]-self.polygon[seg]
            temp            /= np.abs(temp)
            
            self.alpha[seg] = math.asin(np.imag(temp))
            
        
        # Get the central point of the polygon
        self.zc     = np.mean(self.polygon)
            
        # Calculate the area of the polygon with the shoelace formula:
        self.A      = self.get_polygon_area()
            
        # Calculate the coefficients c0, c1, c2 for all segments
        self.L      = np.zeros(self.segments)
        for seg in range(self.segments):
            
            if seg == self.segments-1:
                nextseg = 0
            else:
                nextseg = seg+1
            
            # Save the length of the segment
            self.L[seg]     = np.abs(self.polygon[nextseg]-self.polygon[seg])
        
        # Get strength parameters for each vertex
        self.strength       = strength
        
        # Extract target variables
        self.variables      = variables
        self.priors         = priors
        
        # # Prepare the matrix block containing the effect of this element onto 
        # # itself for future use in solving the linear system. The matrix requires
        # # subtraction of the A_star variable from its diagonal entries for completion
        # self.block          = self.matrix_contribution()
        
        # Check if the prior matches the number of parameters
        if len(self.priors) != len(self.variables):
            raise Exception('Number of priors must match number of unknown variables. Number of priors: '+str(self.priors)+' / Number of unknown variables: '+str(len(self.variables)))
        
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.priors     += [self.priors[idx]]
                self.model.variables  += [var]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']   
                    
    def update(self):
        
        import numpy as np
        
    def get_polygon_area(self):
        import numpy as np
        return 0.5*np.abs(np.dot(np.real(self.polygon),np.roll(np.imag(self.polygon),1))-np.dot(np.imag(self.polygon),np.roll(np.real(self.polygon),1)))
        
    def evaluate_gradient(self,z,detailed = False,derivatives = 'all',override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # # 'detailed' returns the results as a matrix instead of a summed vector
        # if detailed:
        #     grad    = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        # else:
        #     grad    = np.zeros(z.shape, dtype = np.complex)
        
        grad    = np.zeros(z.shape, dtype = np.complex)
        
        # Assemble vectors of Z_plus_1 and Z_minus_1
        Z_plus_1    = []
        Z_minus_1   = []
        log_ratio   = []
        for seg in range(self.segments):
            
            if seg == self.segments-1:
                nextseg = 0
            else:
                nextseg = seg+1
            
            # Get the subtraction and addition
            Z_minus_1.append(2*(copy.copy(z)-self.polygon[nextseg]) /(self.polygon[nextseg]-self.polygon[seg]))
            Z_plus_1.append(2*(copy.copy(z)-self.polygon[seg]) /(self.polygon[nextseg]-self.polygon[seg]))
            log_ratio.append(np.log(Z_minus_1[seg]/Z_plus_1[seg]))
            
        for seg in range(self.segments):
            
            if seg == self.segments-1:
                nextseg = 0
            else:
                nextseg = seg+1
            
                
            dZdz    = 2/(self.polygon[nextseg]-self.polygon[seg])
            
            dHdZ = log_ratio[seg] + (self.polygon[nextseg]-self.polygon[seg])/(z - self.polygon[0])
            
            # Calculate local variables
            Z   = (2*copy.copy(z) - self.polygon[seg] - self.polygon[nextseg])/(self.polygon[nextseg]-self.polygon[seg])
            
            # Get the H function
            # indices      = np.where(np.abs(Z_plus_1[seg]) < 1E-10)[0]
            H   = Z_plus_1[seg]*log_ratio[seg] + 2
            # H[indices]  = 2
            for seg2 in np.arange(seg+1,self.segments,1):
                H   += 2*log_ratio[seg2]
                
            grad    += self.L[seg]**2*(H + (Z-np.conj(Z))*dHdZ)*dZdz
            
        # print(grad)
        
        # Add the pre-factor
        grad    *= self.strength/(32*np.pi*1j)
        
        # Add the second term
        # This equation (8.599) should be divided be 2, not 4, I believe
        # It is derived from equation 8.598, where the factor us 2
        grad    -= self.strength*self.A/(2*np.pi)/(z - self.polygon[0])
        
        
        # grad    += self.strength*self.A/(4*np.pi)/(z - np.mean(self.polygon))
        
        # print(grad)
        
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
        
        return grad
            
    def evaluate(self,z,detailed = False,override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # This follows the approach outlined in Strack 2017, 8.26
        z       = np.asarray(list(z)+[np.mean(self.polygon)+self.influence])
        res     = np.zeros(z.shape,dtype=np.complex)
        
        
        
        # Assemble vectors of Z_plus_1 and Z_minus_1
        Z_plus_1    = []
        log_ratio   = []
        for seg in range(self.segments):
            
            if seg == self.segments-1:
                nextseg = 0
            else:
                nextseg = seg+1
            
            # Get the subtraction and addition
            Z_minus_1 = 2*(z-self.polygon[nextseg]) /(self.polygon[nextseg]-self.polygon[seg])
            Z_plus_1.append(2*(z-self.polygon[seg]) /(self.polygon[nextseg]-self.polygon[seg]))
            log_ratio.append(np.log(Z_minus_1/Z_plus_1[seg]))
            
        
        for seg in range(self.segments):
            
            if seg == self.segments-1:
                nextseg = 0
            else:
                nextseg = seg+1
            
            # Calculate local variables
            # Z   = (z - self.zc[seg])/self.v[seg]
            Z   = (2*z - self.polygon[seg] - self.polygon[nextseg])/(self.polygon[nextseg]-self.polygon[seg])
            
            # Get the H function
            # indices      = np.where(np.abs(Z_plus_1[seg]) < 1E-10)[0]
            H   = Z_plus_1[seg]*log_ratio[seg] + 2
            # H[indices]  = 2
            for seg2 in np.arange(seg+1,self.segments,1):
                H   += 2*log_ratio[seg2]
        
            # Combine that shit
            res += self.L[seg]**2*(Z - np.conj(Z))*H
            
        # Pre-factor
        res     *= self.strength/(32*np.pi*1j)
        
        # And the second term
        res     -= self.strength*self.A/(2*np.pi)*np.log(z-self.polygon[0])
        
        # Correct the offset for the correction factor
        res     -= np.real(res[-1])
        
        # Remove the correction factor
        res     = res[:-1]
            
        return res
        
    def are_vertices_clockwise(self,line):
        
        """
        This function takes an string of 2-D vertices of a polygon, provided as a 
        N x 2 numpy array, and returns a boolean specifying whether the vertices 
        are provided in clock-wise or counter-clock-wise order.
        
        Parameters:
            
            line            - Required  : numpy array of polygon vertices (N by 2)
        """
        
        import numpy as np
        
        signed_area = 0
        for idx in range(line.shape[0]):
        
            x1 = line[idx,0]
            y1 = line[idx,1]
            if idx == line.shape[0]-1:
                x2 = line[0,0]
                y2 = line[0,1]
            else:
                x2 = line[idx+1,0]
                y2 = line[idx+1,1]
        
            signed_area += (x1 * y2 - x2 * y1)
        
        return (np.sign(signed_area) == -1.)
    
    def subdivide_line(self,line,segments):
        
        import numpy as np
        import copy
        
        # If array is one-dimensional, reshape it appropriately
        if len(line.shape) == 1: line = line.reshape((line.shape[0],1))
        
        D   = line.shape[1]
        
        # Calculate the lengths of original segments
        length  = [np.linalg.norm(line[seg,:]-line[seg+1,:]) for seg in range(line.shape[0]-1)]
        
        # Normalize the length of the original segments
        length  /= np.sum(length)
        
        # Calculate the number of new segments we must create, the line already has
        # (#vertices-1) segments. We only require the difference
        new_segments = segments - line.shape[0] + 1
        
        # Calculate where those segments should go
        bins    = np.concatenate(( [0] , np.cumsum(length) ))
        
        # Add Extend the bin length a bit to prevent errors from arithmetic under- or overflow
        bins[0]     -= 1E-10
        bins[-1]    += 1E-10
        
        # Distribute vertices along the segments
        x       = np.linspace(0,1,new_segments)
        num_vertices = []
        for seg in range(line.shape[0]-1):
            if seg == 0:
                num_vertices += [len(np.where(x <= bins[1])[0])]
            else:
                num_vertices += [len(np.where(np.logical_and(
                        x > bins[seg],
                        x <= bins[seg+1]))[0])]
    
        # Subidivide the original segments
        new_vertices = []
        for seg in range(line.shape[0]-1):
            
            temp = None
            for d in range(D):
                if temp is None:
                    temp    = copy.copy(np.linspace(line[seg,d],line[seg+1,d],num_vertices[seg]+2)[1:-1])
                else:
                    temp    = np.column_stack((
                        temp,
                        copy.copy(np.linspace(line[seg,d],line[seg+1,d],num_vertices[seg]+2)[1:-1]) ))
            new_vertices += [copy.copy(temp)]
        
        # Create the seed for the new line
        new_line    = copy.copy(line[0,:].reshape((1,D)) )
        
        # Assemble the new line
        for seg in range(line.shape[0]-1):
            
            # Add the new segments, then the next original vertex
            new_line    = np.row_stack((
                new_line,
                new_vertices[seg],
                line[seg+1,:]))
        
        return new_line
    
    def snap_to_domain(self):
        
        """
        This function takes the user-specified polygon and snaps any vertices 
        outside the model domain onto the domain's edge.
        """
        
        import numpy as np
        
        # Calculate the distances of all edge vertices from the center
        dist = np.abs(self.polygon-self.model.domain_center)
        
        # Any vertex with a distance larger than the domain radius lies outside
        indices = np.where(dist > self.model.domain_radius)[0]
        
        # Go through all outside vertices
        for idx in indices:
            
             #Center the vertex to the model
            temp = self.polygon[idx]-self.model.domain_center
            
            # And collapse its norm to unity
            temp /= dist[idx]
            
            # Then scale it to the boundary
            temp *= self.model.domain_radius
            
            # And translate it back to global variables
            temp += self.model.domain_center
            
            # Then save it to the polygon variables
            self.polygon[idx]       = temp
        
        # It is possible that some vertices have folded onto themselves
        repeat = 0
        while repeat != 2:
            
            # Calculate the interior angles
            self.angles     = np.zeros(self.polygon.shape[0])
            for idx in range(self.polygon.shape[0]):
        
                # Set the previous, current, and next vertex of the polygon
                if idx == self.polygon.shape[0]-1:
                    seg_minus   = idx-1
                    seg_center  = idx
                    seg_plus    = 0
                else:
                    seg_minus   = idx-1
                    seg_center  = idx
                    seg_plus    = idx+1
        
                # This is the routine to calculate the interior angle of a vertex
                newtemp     = np.angle(self.polygon[seg_plus]-self.polygon[seg_center]) - \
                              np.angle(self.polygon[seg_center]-self.polygon[seg_minus])
                
                newtemp     -= np.pi
                
                # Restrict it to the range between -pi and +pi
                while newtemp < -np.pi:  newtemp += 2*np.pi
                while newtemp > +np.pi:  newtemp -= 2*np.pi
                
                # Save the angles to the list
                self.angles[idx] = newtemp
                
            # All vertices whose interior angles are less than five degrees are
            # considered degenerate and are removed
            indices = np.ones(self.polygon.shape[0],dtype=bool)
            indices[np.where(np.abs(self.angles) < np.radians(5))[0]]   = False
            
            # Remove the degenerate vertices
            self.polygon = self.polygon[indices]
            
            if np.sum(indices) == len(indices):
                repeat += 1
        
        # Update the dependent variables
        self.polygon_XY     = np.column_stack((
            np.real(self.polygon)[:,np.newaxis],
            np.imag(self.polygon)[:,np.newaxis] ))
        self.segments       = self.polygon.shape[0]

    def shrink_polygon(self,polygon, offset = 1):
    
        """
        This function shrinks a user-provided polygon.
        
        Parameters:
            
            polygon         - Required  : a 2-D array of polygon vertices
            offset          - Required  : a scalar defining the distance by which we wish to shrink the polygon (default = 1)
        """
    
        import numpy as np
        import copy
        import math
        
        def angle(x1, y1, x2, y2):
            numer = (x1*x2 + y1*y2)
            denom = np.sqrt((x1**2 + y1**2) * (x2**2 + y2**2))
            print(numer)
            print(denom)
            print( math.acos(numer/denom) )
            return math.acos(numer/denom) 
        
        def cross_sign(x1, y1, x2, y2):
            return x1*y2 > x2*y1
        
        # If the polygon is closed, un-close it
        closed = False
        if np.linalg.norm(polygon[0,:]-polygon[-1,:]) < 1E-10:
            polygon = polygon[:-1,:]
            closed  = True
        
        # Make sure polygon is counter-clockwise
        if self.are_vertices_clockwise(np.row_stack((polygon,polygon[0,:]))):
            polygon     = np.flipud(polygon)
            
        polygon_shrinked = copy.copy(polygon)
        
        for idx in range(polygon.shape[0]):
            
            if idx == polygon.shape[0]-1:
                vtx_before  = idx-1
                vtx_center  = idx
                vtx_after   = 0
            else:
                vtx_before  = idx-1
                vtx_center  = idx
                vtx_after   = idx+1
                
            side_before = polygon[vtx_center,:] - polygon[vtx_before,:]
            side_after  = polygon[vtx_after,:] - polygon[vtx_center,:]
            
            side_before /= np.linalg.norm(side_before)
            side_after  /= np.linalg.norm(side_after)
            
            nvec_before = np.asarray([-side_before[1],  side_before[0]])
            nvec_after  = np.asarray([-side_after[1],   side_after[0]])
            
            vtx1_before = polygon[vtx_before,:] + nvec_before*offset
            vtx2_before = polygon[vtx_center,:] + nvec_before*offset
            
            vtx1_after  = polygon[vtx_center,:] + nvec_after*offset
            vtx2_after  = polygon[vtx_after,:]  + nvec_after*offset
            
            p       = vtx1_before
            r       = (vtx2_before-vtx1_before)
            
            q       = vtx1_after
            s       = (vtx2_after-vtx1_after)
            
            if np.cross(r,s) == 0:
                
                # Lines are collinear
                polygon_shrinked[idx,:] = vtx2_before
                
            else:
            
                # Lines are not collinear
                t       = np.cross(q - p,s)/(np.cross(r,s))
                
                # This is the intersection point
                polygon_shrinked[idx,:] = p + t*r
                
        if closed:
            polygon_shrinked = np.row_stack((
                polygon_shrinked,
                polygon_shrinked[0,:]))
            
        return polygon_shrinked

    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z

    def plot(self,facecolor_extract='xkcd:orangish red',
             edgecolor_extract='xkcd:crimson',facecolor_inject='xkcd:cerulean',
             edgecolor_inject='xkcd:cobalt',zorder=12,alpha=0.5,linewidth=5,**kwargs):
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        if self.strength < 0:
            col_face    = facecolor_extract
            col_edge    = edgecolor_extract
        else:
            col_face    = facecolor_inject
            col_edge    = edgecolor_inject
        
        plt.fill(np.real(self.polygon),np.imag(self.polygon),edgecolor=col_edge,
                 facecolor=col_face,alpha=alpha,zorder=zorder,linewidth=linewidth,
                 **kwargs)

#%%

class ElementLineSink:
    
    def __init__(self, model, line, segments = None, influence = None, 
                 strength = 1, variables = [], priors=[]):
        
        """
        This implements a line sink.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        line            - [array]   : either a real N-by-2 matrix or complex vector of length N specifying the vertices of a line string tracing the element's path
        segments        - [scalar]  : this element has a subdivision function; if a finer resolution than the number of segments in 'line' is desired, specify a larger number here; the function will then subdivide 'line' and 'line_ht' so as to create segments of as equal length as possible
        influence       - [scalar]  : radius of zero influence of each line segment; set to twice the model domain_radius if unspecified
        strength        - [scalar]  : this specifies the injection or extraction rate of the element

        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['line_ht']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        """

        import numpy as np
        from scipy.interpolate import interp1d
        import copy
        
        self.model          = model
        model.elementlist.append(self)
        
        self.variables      = variables
        self.priors         = priors
        
        # ---------------------------------------------------------------------
        # Subdivide the provided no flow boundary into #segments pieces
        
        self.line_raw       = copy.copy(line)
        
        if segments is None:
            
            self.segments   = line.shape[0]-1
            
        else:
            self.segments   = segments
            
        if self.segments < self.line_raw.shape[0]-1:
            
            raise Exception('Number of segments '+str(self.segments)+" mustn't be smaller than number of line points "+str(line.shape[0])+'.')
              
        if self.segments > self.line_raw.shape[0]:
            
            # Subdivide the line
            self.line       = self.subdivide_line(line,self.segments)
            self.line_c     = copy.copy(self.line[:,0] + 1j*self.line[:,1])
        else:
            
            self.line       = self.line_raw.copy()
            self.line_c     = self.line[:,0] + 1j*self.line[:,1]
            
        # Also get the normal vector components to each segment
        self.line_nvec  = self.line[:,1] - 1j*self.line[:,0]
        self.line_nvec  = self.line_nvec/np.abs(self.line_nvec)

        # ---------------------------------------------------------------------        
        
        
        
        
        self.strength       = np.ones(self.segments)*strength
        
        if influence is None:
            self.influence  = self.model.domain_radius*2
        else:
            self.influence  = influence
        
        
        self.Zi             = []
        self.offset_outside = []
        self.L              = []
        self.zc             = []
        self.segment_nvec   = []
        self.head_target    = []
        
        for seg in range(self.segments):
            
            self.L              += [np.abs(self.line_c[seg+1] - self.line_c[seg])]
            
            influence_pt        = (self.line_c[seg+1]-self.line_c[seg])*self.influence/self.L[seg] + self.line_c[seg]
            Z                   = (2*influence_pt-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg])
            self.Zi             += [copy.copy(Z)]
            
            self.zc             += [(self.line_c[seg]+self.line_c[seg+1])/2]
            
            # Calculate the normal vector to this segment
            self.segment_nvec   += [(self.line_c[seg]-self.line_c[seg+1])]
            self.segment_nvec[-1]= [np.imag(self.segment_nvec[-1])-1j*np.real(self.segment_nvec[-1])]
            
            part1               = np.nan_to_num((Z+1)*np.log(Z+1))
            part2               = np.nan_to_num((Z-1)*np.log(Z-1))
            self.offset_outside += [self.L[seg] / (4*np.pi) * (part1 - part2)]
            
        # Convert list of segment centers to array
        self.zc             = np.asarray(self.zc)
            
        
        # Check if the prior matches the number of parameters
        if len(self.priors) != len(self.variables):
            raise Exception('Number of priors must match number of unknown variables. Number of priors: '+str(self.priors)+' / Number of unknown variables: '+str(len(self.variables)))
        
        # Go through all elements
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.priors     += [self.priors[idx]]
                self.model.variables  += [var]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']   
                    
    def update(self):
        
        import numpy as np
        
        # self.zc             = self.xc + 1j*self.yc
        # self.L              = np.abs(self.z2 - self.z1)
        
        # influence_pt        = (self.z2-self.z1)*self.influence/self.L + self.z1
        # Z                   = (2*influence_pt-(self.z1+self.z2))/(self.z2-self.z1)
        # part1               = np.nan_to_num((Z+1)*np.log(Z+1))
        # part2               = np.nan_to_num((Z-1)*np.log(Z-1))
        # self.offset_outside = self.L / (4*np.pi) * (part1 - part2)

    def evaluate_gradient(self,z,detailed = False, derivatives = 'all', override_parameters = False):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # 'detailed' returns the results as a matrix instead of a summed vector
        if detailed:
            grad    = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            grad    = np.zeros(z.shape, dtype = np.complex)
        
        for seg in range(self.segments):
        
            Z       = (2*z-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg])
            
            # Now get the gradient d omega(z)/dZ
            
            if not override_parameters:
                temp    = self.strength[seg]*self.L[seg]/4/np.pi*(np.log(Z+1) - np.log(Z-1))
            else:
                temp    = self.L[seg]/4/np.pi*(np.log(Z+1) - np.log(Z-1))
            
            # To get d omega(z)/dz we can use the product rule
            #   d omega(z)/dz = d omega(z)/dZ * dZ/dz
            # hence:
            temp    = temp*2/(self.line_c[seg+1]-self.line_c[seg])
            
            if detailed:
                grad[seg,:]     = copy.copy(temp)
            else:
                grad    += temp
        
        
        
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
        
        
        return grad

    def evaluate(self,z,detailed = False, override_parameters = False):
        
        import copy
        import numpy as np
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        if detailed:
            res     = np.zeros((self.segments,z.shape[0]), dtype = np.complex)
        else:
            res     = np.zeros(z.shape, dtype = np.complex)
        
        for seg in range(self.segments):
            
            # Convert to local coordinates
            Z       = (2*z-(self.line_c[seg]+self.line_c[seg+1]))/(self.line_c[seg+1]-self.line_c[seg])
            
            # Evaluate the complex potential offset by a distance in the 
            if not override_parameters:
                temp    = self.strength[seg]*self.L[seg]/4/np.pi * (\
                          (Z+1)*np.log(Z+1) - \
                          (Z-1)*np.log(Z-1) - \
                          (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                          (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
            else:
                temp    = self.L[seg]/4/np.pi * (\
                          (Z+1)*np.log(Z+1) - \
                          (Z-1)*np.log(Z-1) - \
                          (2/self.L[seg]*self.influence+2)*np.log(2/self.L[seg]*self.influence+2) + \
                          (2/self.L[seg]*self.influence)*np.log(2/self.L[seg]*self.influence))
            
            
            # If evaluated directly at the endpoints, the result would be NaN
            # They should be zero, see Bakker 2009
            temp    = np.nan_to_num(temp)
            
            if detailed:
                res[seg,:]      = copy.copy(temp)
            else:
                res             += temp
            
        return res
    
    def subdivide_line(self,line,segments):
        
        import numpy as np
        import copy
        
        # If array is one-dimensional, reshape it appropriately
        if len(line.shape) == 1: line = line.reshape((line.shape[0],1))
        
        D   = line.shape[1]
        
        # Calculate the lengths of original segments
        length  = [np.linalg.norm(line[seg,:]-line[seg+1,:]) for seg in range(line.shape[0]-1)]
        
        # Normalize the length of the original segments
        length  /= np.sum(length)
        
        # Calculate the number of new segments we must create, the line already has
        # (#vertices-1) segments. We only require the difference
        new_segments = segments - line.shape[0] + 1
        
        # Calculate where those segments should go
        bins    = np.concatenate(( [0] , np.cumsum(length) ))
        
        # Add Extend the bin length a bit to prevent errors from arithmetic under- or overflow
        bins[0]     -= 1E-10
        bins[-1]    += 1E-10
        
        # Distribute vertices along the segments
        x       = np.linspace(0,1,new_segments)
        num_vertices = []
        for seg in range(line.shape[0]-1):
            if seg == 0:
                num_vertices += [len(np.where(x <= bins[1])[0])]
            else:
                num_vertices += [len(np.where(np.logical_and(
                        x > bins[seg],
                        x <= bins[seg+1]))[0])]
    
        # Subidivide the original segments
        new_vertices = []
        for seg in range(line.shape[0]-1):
            
            temp = None
            for d in range(D):
                if temp is None:
                    temp    = copy.copy(line[seg,d] + (line[seg+1,d]-line[seg,d]) * np.linspace(0,1,num_vertices[seg]+2)[1:-1])
                else:
                    temp    = np.column_stack((
                        temp,
                        copy.copy(line[seg,d] + (line[seg+1,d]-line[seg,d]) * np.linspace(0,1,num_vertices[seg]+2)[1:-1])))
                
            new_vertices += [copy.copy(temp)]
        
        # Create the seed for the new line
        new_line    = copy.copy(line[0,:].reshape((1,D)) )
        
        # Assemble the new line
        for seg in range(line.shape[0]-1):
            
            # Add the new segments, then the next original vertex
            new_line    = np.row_stack((
                new_line,
                new_vertices[seg],
                line[seg+1,:]))
        
        return new_line
    
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,color_extract='xkcd:orangish red',color_inject='xkcd:cerulean',
             zorder=12,linewidth=5,**kwargs):
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        if self.strength < 0:
            col     = color_extract
        else:
            col     = color_inject
        
        plt.plot(np.real(self.line_c),np.imag(self.line_c),color=col,
                 zorder=zorder,linewidth=linewidth,**kwargs)
    
#%%
    
class ElementWell:
    
    def __init__(self, model, zc, rw, influence = None, head_change = -1, strength = 1,
                 drawdown_specified = False, variables = [], priors = []):
        
        """
        This implements an injection or extraction well.
        
        Parameters:
            
        model           - [object]  : the model object to which this element is added
        zw              - [vector]  : either a complex scalar or a real vector of length 2 specifying the xy coordinates of the well
        rw              - [scalar]  : a real scalar specifying the screen radius of the well in [length units]
        strength        - [scalar]  : extraction or injection rate at this well in [length units]^3/[time units]
        head_change     - [scalar]  : alternative to strength, induces the prescribed drawdown at the well; only used if drawdown_specified is True
        drawdown_specified - [boolean]  : flag for whether the well's strength is determined through a prescribed head_change; defaults to False
        
        If MCMC is used, we further require:
            
        variable        - [list]    : list of variables which are inferred by MCMC, example: ['line_ht']; leave empty if unused
        priors          - [list]    : list of dictionaries, one for each unknown 'variable'; each dictionary must contain the name of distribution (in scipy.stats) and the relevant parameters as keys
        """
        
        import numpy as np
        
        self.model      = model
        model.elementlist.append(self)
        
        self.variables      = variables
        self.priors         = priors
        
        if influence is None:
            # If no influence radius is specified, set it to twice the model radius
            self.influence  = 2*self.model.domain_radius
        else:
            # Otherwise, set it to the user-defined value
            self.influence  = influence
            
        # The well's strength defines its effect on the flow field; this is
        # overwritten later on to achieve the desired head_change which depends
        # on the aquifer parameters
        self.strength       = strength
        
        # This is the well's position in terms of complex coordinates
        self.zc             = zc
        if not np.isscalar(self.zc):
            self.zc         = self.zc[0] + 1j*self.zc[1]
        
        # The well radius is specified in canonical units
        self.rw             = rw
        
        # Check if drawdown specified
        self.drawdown_specified = drawdown_specified
        
        if self.drawdown_specified:
            
            # Get the head change variable
            self.head_change    = head_change
        
            # Adjust the strength so that the desired drawdown is achieved
            self.set_potential_target()
        
        # Check if the prior matches the number of parameters
        if len(self.priors) != len(self.variables):
            raise Exception('Number of priors must match number of unknown variables. Number of priors: '+str(self.priors)+' / Number of unknown variables: '+str(len(self.variables)))
        
        # Go through all elements
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.priors     += [self.priors[idx]]
                self.model.variables  += [var]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']  
                    
    def update(self):
        
        if self.drawdown_specified:
            self.set_potential_target()
        
    def evaluate_gradient(self,z,derivatives = 'all'):
        
        import numpy as np
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Find the indices of wells      
        dist            = np.abs(z-self.zc)
        idx_inside      = np.where(dist < self.rw)[0]
        idx_outside     = np.where(dist > self.influence)[0]
        idx_valid       = [i for i in np.arange(len(z),dtype=int) if (i not in idx_inside and i not in idx_outside)]
        
        # Correct the coordinates ---------------------------------------------
        # Set the well center to the origin of the complex plane
        zs              = z.copy()-self.zc

        # Pre-allocate an array for the gradient
        grad            = np.zeros(zs.shape,    dtype=np.complex)
        
        # Calculate the gradient
        # grad[idx_valid] = self.strength/(zs[idx_valid]*np.log(self.rw/self.influence))
        grad[idx_valid] = -self.strength/(2*np.pi)/zs[idx_valid]
        
        # If partial derivatives are demanded, calculate them
        if derivatives == 'phi': # phi corresponds to the hydraulic potential
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi': # psi corresponds to the stream function
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all': # all returns the complex derivative
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
        
        return grad

    def evaluate(self,z):
        
        import numpy as np
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Effect at influence range
        temp    = -self.strength/(2*np.pi)*np.log(self.influence)
        
        # Find the indices of wells      
        dist            = np.abs(z-self.zc)
        idx_inside      = np.where(dist < self.rw)[0]
        
        # Correct the coordinates ---------------------------------------------
        # Center the evaluation points on the well
        zs  = z.copy()-self.zc
        
        # Snap points inside the well to the well edge
        zs[idx_inside] = self.rw + 0j    
        
        # Calculate the complex potential
        res = -self.strength/(2*np.pi)*np.log(zs) - temp

        return res
    
    def set_potential_target(self):
        
        """
        We define the drawdown in terms of head, but for the calculations we 
        require it in terms of potential.
        """
        
        import copy
        import numpy as np
        
        # Get the hydraulic conductivity
        for e in self.model.elementlist:
            if isinstance(e, ElementMoebiusBase) or isinstance(e, ElementUniformBase):
                temp_k = e.k
        
        for e in self.model.elementlist:
            if isinstance(e, ElementInhomogeneity):
                if e.are_points_inside_polygon(self.zc):
                    temp_k = e.k
        
        # Create a list of hydraulic potential targets
        self.strength = copy.copy(self.head_change)
        if self.model.aquifer_type == 'confined':
            # Strack 1989, Eq. 8.6
            self.strength = temp_k*self.model.H*self.strength - \
                0.5*temp_k*self.model.H**2
        elif self.model.aquifer_type == 'unconfined':
            # Strack 1989, Eq. 8.7
            self.strength = 0.5*temp_k*self.strength**2
        elif self.model.aquifer_type == 'convertible':
            # Find out which points are confined and which are unconfined
            index_conf      = np.where(self.strength >= self.model.H)[0]
            index_unconf    = np.where(self.strength < self.model.H)[0]
            # Account for the confined points
            # confined:     Strack 1989, Eq. 8.6
            self.strength[index_conf] = \
                temp_k[index_conf]*self.model.H*self.strength[index_conf] - \
                0.5*temp_k[index_conf]*self.model.H**2
            # unconfined:   Strack 1989, Eq. 8.7
            self.strength[index_unconf] = \
                0.5*temp_k[index_unconf]*self.strength[index_unconf]**2
                
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,color_extract='xkcd:orangish red',color_inject='xkcd:cerulean',
             zorder=15,s=100,edgecolor='xkcd:dark grey',linewidth=2,**kwargs):
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        if self.strength < 0:
            col     = color_extract
            marker  = '^'
        else:
            col     = color_inject
            marker  = 'v'
        
        plt.scatter(np.real(self.zc),np.imag(self.zc),s=s,color=col,
                 marker=marker,zorder=zorder,edgecolor=edgecolor,
                 linewidth=linewidth,**kwargs)
    
#%%

class ElementMoebiusOverlay:
    
    def __init__(self,model,r=None,a=None,b=None,c=None,d=None,head_min=0,head_max=1,
                 variables=[],priors=[],angular_limit=1):
        
        """
        Similar to the Möbius base, but an additive overlay. Unfinished.

        """
        
        import numpy as np
        
        # Append the base to the elementlist
        self.model          = model
        model.elementlist.append(self)
        
        # Define an angular limit. This is designed to keep the Möbius control 
        # points from getting arbitrarily close to each other; defined in radians
        self.angular_limit  = angular_limit
        
        # Set Moebius values
        self.r              = r
        self.a              = a
        self.b              = b
        self.c              = c
        self.d              = d
        
        # Set potential scaling variables
        self.head_min       = head_min
        self.head_max       = head_max
        
        # The model requires the base flow in terms of hydraulic potential (phi)
        # The function head_to_potential extracts the following variables:
        #   phi_min         hydraulic potential corresponding to head_min
        #   phi_max         hydraulic potential corresponding to head_max
        self.head_to_potential()
        
        
        # Check input for validity
        self.check_input()
        
        # Define the original control points in the Moebius base disk
        self.z0     = np.asarray(
            [np.complex(np.cos(-0.25*np.pi),np.sin(-0.25*np.pi)),
             np.complex(np.cos(0.25*np.pi),np.sin(0.25*np.pi)),
             np.complex(np.cos(0.75*np.pi),np.sin(0.75*np.pi))])
        
        # If only rotation is specified, get the Moebius coefficients
        if self.r is not None and (self.a is None and self.b is None and \
                                   self.c is None and self.d is None):
            # Find Moebius coefficients
            self.find_moebius_coefficients()
        
        self.variables      = variables
        self.priors         = priors
        
        self.Ke             = 1.854
        
        if len(self.variables) > 0:
            # There are some model variables specified
            for idx,var in enumerate(self.variables):
                self.model.num_params += 1
                exec("self.model.params += [self.%s]" % var)
                self.model.variables  += [var]
                self.model.priors     += [self.priors[idx]]
                if 'name' in list(self.priors[idx].keys()):
                    self.model.param_names  += [self.priors[idx]['name']]   
                else:                       
                    self.model.param_names  += ['unknown']  
        
    def update(self):
        
        # If this model is updated, make sure to repeat any initialization
        # Find Moebius coefficients
        self.find_moebius_coefficients()
            
    def check_input(self):
        
        import numpy as np
        
        # See if either control point rotations or a full set of Moebius
        # coefficients are specified
        if self.r is None and (self.a is None or self.b is None or \
                               self.c is None or self.d is None):
            raise Exception('Either control point rotations r or Moebius coefficients a, b, c, and d must be specified.')
        
        # Check if phi_min is smaller than phi_offset, switch if necessary
        if self.phi_min > self.phi_max:
            raise Exception('Minimum potential phi_min is larger than maximum potential phi_max.')
    
        # Check if the control points fulfill the minimum angular spacing
        r   = np.degrees(self.r)
        if np.abs((r[0]-r[1] + 180) % 360 - 180) < self.angular_limit or \
           np.abs((r[1]-r[2] + 180) % 360 - 180) < self.angular_limit or \
           np.abs((r[2]-r[0] + 180) % 360 - 180) < self.angular_limit:
            raise Exception('Control points '+str(self.r)+' are too close to each other. Define different control points or adjust the angular limit: '+str(self.angular_limit))
    
    def evaluate(self,z):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Coordinates in canonical space are the start values
        z_canonical     = copy.copy(z)
        
        # Scale the canonical disk to unity canonical disk
        z = (z - self.model.domain_center)/self.model.domain_radius
        
        # Map from canonical disk to Möbius base
        z = self.moebius(z,inverse=True)
        
        # Map from Möbius base to unit square
        z = self.disk_to_square(z)
        
        # Rescale the complex potential
        z   = (z+1)/2 * (self.phi_max-self.phi_min) + self.phi_min
            
        return z
    
    def evaluate_gradient(self,z,derivatives = 'all'):
        
        import numpy as np
        import copy
        
        # Complexify the evaluation points, if they aren't already complex
        z       = self.complexify(z)
        
        # Map from the canonical disk to Möbius base
        z_mb = (copy.copy(z) - self.model.domain_center)/self.model.domain_radius
        
        # dz_mb / dz_c
        grad_4  = 1/self.model.domain_radius
        
        # Map from Möbius base to unit disk
        z_ud = self.moebius(copy.copy(z_mb),inverse=True)
        
        # dz_ud / dz_mb
        grad_3  = (self.a*self.d-self.b*self.c)/(self.c*z_mb-self.a)**2
        
        grad_2  = 2/(self.Ke*np.sqrt(z_ud**4+1))
        
        grad_1  = (self.phi_max-self.phi_min)/2
                    
        grad    = grad_1*grad_2*grad_3*grad_4
            
        if derivatives == 'phi':
            dudx    = np.real(grad)
            dudy    = -np.imag(grad)
            grad    = dudx + 1j*dudy
        elif derivatives == 'psi':
            dvdx = -np.imag(np.conjugate(grad))
            dvdy = np.real(np.conjugate(grad))
            grad    = dvdx + 1j*dvdy
        elif derivatives != 'all':
            raise Exception("'derivatives' has to be either 'all' (complex derivative), " + \
                "'phi' (hydraulic potential partial derivatives), or 'psi' " + \
                "(flow line partial derivatives)")
            
        return grad
    
    def complex_integral(self,func,a,b):
        
        """
        This implements the Gauss-Kronrod integration for complex-valued functions.
        We use this to evaluate the Legendre incomplete elliptic integral of the 
        first kind, since it is about ten times as fast as using mpmath's ellipf
        function. Since this integration is a major computational bottleneck of
        this function, we stick with this approach.
        
        The equations below are adapted from: 
        https://stackoverflow.com/questions/5965583/use-scipy-integrate-quad-to-integrate-complex-numbers
        """
        
        import scipy
        from scipy import array
        
        def quad_routine(func, a, b, x_list, w_list):
            c_1 = (b-a)/2.0
            c_2 = (b+a)/2.0
            eval_points = map(lambda x: c_1*x+c_2, x_list)
            func_evals = list(map(func, eval_points))    # Python 3: make a list here
            return c_1 * sum(array(func_evals) * array(w_list))
        
        def quad_gauss_7(func, a, b):
            x_gauss = [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, 0.405845151377397, 0.741531185599394, 0.949107912342759]
            w_gauss = array([0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277,0.129484966168870])
            return quad_routine(func,a,b,x_gauss, w_gauss)
        
        def quad_kronrod_15(func, a, b):
            x_kr = [-0.991455371120813,-0.949107912342759, -0.864864423359769, -0.741531185599394, -0.586087235467691,-0.405845151377397, -0.207784955007898, 0.0, 0.207784955007898,0.405845151377397, 0.586087235467691, 0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813]
            w_kr = [0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728, 0.204432940075298, 0.190350578064785, 0.169004726639267, 0.140653259715525,  0.104790010322250, 0.063092092629979, 0.022935322010529]
            return quad_routine(func,a,b,x_kr, w_kr)
        
        class Memorize:                     # Python 3: no need to inherit from object
            def __init__(self, func):
                self.func = func
                self.eval_points = {}
            def __call__(self, *args):
                if args not in self.eval_points:
                    self.eval_points[args] = self.func(*args)
                return self.eval_points[args]
        
        def quad(func,a,b):
            ''' Output is the 15 point estimate; and the estimated error '''
            func = Memorize(func) #  Memorize function to skip repeated function calls.
            g7 = quad_gauss_7(func,a,b)
            k15 = quad_kronrod_15(func,a,b)
            # I don't have much faith in this error estimate taken from wikipedia
            # without incorporating how it should scale with changing limits
            return [k15, (200*scipy.absolute(g7-k15))**1.5]
        
        return quad(func,a,b)
    
    def angle_to_unit_circle(self):
        
        import numpy as np
        
        # Angle must be provided in radians, counter-clockwise from 3 o'clock
        return np.cos(self.r)+1j*np.sin(self.r)
    
    def find_moebius_coefficients(self):
        
        import numpy as np
        
        # Find the images of the z0 control points
        w0  = self.angle_to_unit_circle()
        
        # Then calculate the four parameters for the corresponding Möbius map
        self.a = np.linalg.det(np.asarray(
            [[self.z0[0]*w0[0],     w0[0],          1],
             [self.z0[1]*w0[1],     w0[1],          1],
             [self.z0[2]*w0[2],     w0[2],          1]]))
        
        self.b = np.linalg.det(np.asarray(
            [[self.z0[0]*w0[0],     self.z0[0],     w0[0]],
             [self.z0[1]*w0[1],     self.z0[1],     w0[1]],
             [self.z0[2]*w0[2],     self.z0[2],     w0[2]]]))
        
        self.c = np.linalg.det(np.asarray(
            [[self.z0[0],           w0[0],          1],
             [self.z0[1],           w0[1],          1],
             [self.z0[2],           w0[2],          1]]))
        
        self.d = np.linalg.det(np.asarray(
            [[self.z0[0]*w0[0],     self.z0[0],     1],
             [self.z0[1]*w0[1],     self.z0[1],     1],
             [self.z0[2]*w0[2],     self.z0[2],     1]]))
        
        return
    
    def moebius(self,z,inverse=False):
        
        if not inverse:
            z = (self.a*z+self.b)/(self.c*z+self.d)
        else:
            z = (-self.d*z+self.b)/(self.c*z-self.a)
        
        return z
    
    def square_to_disk(self,z,k='default'):
        
        import numpy as np
        from mpmath import mpc,mpmathify,ellipfun
        
        if k == 'default': k = 1/mpmathify(np.sqrt(2))
        
        Ke = 1.854
        cn = ellipfun('cn')
    
        if type(z) is complex: 
            z = np.asarray([z])
        zf  = np.ndarray.flatten(z)
        w   = np.zeros(zf.shape)*1j
        
        pre_factor  = mpc(1,-1)/mpmathify(np.sqrt(2))
        mid_factor  = Ke*(mpc(1,1)/2)
    
        for idx,entry in enumerate(zf): # Go through all complex numbers
            
            # Calculate result
            temp = pre_factor*cn(
                u = mid_factor*entry-Ke,
                k = k)
            
            # Then place it into the array
            w[idx]  = np.complex(temp.real,temp.imag)
            
        # Now reshape the array back to its original shape
        z = w.reshape(z.shape).copy()
        
        return z
    
    def disk_to_square(self,z,k='default'):
        
        import numpy as np
    
        Ke = 1.854

        if type(z) is complex: 
            z = np.asarray([z])
        zf  = np.ndarray.flatten(z)
        w   = np.zeros(zf.shape)*1j

        # Using the Gauss-Kronrod integration is about 10 times faster than 
        # using the mpmath.ellipf function
        if k == 'default': k = 1/np.sqrt(2)
        m = k**2
        pre_factor  = (1-1j)/(-Ke)
        mid_factor  = (1+1j)/np.sqrt(2)
        
        temp = [pre_factor*self.complex_integral(
            func    = lambda t: (1-m*np.sin(t)**2)**(-0.5), 
            a       = 0, 
            b       = i)[0] + 1 - 1j for i in np.arccos(zf*mid_factor)]
    
        w = np.asarray(temp)
            
        # Now reshape the array back to its original shape
        z = w.reshape(z.shape).copy()
        
        return z
    
    def are_points_clockwise(self):
        
        import numpy as np
        
        verts = np.zeros((3,2))
        
        verts[0,:] = np.asarray([np.cos(self.r[0]),np.sin(self.r[0])])
        verts[1,:] = np.asarray([np.cos(self.r[1]),np.sin(self.r[1])])
        verts[2,:] = np.asarray([np.cos(self.r[2]),np.sin(self.r[2])])
        
        signed_area = 0
        for vtx in range(verts.shape[0]):
            x1 = verts[vtx,0]
            y1 = verts[vtx,1]
            if vtx == verts.shape[0]-1: # Last vertex
                x2 = verts[0,0]
                y2 = verts[0,1]
            else:
                x2 = verts[vtx+1,0]
                y2 = verts[vtx+1,1]
            signed_area += (x1 * y2 - x2 * y1)/2
            
        return (signed_area < 0)
       
    def head_to_potential(self):
        
        # Extract the hydraulic conductivity from the base element
        k = self.model.elementlist[0].k
        
        for idx,h in enumerate([self.head_min,self.head_max]):
        
            if self.model.aquifer_type == 'confined' or (self.model.aquifer_type == 'convertible' and h >= self.model.H):
                # Strack 1989, Eq. 8.6
                pot = k*self.model.H*h - 0.5*k*self.model.H**2
                    
            elif self.model.aquifer_type == 'unconfined' or (self.model.aquifer_type == 'convertible' and h < self.model.H):
                # Strack 1989, Eq. 8.7
                pot = 0.5*k*h**2
                
            if idx == 0:
                self.phi_min    = pot
            elif idx == 1:
                self.phi_max    = pot
                
    def complexify(self,z):
        
        """
        This function takes the provided line or polygon and converts it into
        a complex-valued vector, if it isn't already provided as one.'
        """
        
        import numpy as np
        
        if not np.iscomplex(z).any():
            if len(z.shape) != 2 or z.shape[1] != 2:
                raise Exception('Shape format not understood. Provide shape vertices either as a complex vector, or as a N-by-2 real numpy array.')
            else:
                z   = z[:,0] + 1j*z[:,1]
                
        return z
    
    def plot(self,label_offset = 1.1,fontsize=12,fontcolor='xkcd:dark grey',
             pointcolor='xkcd:dark grey',pointsize=50,horizontalalignment='center',
             verticalalignment='center',color_low = 'xkcd:cerulean',
             color_high = 'xkcd:orangish red',**kwargs):
        
        """
        This function plots the Möbius control/reference points on the unit disk.
        """
        
        import numpy as np
        import matplotlib.pyplot as plt
        import math
        
        # Get the coordinates of the control points
        z_A     = (1-1j)/np.abs(1-1j)
        z_A     = self.moebius(z_A,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_A     = np.asarray([np.real(z_A),np.imag(z_A)])
        
        z_B     = (1+1j)/np.abs(1+1j)
        z_B     = self.moebius(z_B,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_B     = np.asarray([np.real(z_B),np.imag(z_B)])
        
        z_C     = (-1+1j)/np.abs(-1+1j)
        z_C     = self.moebius(z_C,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_C     = np.asarray([np.real(z_C),np.imag(z_C)])
        
        z_D     = (-1-1j)/np.abs(-1-1j)
        z_D     = self.moebius(z_D,inverse=False)*self.model.domain_radius + self.model.domain_center
        z_D     = np.asarray([np.real(z_D),np.imag(z_D)])
        
        a_low   = np.linspace(math.atan2(z_C[1],z_C[0]),math.atan2(z_D[1],z_D[0]),360)
        if abs(a_low[0]-a_low[-1]) > np.pi:
            a_low = np.concatenate((
                np.linspace(np.min(a_low),-np.pi,360),
                np.linspace(np.pi,np.max(a_low),360) ))
        
        a_high  = np.linspace(math.atan2(z_A[1],z_A[0]),math.atan2(z_B[1],z_B[0]),360)
        if abs(a_high[0]-a_high[-1]) > np.pi:
            a_high = np.concatenate((
                np.linspace(np.min(a_high),-np.pi,360),
                np.linspace(np.pi,np.max(a_high),360) ))
            
        plt.plot(np.cos(a_low)*self.model.domain_radius + self.model.domain_center,
                 np.sin(a_low)*self.model.domain_radius + self.model.domain_center,
                 color = color_low,linewidth=2)
        plt.plot(np.cos(a_high)*self.model.domain_radius + self.model.domain_center,
                 np.sin(a_high)*self.model.domain_radius + self.model.domain_center,
                 color = color_high,linewidth=2)
        
        plt.scatter(z_A[0],z_A[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        plt.scatter(z_B[0],z_B[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        plt.scatter(z_C[0],z_C[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        plt.scatter(z_D[0],z_D[1],s=pointsize,color=pointcolor,zorder=11,**kwargs)
        
        plt.text(z_A[0]*label_offset,z_A[1]*label_offset,'A',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        plt.text(z_B[0]*label_offset,z_B[1]*label_offset,'B',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        plt.text(z_C[0]*label_offset,z_C[1]*label_offset,'C',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)
        plt.text(z_D[0]*label_offset,z_D[1]*label_offset,'D',fontsize=fontsize,
                 horizontalalignment=horizontalalignment,verticalalignment=verticalalignment,
                 color=fontcolor,**kwargs)

#%%
    
def equidistant_points_in_circle(rings = 3, radius = 1, offset = 0+0j):
    
    """
    This function creates equidistant points on a number of specified rings
    inside a unit disk.
    
    Parameters:
        
    rings           - [scalar]  : number of rings on which equidistant points are placed; the more rings, the more points
    radius          - [scalar]  : radius by which the unit disk is scaled
    """
    
    import numpy as np
    import math
    
    # If the offset is complex, convert it to a real vector of length 2
    if np.iscomplex(offset).any():
        if not np.isscalar(offset):
            raise Exception('Shape format not understood. Provide the offset either as a complex scalar, or as a real numpy array of shape (2,).')
        else:
            offset = np.asarray([np.real(offset),np.imag(offset)])
    if np.isscalar(offset):
        offset  = np.zeros(2)
    
    # Pre-allocate lists for the X and Y coordinates
    x = []
    y = []
    
    # Then go through each ring
    for k in range(rings):
        
        if k > 0:
            pts = round(np.pi/math.asin(1/(2*k)))
        else:
            pts = 1
        
        theta = np.linspace(0, 2*np.pi, pts)
    
        rad = k/(rings-1)
    
        x += list(np.sin(theta)*rad)
        y += list(np.cos(theta)*rad)
    
    # Combine both lists to a common array
    XY = np.column_stack((
        np.asarray(x),
        np.asarray(y)))
    
    # And scale it, if desired
    XY  *= radius
    
    # Apply the offset
    XY[:,0]     += offset[0]
    XY[:,1]     += offset[1]
    
    return XY
    

