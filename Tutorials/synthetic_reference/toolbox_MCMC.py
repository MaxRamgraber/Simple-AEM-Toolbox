class MCMC:
    
    def __init__(self, model, observations = [], chain_length = 100, likelihood_dictionary = None,
                 adapt_proposal = True, adapt_delay = 100,adapt_frequency = 10,
                 adapt_factor = 0.1, acceptance_target = None):

        import copy
        
        """
        This class initializes an instance of Markov Chain Monte Carlo routine.
        
        Parameters:
            model           - [class]   : The model instance for which this element is defined
            observations    - [list]    : a list containing dictionaries which specify observation locations, with keys 'name' [string], 'location' [complex], and 'head' [scalar]
            chain_length    - [scalar]  : number specifying the target length of the Markov Chain; larger values improve the statistical reliability
            likelihood_dictionary - [dict]  : a dictionary specifying the likelihood function; should contain keys 'distribution' and any necessary keywords for the distribution
            adapt_proposal  - [boolean] : a boolean specifying whether the MCMC proposal distribution is adjusted or not
            adapt_frequency - [scalar]  : number specifying how often (every * MCMC steps) the proposal distribution is adapted
            adapt_factor    - [scalar]  : starting value for the adaptation factor f; this value scales the covariance of the chain's unique values to yield the proposal covariance
            acceptance_target   - [scalar]  : if proposal adaptation is used, this defines the desired target acceptance ratio of MCMC proposals. Should be around 0.238
        """
        
        # Pass the AEM model
        self.model                  = model
        
        # Define the information required to calculate likelihoods
        self.observations           = observations
        self.likelihood_dictionary  = likelihood_dictionary
        
        # Define the length of the Markov Chain
        self.chain_length           = chain_length
        
        # Define the variables required for adaptive proposals
        self.adapt_proposal         = adapt_proposal
        self.adapt_delay            = adapt_delay
        self.adapt_frequency        = adapt_frequency
        self.adapt_factor           = adapt_factor
        self.acceptance_target      = acceptance_target
        
        # Assemble the proposal distribution dictionary
        self.proposal_dictionary    = []
        for entry in self.model.proposals:
            self.proposal_dictionary.append(copy.deepcopy(entry))
        for e in self.model.elementlist:
            for entry in e.proposals:
                self.proposal_dictionary.append(copy.deepcopy(entry))
        
        self.check_input()
        

    def check_input(self):
        
        if self.proposal_dictionary is None:
            raise Exception('No proposal_distribution specified.')
        
        # if self.model.r is None and 'r' in self.model.variables:
        #     raise Exception("To consider the control point rotation as an unknown variable, the model's 'r' variable must be specified." )

        # if len(self.model.variables) == 0:
        #     raise Exception("MCMC requires the designation of at least one unknown variable to infer.")
            
        elif len(self.model.variables) > 0:
            # There are some model variables specified
            listvar     = list(self.model.__dict__.keys())
            for var in self.model.variables:
                if var not in listvar:
                    raise Exception("Unknown variable '"+var+"' not among the base model parameters: "+str(listvar))

        for e in self.model.elementlist:
            if len(e.variables) > 0:
                # There are some model variables specified
                listvar     = list(e.__dict__.keys())
                for var in e.variables:
                    if var not in listvar:
                        raise Exception("Unknown variable '"+var+"' not among the element's parameters: "+str(listvar))

        if self.likelihood_dictionary is None:
            raise Exception('No likelihood specify. Define a likelihood_dictionary.')
            
    def start(self):
        
        import numpy as np
        import copy
        from toolbox_AEM import ElementMoebiusBase,ElementMoebiusOverlay
        
        # Prepare the start of the chain --------------------------------------
        
        # Get the original logprior
        logprior                    = self.model.logprior(
            params                  = self.model.params,
            priors                  = self.model.priors)
        
        # Get the original loglikelihood
        loglikeli, residuals        = self.model.loglikelihood(
            observations            = self.observations,
            likelihood_dictionary   = self.likelihood_dictionary,
            predictions             = None)
        
        # Combine both to obtain the posterior
        self.logposterior_original  = copy.copy(np.sum(logprior)) + copy.copy(np.sum(loglikeli))
        
        # Initiate the chain
        self.chain                  = [copy.copy(self.model.params)]
        self.chain_logposterior     = [copy.copy(self.logposterior_original)]
        self.chain_residuals        = [copy.copy(residuals)]
        self.chain_duplicates       = [1]
        
        # Initialize a cycle counter used to reduce the adaptation factor for
        # the proposal adaptation
        cycle_counter               = 0
        
        # Find the base element
        MoebiusBase_index           = None
        for idx,e in enumerate(self.model.elementlist):
            if isinstance(e, ElementMoebiusBase) and 'r' in e.variables:
                MoebiusBase_index   = idx
                
        # Find any Möbius overlay element
        MoebiusOverlay_index        = None
        for idx,e in enumerate(self.model.elementlist):
            if isinstance(e, ElementMoebiusOverlay) and 'r' in e.variables:
                if MoebiusOverlay_index is None:
                    MoebiusOverlay_index    = [idx]
                else:
                    MoebiusOverlay_index    += [idx]
        
        # Start the chain -----------------------------------------------------
        
        for iteration in range(self.chain_length - 1):
            
            # Adapt the proposal density --------------------------------------
            
            if self.adapt_proposal and iteration >= self.adapt_delay and \
            (iteration == self.adapt_delay or iteration%self.adapt_frequency == 0):
        
                # Prepare the proposal distribution adaptation if there is an
                # acceptance rate target
                if self.acceptance_target is not None:
                    
                    cycle_counter   += 1
                    
                    # Calculate the acceptance rate
                    innovation = []
                    for i in range(len(self.chain_duplicates)):
                        innovation += [self.chain_duplicates[-i-1]]
                        if np.sum(innovation) >= self.adapt_frequency:
                            break
                    
                    acceptance_rate = len(innovation)/self.adapt_frequency
                    
                    # Prevent the adaptation factor from collapsing to zero
                    if acceptance_rate == 0: acceptance_rate = 0.01
                    
                    # Calculate the ratio between acceptance rate and target
                    update_factor   = acceptance_rate/self.acceptance_target
                    
                    # Subtract 1 to center it, then scale it by 0.1, the add 1 to shift it back
                    update_factor   -= 1
                    update_factor   = update_factor*(1.1**(-cycle_counter))
                    update_factor   += 1
                           
                    # Adjust the adapt_factor
                    self.adapt_factor   = self.adapt_factor * update_factor
                
                # Adjust the proposal distribution
                self.adapt_proposal_density()  
            
            # Mutate the current sample ---------------------------------------

            # Propose a new sample (corresponding to a new model instance)
            self.proposed_model = copy.deepcopy(self.model)
            
            # Add a new proposal
            self.proposed_model.params = copy.copy(self.proposal())
            
            # Update the proposed model
            self.proposed_model.update()
            
            if MoebiusBase_index is not None:
                if self.model.elementlist[MoebiusBase_index].are_points_clockwise():
                    raise Exception
                    
            if MoebiusOverlay_index is not None:
                for idx in MoebiusOverlay_index:
                    if self.model.elementlist[idx].are_points_clockwise():
                        raise Exception
            
            # Get the logposterior of the proposal ----------------------------
            
            # Get the proposed logprior
            logprior                    = self.proposed_model.logprior(
                params                  = self.proposed_model.params,
                priors                  = self.proposed_model.priors)
            
            # Check if the sample violates any limits
            if logprior is None:
                reject_immediately  = True
            else:
                reject_immediately  = False
                
            # The Möbius base has special rejection requirements
            if MoebiusBase_index is not None:
                
                # Check if the rotation of the control points is counter-clockwise
                if self.proposed_model.elementlist[MoebiusBase_index].are_points_clockwise():
                    reject_immediately  = True
                    
                # Check if the control points fulfill the minimum angular spacing
                r               = np.degrees(self.proposed_model.elementlist[MoebiusBase_index].r)
                angular_limit   = np.degrees(self.proposed_model.elementlist[MoebiusBase_index].angular_limit)
                if np.abs((r[0]-r[1] + 180) % 360 - 180) < angular_limit or \
                   np.abs((r[1]-r[2] + 180) % 360 - 180) < angular_limit or \
                   np.abs((r[2]-r[0] + 180) % 360 - 180) < angular_limit:
                    reject_immediately  = True
                    
            # The Möbius base has special rejection requirements
            if MoebiusOverlay_index is not None:
                
                for idx in MoebiusOverlay_index:
                
                    # Check if the rotation of the control points is counter-clockwise
                    if self.proposed_model.elementlist[idx].are_points_clockwise():
                        reject_immediately  = True
                        
                    # Check if the control points fulfill the minimum angular spacing
                    r               = np.degrees(self.proposed_model.elementlist[idx].r)
                    angular_limit   = np.degrees(self.proposed_model.elementlist[idx].angular_limit)
                    if np.abs((r[0]-r[1] + 180) % 360 - 180) < angular_limit or \
                       np.abs((r[1]-r[2] + 180) % 360 - 180) < angular_limit or \
                       np.abs((r[2]-r[0] + 180) % 360 - 180) < angular_limit:
                        reject_immediately  = True
                
                
            # Combine both to obtain the posterior if the proposal was valid
            if not reject_immediately:
                
                # Get the original loglikelihood, only if not immediately rejected
                loglikeli, residuals = self.proposed_model.loglikelihood(
                    observations            = self.observations,
                    likelihood_dictionary   = self.likelihood_dictionary,
                    predictions             = None)
                
                self.logposterior_proposal  = copy.copy(np.sum(logprior)) + copy.copy(np.sum(loglikeli))
                
            # Decide whether to accept or reject the proposal -----------------
            
            # Accept or reject the sample based on
            if not reject_immediately:
                
                print('\r'+str(iteration)+': '+str(np.exp(self.logposterior_proposal-self.logposterior_original)),end='\r')
                
                # Throw the dice
                if np.random.uniform(0,1,1)[0] <= np.exp(self.logposterior_proposal-self.logposterior_original):
                    
                    # The proposal has been accepted
                    accept  = True
                    
                else:
                    
                    # The proposal has been rejected
                    accept  = False
            
            else: 
                
                accept = False # The proposal has been rejected immediately

            
            # Append results to chain -----------------------------------------
            
            if accept:
                
                # Append the proposal to the chain
                self.chain                  += [copy.deepcopy(self.proposed_model.params)]
                self.chain_logposterior     += [copy.deepcopy(self.logposterior_proposal)]
                self.chain_residuals        += [copy.deepcopy(residuals)]
                self.chain_duplicates       += [1]
            
                # Overwrite the original entries
                self.model                  = copy.deepcopy(self.proposed_model)
                self.logposterior_original  = copy.deepcopy(self.logposterior_proposal)
                
            else:
                
                # Append the result to the chain
                self.chain_duplicates[-1]   += 1
                
        print('MCMC run finished.')
                
    def adapt_proposal_density(self):
        
        import numpy as np
        import scipy.stats
        
        # Go through all parameters
        for idx,prop in enumerate(self.proposal_dictionary):
            
            # Go through all entries in the chain
            array   = []
            for entry in self.chain:
                array.append(entry[idx])
            array   = np.asarray(array)
            
            if len(array.shape) > 1: # This variable is multi-dimensional
                
                if 'circular' in list(prop.keys()):
                    
                    if prop['circular'] is True:
                        
                        std     = scipy.stats.circstd(
                            samples     = array,
                            low         = - np.pi,
                            high        = np.pi,
                            axis        = 0)
                        
                    else:
                        
                        std     = np.std(
                            a           = array,
                            axis        = 0)
                        
                else:
                    
                    std     = np.std(
                        a           = array,
                        axis        = 0)
                    
                # Calculate the covariance matrix
                cov     = np.identity(array.shape[1])
                np.fill_diagonal(cov,(std*self.adapt_factor)**2)
                
                # Write it into the dictionary
                self.proposal_dictionary[idx]['cov'] = cov.copy()
                
            else:
                
                std     = np.std(
                    a           = array,
                    axis        = 0)
                    
                # Write it into the dictionary
                self.proposal_dictionary[idx]['scale'] = std*self.adapt_factor
            
            
            
    def proposal(self):
        
        import scipy.stats
        import copy
        import math
        import numpy as np
        
        params_proposal = copy.deepcopy(self.model.params)
        
#        print(params_proposal)
        
        # Go through all parameters
        for idx,val in enumerate(params_proposal):
            
            # The proposal distribution is a univariate Gaussian
            if  self.proposal_dictionary[idx]['distribution'] == 'norm' or \
                self.proposal_dictionary[idx]['distribution'] == 'normal':
                      
                # The mean isn't strictly necessary to define, make sure it's there
                if 'loc' not in list(self.proposal_dictionary[idx].keys()):
                    self.proposal_dictionary[idx]['loc'] = 0
                   
                addition = scipy.stats.norm.rvs(
                    loc     = self.proposal_dictionary[idx]['loc'],
                    scale   = self.proposal_dictionary[idx]['scale'],
                    size    = 1)[0]
            
            # The proposal distribution is a multivariate Gaussian
            elif    self.proposal_dictionary[idx]['distribution'] == 'multivariate_normal' or \
                    self.proposal_dictionary[idx]['distribution'] == 'multivariate normal':
            
                # The mean isn't strictly necessary to define, make sure it's there
                if 'mean' not in list(self.proposal_dictionary[idx].keys()):
                    self.proposal_dictionary[idx]['mean'] = np.zeros(self.proposal_dictionary[idx]['cov'].shape[0])
                   
                addition = scipy.stats.multivariate_normal.rvs(
                    mean    = self.proposal_dictionary[idx]['mean'],
                    cov     = self.proposal_dictionary[idx]['cov'],
                    size    = 1)
                
            else: raise Exception('Proposal distribution type not understood: ' + \
                  str(self.proposal_dictionary[idx]['distribution']) + ' \ ' + \
                  "implemented types: 'norm', 'multivariate_normal'")
                
            # Check if the user specified any converter
            converter_used = False
            if 'converter' in list(self.proposal_dictionary[idx].keys()) and 'deconverter' in list(self.proposal_dictionary[idx].keys()):
                
                # Activate the logarithmic boolean and save the base
                converter_used  = True
                
                # Extract the converter and the deconverter
                converter       = self.proposal_dictionary[idx]['converter']
                deconverter     = self.proposal_dictionary[idx]['deconverter']
                
                # And convert the variable
                params_proposal[idx] = converter(params_proposal[idx])
                
            # Now add the proposal
            params_proposal[idx] += copy.copy(addition)
            
            # If a converter was used, also use the deconverter
            if converter_used:
                
                # Deconvert the variable
                params_proposal[idx] = deconverter(params_proposal[idx])
                        
            # # Add the mutation to the original parameters
            # # To do so, first check whether the variable in question is logarithmic
            # if 'logarithmic base' in list(self.proposal_dictionary[idx].keys()):
                
            #     # Convert the variable to logarithmic space
            #     params_proposal[idx] = math.log(
            #         params_proposal[idx], 
            #         self.proposal_dictionary[idx]['logarithmic base'])
                
            #     # Add the proposal
            #     params_proposal[idx] += copy.copy(addition)
                
            #     # Then convert it back to canonical space
            #     params_proposal[idx] = self.proposal_dictionary[idx]['logarithmic base']**params_proposal[idx]
                
            # else:
                
            #     # If the variable is not logarithmic, simply add the noise in canonical space
            #     params_proposal[idx] += copy.copy(addition)
            
            # If this parameter is ciruclar (has an entry 'circular' : True),
            # constrain it to the range between -pi and +pi
            if 'circular' in list(self.proposal_dictionary[idx].keys()):
                
                # In case someone wants to define 'circular' as False
                if self.proposal_dictionary[idx]['circular'] is True:
                    
                    if len(params_proposal[idx]) > 1:
                        
                        for ix,entry in enumerate(params_proposal[idx]):
                            
                            if entry > np.pi:
                                params_proposal[idx][ix] -= 2*np.pi
                            
                            if entry < -np.pi:
                                params_proposal[idx][ix] += 2*np.pi
                    else:
                    
                        if params_proposal[idx] > np.pi:
                            params_proposal[idx] -= 2*np.pi
                            
                        if params_proposal[idx] < -np.pi:
                            params_proposal[idx] += 2*np.pi
                            
#        print(params_proposal)
        
        return params_proposal
    
    def evaluate_chain(self, z, cutoff = 0, subsampling = 1, 
        subsampling_type = 'random', mode = 'head', derivatives = 'all'):
        
        import copy
        import numpy as np
        from random import sample
        
        # Convert the evaluation points to complex variables, if necessary
        z           = self.complexify(z)
        
        # Create copies of the basic stuff
        z_base      = copy.deepcopy(z)
        model_base  = copy.deepcopy(self.model)
        
        chain_base              = copy.deepcopy(self.chain)
        chain_logposterior_base = copy.deepcopy(self.chain_logposterior)
        chain_residuals_base    = copy.deepcopy(self.chain_residuals)
        chain_duplicates_base   = copy.deepcopy(self.chain_duplicates)
        
        chain_base              = chain_base[cutoff:]
        chain_logposterior_base = chain_logposterior_base[cutoff:]
        chain_residuals_base    = chain_residuals_base[cutoff:]
        chain_duplicates_base   = chain_duplicates_base[cutoff:]
        
        # Apply subsampling, if any was specified. Two subsampling_types are
        # accepted:
        #
        #   'random'        : In this case, 'subsampling' specifies a percentage
        #                     between 0 and 1 which is kept, so 0.1 corresponds
        #                     to 10% of the samples. The chosen samples are
        #                     selected randomly.
        #
        #   'incremental'   : Works similarly to 'random', but instead of random
        #                     sampling, we sample in regular increments.
        
        if subsampling_type == 'random':
            
            indices     = []
            for i in [[idx]*num for idx,num in enumerate(chain_duplicates_base)]:
                indices += i
            num_samples = int(np.ceil(len(indices)*subsampling))
            indices     = sample(indices, k = num_samples)
            indices     = np.sort(indices)
            
            indices, chain_duplicates_base = np.unique(
                indices,
                return_counts = True)
            
        elif subsampling_type == 'incremental':
            
            indices     = []
            for i in [[idx]*num for idx,num in enumerate(chain_duplicates_base)]:
                indices += i
            num_samples = int(np.ceil(len(indices)*subsampling))
            indices     = np.asarray(indices)[np.linspace(0,len(chain_base),num_samples).astype(int)]
            indices     = np.sort(indices)
            
            indices, chain_duplicates_base = np.unique(
                indices,
                return_counts = True)
            
        else:
            
            raise Exception("subsampling_type not understood. It should be either 'random' or 'incremental'.")
            
        chain_base              = [entry for idx,entry in enumerate(chain_base) if idx in indices]
        chain_logposterior_base = [entry for idx,entry in enumerate(chain_logposterior_base) if idx in indices]
        chain_residuals_base    = [entry for idx,entry in enumerate(chain_residuals_base) if idx in indices]
        
        results                 = []
        
        for idx,val in enumerate(chain_base):
            
            iteration   = idx+1
            total       = len(chain_base)
            percent = ("{0:." + str(1) + "f}").format(100 * (iteration / float(total)))
            filledLength = int(50 * iteration // total)
            bar = '█' * filledLength + '-' * (50 - filledLength)
            print('\r%s |%s| %s%% %s' % ('Progress:', bar, percent, ''), end = '\r')
            
            # Propose a new sample (corresponding to a new model instance)
            model = copy.deepcopy(model_base)
            
            # Add a new proposal
            model.params = copy.deepcopy(val)
            
            # Update the proposed model
            model.update()
            
            # Extract the complex potential at the query points
            results     += [copy.copy(model.evaluate(copy.copy(z_base),mode=mode,derivatives=derivatives))]
            
        self.evaluation_dictionary = {
            'chain'             : copy.deepcopy(chain_base),
            'chain_logposterior': copy.deepcopy(chain_logposterior_base),
            'chain_residuals'   : copy.deepcopy(chain_residuals_base),
            'chain_duplicates'  : copy.deepcopy(chain_duplicates_base),
            'indices'           : copy.deepcopy(indices),
            'results'           : copy.deepcopy(results),
            'z'                 : copy.deepcopy(z_base)}
        
        return self.evaluation_dictionary
    
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
    
#%%

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    
    import numpy as np
    
    average = np.average(values, weights=weights,axis = 0)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights,axis = 0)
    return (average, np.sqrt(variance))
