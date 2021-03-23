class MCMC:
    
    def __init__(self, model, observations = [], chain_length = 100, 
                 ensemble = None, likelihood_dictionary = None, 
                 walkers = 3, epsilon = 1E-6, cutoff = 50, flat = True,
                 clean_chains = True):

        import copy
        import numpy as np
        
        """
        This class initializes an instance of Differential Evolution Markov Chain Monte Carlo routine.
        
        Parameters:
            model               	- [class]   : The model instance for which this element is defined
            observations            - [list]    : a list containing dictionaries which specify observation locations, with keys 'name' [string], 'location' [complex], and 'head' [scalar]
            chain_length            - [scalar]  : number specifying the target length of the Markov Chain; larger values improve the statistical reliability
            ensemble                - [array]   : a N-by-D matrix, where N is the number of walkers and D is the number of coefficients; calculated based on prior information, if not provided
            likelihood_dictionary   - [dict]    : a dictionary specifying the likelihood function; should contain keys 'distribution' and any necessary keywords for the distribution
            walkers                 - [integer] : number of walkers used for the Differential Evolution Markov Chain
            epsilon                 - [scalar]  : standard deviation for small perturbation term
            cutoff                  - [integer] : number of chain steps to be discarded as burn-in
            flag                    - [boolean] : flag for whether results should be returned as the individual walkers' chains (False) or concatenated into a single term (True)
            clean_chains            - [boolean] : flag for whether the final walkers' chains should automatically remove all samples prior to convergence
        """
        
        # Pass the AEM model
        self.model                  = model
        
        # Define the information required to calculate likelihoods
        self.observations           = observations
        self.likelihood_dictionary  = likelihood_dictionary
        
        # Define the length of the Markov Chain
        self.chain_length           = chain_length
        
        # Define the MCMC cutoff
        self.cutoff                 = cutoff
        
        # Define the variables required for adaptive proposals
        self.walkers                = walkers
        self.epsilon                = epsilon
        
        self.clean_chains           = clean_chains
        self.flat                   = flat
                
        # Calculate the number of flattened parameters
        self.C                      = np.sum([1 if np.isscalar(x) else len(x) for x in self.model.params])
        
        # Draw an initial ensemble if none was specified
        if ensemble is None:
            self.draw_initial_ensemble()
        else:
            self.ensemble   = ensemble
        
        self.check_input()
        
        # Create a flag for whether the chain was already run or not
        self.chain_initialized      = False
        

    def check_input(self):
        
        """
        This function checks for possible input errors in the specification of
        the MCMC object.
        """
        
        # if self.model.r is None and 'r' in self.model.variables:
        #     raise Exception("To consider the control point rotation as an unknown variable, the model's 'r' variable must be specified." )

        if len(self.model.variables) == 0:
            raise Exception("MCMC requires the designation of at least one unknown variable to infer.")
            
        elif len(self.model.variables) > 0:
            # There are some model variables specified
            listvar     = list(self.model.__dict__.keys())
            for var in self.model.variables:
                if var not in listvar:
                    raise Exception("Specified unknown variable '"+var+"' not among the base model parameters: "+str(listvar))

        for e in self.model.elementlist:
            if len(e.variables) > 0:
                # There are some model variables specified
                listvar     = list(e.__dict__.keys())
                for var in e.variables:
                    if var not in listvar:
                        raise Exception("Unknown variable '"+var+"' not among the element's parameters: "+str(listvar))

        if self.likelihood_dictionary is None:
            raise Exception('No likelihood specify. Define a likelihood_dictionary.')
            
        if self.walkers <= self.C:
            print('WARNING: Number of walkers ('+str(self.walkers)+') is smaller or equal to number of parameter space dimensions ('+str(self.C)+'). This risks confining exploration of parameter space to a subspace.')
            
        if self.walkers < 3:
            raise Exception('The differential requires at least three walkers. Currently specified are '+str(self.walkers)+' walkers.')
            
        if self.ensemble.shape != (self.walkers,self.C):
            raise Exception('Shape of ensemble '+str(self.ensemble.shape)+' does not match expected '+str(self.walkers)+'-'+str(self.C)+' shape.')
            
    def draw_initial_ensemble(self):
        
        """
        If no initial ensemble was provided, the algorithm initializes its own
        starting ensemble based on the prior information provided.
        """
        
        import numpy as np
        import copy
        import scipy.stats
        
        # Initialize an ensemble from the prior
        self.ensemble   = np.zeros((self.walkers,self.C))
        
        # Draw samples for each walker
        for walker in range(self.walkers):
            
            # Repeat the sampling until a valid set of samples is found.
            repeat      = True
            while repeat:
            
                params_draw     = []
                
                # Go through each prior
                for idx,entry in enumerate(self.model.priors):
                    
                    # Create a copy of the dictionary, as we must delete some entries
                    entry_local = copy.deepcopy(entry)
                    
                    # Are there any limits?
                    if 'limits' in list(entry_local.keys()):
                        # If yes, save them and delete the dictionary key
                        limits  = entry_local['limits']
                        del entry_local['limits']
                    else:
                        # If not, remember that
                        limits  = None
                        
                    # Are there any converters and deconverters?
                    if 'converter' in list(entry_local.keys()):
                        # If yes, save them and delete the dictionary key
                        converter   = entry_local['converter']
                        deconverter = entry_local['deconverter']
                        del entry_local['converter'], entry_local['deconverter']
                    else:
                        # If not, remember that
                        converter   = None
                        deconverter = None
                        
                    # Remove any name
                    if 'name' in list(entry_local.keys()):
                        del entry_local['name']
                        
                    # Store the distribution, then delete the key
                    if 'distribution' in list(entry_local.keys()):
                        distribution    = entry_local['distribution']
                        del entry_local['distribution']
                        
                    # Create the function string for drawing the prior ensemble,
                    # then execute it
                    exec('params_draw.append(scipy.stats.'+distribution+'.rvs(**entry_local))')
                    
                    # If the variable we drew was the Möbius control point rotation,
                    # make sure that it is between -pi and +pi
                    if self.model.variables[idx] == 'r':
                        params_draw[-1][np.where(params_draw[-1] > np.pi)]  -= 2*np.pi
                        params_draw[-1][np.where(params_draw[-1] < -np.pi)] += 2*np.pi
                    
                for idx,entry in enumerate(self.model.priors):
                        
                    # Are there any converters and deconverters?
                    if 'converter' in list(entry.keys()):
                        params_draw[idx]    = entry['deconverter'](params_draw[idx])
                        
                        
                # Check if the model predicts non-zero water tables
                self.model.params = copy.deepcopy(params_draw)
                
                # Update the proposed model
                self.model.update()
                
                # Get the original logprior; this conveniently checks whether 
                # the limits were honored, returning "None" if they were not
                logprior                    = self.model.logprior(
                    params                  = self.model.params,
                    priors                  = self.model.priors)
                
                # Get the original loglikelihood
                loglikeli, residuals        = self.model.loglikelihood(
                    observations            = self.observations,
                    likelihood_dictionary   = self.likelihood_dictionary,
                    predictions             = None)
                
                # Logprior can return "None" if any of the limits are violated,
                # Loglikely can return None if a cell falls dry or artefacts 
                # occur
                if logprior is not None and loglikeli is not None:
                    # If no failure is detected, accept this parameter set as 
                    # valid
                    repeat  = False

            # Flatten the parameters
            params_draw = self.unpack_parameters(np.asarray(params_draw))
            
            # Store the results
            self.ensemble[walker,:] = copy.copy(params_draw)
           
    def start(self):
        
        """
        This function is the heart piece of this toolbox, setting up and running
        the MCMC chain.        
        """
        
        import numpy as np
        import copy
        import scipy.stats
        from toolbox_AEM import ElementMoebiusBase,ElementMoebiusOverlay
        
        import warnings
        warnings.filterwarnings('ignore')
        
        # Prepare the start of the chain --------------------------------------
        
        # This is the chain of unique parameters
        self.chain                  = [[]]*self.walkers
        
        # This is the chain of all parameters, storing duplicates
        self.chain_full             = [copy.copy(self.ensemble)]
        
        # This is the chain for the unnormalized logposteriors
        self.chain_logposterior     = [[1]]*self.walkers
        
        # This is the chain for the observation mismatches
        self.chain_residuals        = [[]]*self.walkers
        
        # This is the chain for the duplicates of unique parameter sets
        self.chain_duplicates       = [[1]]*self.walkers
        
        # Initialize the chain ------------------------------------------------
        for nw in range(self.walkers):
            
            # Add a new proposal
            self.model.params = copy.deepcopy(self.pack_parameters(self.ensemble[nw,:]))
            
            # Update the proposed model
            self.model.update()
            
            # Get the original logprior
            logprior                    = self.model.logprior(
                params                  = self.model.params,
                priors                  = self.model.priors)
            
            # Get the original loglikelihood
            loglikeli, residuals        = self.model.loglikelihood(
                observations            = self.observations,
                likelihood_dictionary   = self.likelihood_dictionary,
                predictions             = None)
            
            # Check whether the ensemble is valid
            if logprior is None or loglikeli is None:
                raise Exception("Either the logprior ("+str(logprior)+") or the loglikelihood ("+str(loglikeli)+") of walker "+str(nw)+" returned an error flag ('None') during initialization. Make sure that the ensemble only contains valid particles. If in doubt, not specifying an ensemble will cause this toolbox to generate its own valid ensemble based on prior information.")
            
            # Combine both to obtain the posterior
            logposterior_original  = copy.copy(np.sum(logprior)) + copy.copy(np.sum(loglikeli))
            
            # Append results
            self.chain[nw]              = [copy.copy(self.ensemble[nw,:])]
            self.chain_logposterior[nw] = [logposterior_original]
            self.chain_residuals[nw]    = [copy.copy(residuals)]
            
        # Store the latest ensemble
        self.chain_full             = [copy.copy(self.ensemble)]
        
        # Set the chain as initialized
        self.chain_initialized      = True
        
        # Initialize an iteration counter
        iteration_counter   = 0
        
        # Print the header for the MCMC progress bar
        print('--------------------------------------------------------------')
        print('iteration counter / acc. rate (last 100 steps) / acc. (overall) : posterior density ratio')
        
        # Start the chain -----------------------------------------------------
        for iteration in range(self.chain_length - 1 + self.cutoff):
            
            # Increment the iteration counter
            iteration_counter   += 1
            
            # Once we reach the cutoff, reset the chains ----------------------
            if iteration == self.cutoff:
                
                # Reset the chains for each walker
                for nw in range(self.walkers):
                    self.chain[nw]              = [self.chain[nw][-1]]
                    self.chain_logposterior[nw] = [self.chain_logposterior[nw][-1]]
                    self.chain_residuals[nw]    = [self.chain_residuals[nw][-1]]
                    self.chain_duplicates[nw]   = [1]
                
                # Reset the full storage chain
                self.chain_full             = [self.chain_full[-1]]
                
                # Reset the iteration counter
                iteration_counter           = 1
                
            # Pre-allocate space for a new ensemble entry
            new_ensemble    = np.zeros(self.ensemble.shape)
                
            # Make a jump for each walker -------------------------------------
            for nw in range(self.walkers):
                # For each walker, select a pair of other walkers to obtain the
                # direction of the jump vector. First, find all other walkers.
                other_walkers   = [x for x in np.arange(self.walkers) if x != nw]
                
                # Choose two other walkers without replacement
                walker_choices  = np.random.choice(other_walkers,size=2,replace=False)
                
                # # Select two samples in those walkers' chains
                # sample_choices  = np.random.choice(np.arange(iteration_counter),size=2,replace=True)
                
                # Obtain the direction vector between those two samples
                # jump_direction  = \
                #     self.chain_full[sample_choices[0]][walker_choices[0],:] - \
                #     self.chain_full[sample_choices[1]][walker_choices[1],:]
                jump_direction  = \
                    self.chain_full[-1][walker_choices[0],:] - \
                    self.chain_full[-1][walker_choices[1],:]
                
                # Draw a stepsize from a beta distribution, then make the proposal
                proposal        = \
                    copy.copy(self.chain_full[-1][nw,:]) + \
                    jump_direction*scipy.stats.beta.rvs(a=1,b=4)
                
                # Add some random noise to prevent degeneracy
                proposal        += scipy.stats.norm.rvs(
                    loc         = 0,
                    scale       = self.epsilon,
                    size        = proposal.shape)
                
                # Create a copy of the model, add the proposal ----------------
                
                # Propose a new sample (corresponding to a new model instance)
                self.proposed_model = copy.deepcopy(self.model)
                
                # Add a new proposal
                self.proposed_model.params = copy.deepcopy(self.pack_parameters(proposal))
                
                # Update the proposed model
                self.proposed_model.update()
                
                # Get the logposterior of the proposal ------------------------
                
                # Get the proposed logprior
                logprior                    = self.proposed_model.logprior(
                    params                  = self.proposed_model.params,
                    priors                  = self.proposed_model.priors)
                
                # Check if the sample violates any limits
                if logprior is None:
                    reject_immediately  = True
                else:
                    reject_immediately  = False
                    
                # Combine both to obtain the posterior if the proposal was valid
                if not reject_immediately:
                    
                    # Get the original loglikelihood, only if not immediately rejected
                    loglikeli, residuals = self.proposed_model.loglikelihood(
                        observations            = self.observations,
                        likelihood_dictionary   = self.likelihood_dictionary,
                        predictions             = None)
                    
                    if loglikeli is None:
                        reject_immediately  = True
                    
                # Get the original logposterior density
                logposterior_original   = self.chain_logposterior[nw][-1]
                
                # If both prior and likelihood were evaluated successfully, 
                # calculate the logposterior of the proposal
                if not reject_immediately:
            
                    # Get the proposal logposterior
                    logposterior_proposal  = \
                        copy.copy(np.sum(logprior)) + copy.copy(np.sum(loglikeli))
                
                # Decide whether to accept or reject the proposal -----------------
                
                # If the sample is valid, decide whether to reject it or not
                if not reject_immediately:
                    
                    # Prepare a vector to calculate the overall acceptance rate
                    # of this walker's chain so far
                    acc_full        = []
                    for val in self.chain_duplicates[nw]:
                        acc_full    += [1]
                        acc_full    += [0]*(val-1)
                        
                    # If something goes horribly, horribly wrong, and we somehow
                    # did not catch the error before, terminate here.
                    # If this ever happens, smear goat's blood on the side of
                    # your computer, ritually purify this code under a full moon,
                    # then lie in a corner and cry.
                    if np.isnan(np.exp(logposterior_proposal-logposterior_original)):
                        print(logposterior_proposal)
                        print(logposterior_original)
                        print(loglikeli)
                        print(logprior)
                        raise Exception
                    
                    # Print current chain information
                    #   - the current iteration
                    #   - the acceptance rate of this walker for the last 100 steps
                    #   - the acceptance rate of this walker overall
                    #   - the ratio between the proposal and original posterior density
                    string  = '\r'+str(iteration_counter)+' '+\
                        str(np.sum(acc_full[-100:])/np.min([len(acc_full),100]))[:4]+' '+\
                        str(len(self.chain_duplicates[nw])/np.sum(self.chain_duplicates[nw]))[:4]+\
                        ': '+str(np.exp(logposterior_proposal-logposterior_original))
                    print(string,end='\r')

                    # Throw the dice, did we accept?
                    if np.random.uniform(0,1,1)[0] <= np.exp(logposterior_proposal-logposterior_original):
                        
                        # The proposal has been accepted
                        accept  = True
                        
                    else:
                        
                        # The proposal has been rejected
                        accept  = False
                
                else: 
                    
                    accept = False # The proposal has been rejected immediately
    
                
                # Append results to chain -----------------------------------------
                
                if accept: # Yes, we did accept the proposal. Woohoo!
                    
                    # Add the proposal to the new ensemble
                    new_ensemble[nw,:]          = copy.copy(proposal)
                    
                    self.chain[nw]              += [copy.copy(proposal)]
                    self.chain_logposterior[nw] += [logposterior_proposal]
                    self.chain_residuals[nw]    += [copy.copy(residuals)]
                    self.chain_duplicates[nw]   += [1]
                    
                    self.model                  = copy.deepcopy(self.proposed_model)
                    
                else: # No, we didn't accept the proposal. Boohoo! (Disclaimer: That's not actually a bad thing.)
                    
                    # Add the previous result to the new ensemble once more
                    new_ensemble[nw,:]          = copy.copy(self.chain_full[-1][nw,:])
                    
                    # This entry is a duplicate of the last
                    self.chain_duplicates[nw][-1]   += 1
                    
                    
            # Append the new ensemble to the chain ----------------------------
            self.chain_full     += [copy.copy(new_ensemble)]
        
        # Declare victory
        print('                                                     ')
        print('MCMC run finished.')
        
        # Clean up the chains, if desired
        if self.clean_chains:
            self.clean_up_chains()
        
        # Unify the chains, if desired
        if self.flat:
            self.flatten_chains()
        
    def unpack_parameters(self,params):
        
        """
        The model requires the parameters in a special shape. This function takes
        the packed parameters and flattens them into a vector.
        """
        
        import numpy as np
        import copy
        
        params_flat = np.zeros(self.C)
        
        counter     = 0
        for entry in params:
            if np.isscalar(entry):
                params_flat[counter]    = entry
                counter                 += 1
            else:
                length                  = len(entry)
                params_flat[counter:counter+length] = \
                    copy.copy(np.asarray(entry))
                counter                 += length
        
        return params_flat
    
    def pack_parameters(self,params_flat):
        
        """
        The model requires the parameters in a special shape. This function takes
        the flattened parameters and packs them into the shape required by the
        model.
        """
        
        import numpy as np
        import copy
        
        params      = []
        
        counter     = 0
        for entry in self.model.params:
            if np.isscalar(entry):
                params                  += [params_flat[counter]]
                counter                 += 1
            else:
                length                  = len(entry)
                params                  += [copy.copy(params_flat[counter:counter+length])]
                counter                 += length
                
        return params

    def evaluate_chain(self, z, subsampling = 1, 
        subsampling_type = 'random', mode = 'head', derivatives = 'all'):
        
        """
        Once you have obtained a valid MCMC chain, you can use this function to
        evaluate the model at points of your choosing. To help you, this function
        might do the trick:
            
        Parameters:
            z               	   - [vector]  : a vector of complex points inside the circular model domain where we wish to evaluate the model
            subsampling            - [scalar]  : a flat between 0 and 1 denoting the fraction of the chain we wish to evaluate. 1 means that all chain entries are evaluate, 0.1 means only 10% of the chain entries are evaluated, etc.
            subsampling_type       - [string]  : a flag which decides how we subsample the chain; 'random' means samples are selected at random, 'incremental' means samples are chosen in regular increments
            mode                   - [string]  : a flag which determines for what we evaluate; corresponds to 'mode' in the AEM toolbox' 'evaluate' function.
            derivatives            - [string]  : a flag which determines for what we evaluate; corresponds to 'derivatives' in the AEM toolbox' 'evaluate' function.

        """
        
        import copy
        import numpy as np
        from random import sample
        
        # Convert the evaluation points to complex variables, if necessary
        z           = self.complexify(z)
        
        # Create copies of the basic stuff
        z_base      = copy.deepcopy(z)
        model_base  = copy.deepcopy(self.model)
        
        # Copy the relevant chains
        chain_base              = copy.deepcopy(self.chain)
        chain_logposterior_base = copy.deepcopy(self.chain_logposterior)
        chain_residuals_base    = copy.deepcopy(self.chain_residuals)
        chain_duplicates_base   = copy.deepcopy(self.chain_duplicates)
        
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
        failure                 = []
        
        # Go through all walkers to find the progressbar maximum
        progressbar_maximum     = len(chain_base)
            
        progressbar_counter     = 0
            
        
        for idx,val in enumerate(chain_base):
            
            progressbar_counter += 1

            percent = ("{0:." + str(1) + "f}").format(100 * (progressbar_counter / float(progressbar_maximum)))
            filledLength = int(50 * progressbar_counter // progressbar_maximum)
            bar = '█' * filledLength + '-' * (50 - filledLength)
            print('\r%s |%s| %s%% %s' % ('Progress:', bar, percent, ''), end = '\r')
            
            # Propose a new sample (corresponding to a new model instance)
            model = copy.deepcopy(model_base)
            
            # Add a new proposal
            model.params = copy.deepcopy(self.pack_parameters(copy.copy(val)))
            
            # Update the proposed model
            model.update()
            
            result,error_flag = model.evaluate(
                copy.copy(z_base),
                mode=mode,
                derivatives=derivatives,
                return_error_flag=True)
            
            # Extract the complex potential at the query points
            if not error_flag:
                results     += [copy.copy(result)]
                failure     += [False]
            else:
                results     += [np.nan]
                failure     += [True]
         
        # Throw out all failed parts
        chain_base  = [entry for idx,entry in enumerate(chain_base) if not failure[idx]]
        chain_logposterior_base = [entry for idx,entry in enumerate(chain_logposterior_base) if not failure[idx]]
        chain_residuals_base    = [entry for idx,entry in enumerate(chain_residuals_base) if not failure[idx]]
        chain_duplicates_base   = [entry for idx,entry in enumerate(chain_duplicates_base) if not failure[idx]]
        results  = [entry for idx,entry in enumerate(results) if not failure[idx]]
        indices  = [entry for idx,entry in enumerate(indices) if not failure[idx]]
            

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
    
    def clean_up_chains(self, judgement_fraction = 0.1, threshold_offset = 0.1):
        
        import numpy as np
        import copy
        
        # Determine the evaluation period
        judgement_period    = int(self.walkers*judgement_fraction)
        
        # Pre-allocate space for the judgement vector
        judgement_matrix    = np.zeros((self.walkers,judgement_period))
        
        # Go through all walkers
        for nw in range(self.walkers):
        
            # Create a cumulative sum of duplicates from the end
            duplis  = np.flip(copy.copy(self.chain_duplicates[nw]))
            duplis_cumsum  = np.cumsum(duplis)
            
            # Find the first index where the cumulative sum exceeds the judgement period
            idx     = np.where(duplis_cumsum >= judgement_period)[0][0]
            
            # Then truncate the duplicates to this value
            over    = judgement_period - duplis_cumsum[idx]
            duplis  = duplis[:idx+1]
            duplis[-1]  += over
            
            # Write the loglikeihood entries into the storage matrix
            start   = 0
            for idx,end in enumerate(np.cumsum(duplis)):
                judgement_matrix[nw,start:end]   = \
                    self.chain_logposterior[nw][-idx-1]
                start   = end
        
        # Find the average logposterior value for identifying convergence
        convergence_threshold   = np.mean(judgement_matrix)
        
        # Offset the convergence threshold, if desired (default: 10%) 
        convergence_threshold   -= threshold_offset*abs(convergence_threshold)
        
#        print(convergence_threshold)
        
        # Now go through every chain and cut mask entries which did not converge
        for nw in range(self.walkers):
            
            converged   = False
            for idx in range(len(self.chain_logposterior[nw])):
        
                # If the logposterior crossed the convergence threshold, consider it converged
                if not converged and self.chain_logposterior[nw][idx] >= convergence_threshold:
                    converged   = True
                    
                if not converged:
                    self.chain_logposterior[nw][idx]    = np.nan
                
        return
    
    def flatten_chains(self, chain_first = True):
        
        import copy
        import numpy as np
        
        # Create copies of the chains for flattening
        copy_chain                  = copy.deepcopy(self.chain)
        copy_chain_logposterior     = copy.deepcopy(self.chain_logposterior)
        copy_chain_residuals        = copy.deepcopy(self.chain_residuals)
        copy_chain_duplicates       = copy.deepcopy(self.chain_duplicates)
        
        # Empty the chains
        self.chain                  = []
        self.chain_logposterior     = []
        self.chain_residuals        = []
        self.chain_duplicates       = []
        
        if chain_first:
            
            # Go through all walkers
            for nw in range(self.walkers):
                
                # Go through all entries of this walker
                for idx in range(len(copy_chain_logposterior[nw])):
                    
                    # If this entry hasn't been marked for non-convergence, append it
                    if not np.isnan(copy_chain_logposterior[nw][idx]):
                        
                        # Unify the individual chains
                        self.chain              += [copy.copy(copy_chain[nw][idx])]
                        self.chain_logposterior += [copy.copy(copy_chain_logposterior[nw][idx])]
                        self.chain_residuals    += [copy.copy(copy_chain_residuals[nw][idx])]
                        self.chain_duplicates   += [copy.copy(copy_chain_duplicates[nw][idx])]
            
        else: # Time first
            
            # To concatenate the chains in the correct order, we should first create
            # a cumulative sum of all duplicates
            duplicates_cumsum   = [np.cumsum(copy_chain_duplicates[nw]) for nw in range(self.walkers)]
            
            # Find the chain lengths for each walker
            chain_end           = np.asarray([entry[-1] for entry in duplicates_cumsum])
            
            # Initialize a counter for the chain length
            chain_counter   = 0
            
            # Initialize a repeater
            repeat          = True
            
            # Assemble the chain
            while repeat:
                
                # Go through all walkers
                for nw in range(self.walkers):
                    
                    index   = np.where(duplicates_cumsum[nw] >= chain_counter)[0]
                    
                    # The chain_counter hasn't reached this chain's end yet
                    if len(index) > 0:
                        
                        index   = index[0]
                        
                        if not np.isnan(copy_chain_logposterior[nw][index]):
                        
                            self.chain              += [copy.copy(copy_chain[nw][index])]
                            self.chain_logposterior += [copy.copy(copy_chain_logposterior[nw][index])]
                            self.chain_residuals    += [copy.copy(copy_chain_residuals[nw][index])]
                            self.chain_duplicates   += [copy.copy(copy_chain_duplicates[nw][index])]
                            
                            copy_chain_logposterior[nw][index]  = np.nan
                   
                if (chain_counter > chain_end).all():
                    repeat  = False
                   
                # Increment the chain counter
                chain_counter   += 1
            
        return
    
    def plot_chain(self):
        
        import matplotlib.pyplot as plt
        
        if not self.chain_initialized:
            raise Exception('Chain has not been evaluated yet.')
            
        plt.figure()
            
        if self.flat:
            
            all_logpost     = []
            for idx,val in enumerate(self.chain_duplicates):
                all_logpost     += [self.chain_logposterior[idx]]*val
            plt.plot(all_logpost)
            plt.xlabel('total chain length')
            plt.ylabel('unnormalized logposterior density')
            
        else:
            
            for nw in self.walkers:
            
                all_logpost     = []
                for idx,val in enumerate(self.chain_duplicates[nw]):
                    all_logpost     += [self.chain_logposterior[nw][idx]]*val
                plt.plot(all_logpost)
                
            plt.xlabel('chain length')
            plt.ylabel('unnormalized logposterior density')
        
        
    def plot_parameter_uncertainty(self,nbins = 20,**kwargs):
        
        """
        This function takes as input diagonal scatter data of dimensions [values,parameters]
        and plots them pairwise. Scatter points are assigned colors linearly along
        the chain, visualizing its course.
        """
        
        import numpy as np
        import scipy.stats
        import matplotlib.pyplot as plt
        import copy
        
        # Extract the labels; this is a bit of a hack, 
        labels  = []
        for idx,entry in enumerate(self.model.priors):
            
            # Check how often this variable occurs
            if not np.isscalar(self.model.params[idx]):
                occurence   = len(self.model.params[idx])
            else:
                occurence   = 1
            
            # Append this label that many times to the list
            if 'name' in list(entry.keys()):
                labels  += [entry['name']]*occurence
            else:
                labels  += [None]*occurence
                    
                
        # Convert the chain to a chain_length-by-C array
        data    = np.asarray(copy.deepcopy(self.chain))
        
        for row in range(self.C):
            
            for col in np.arange(row,self.C,1):
                
                if row == col:
                    
                    plt.subplot(self.C,self.C,self.C*row + col + 1)
                    
                    if labels[row] is not None:
                        plt.text(0.5,0.5,labels[row])
                    
                    plt.axis('off')
                    
                else:
                    
                    # Get a kernel density estimate
                    k = scipy.stats.kde.gaussian_kde(np.column_stack((data[:,col],data[:,row])).T)
                    
                    # Get the KDE grid if desired
                    xi, yi = np.mgrid[data[:,col].min():data[:,col].max():nbins*1j, data[:,row].min():data[:,row].max():nbins*1j]

                    # Interpolate the density
                    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
                    
                    # -------------------------------------------------------------
                    # Plot the upper triangular plot
                    # -------------------------------------------------------------
                    
                    plt.subplot(self.C,self.C,self.C*row + col + 1)
    
                    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), **kwargs)
                    # plt.contour(xi, yi, zi.reshape(xi.shape) ,levels=5)
                    
                    plt.gca().yaxis.tick_right()
                    plt.gca().xaxis.tick_top()
                    
                    if col == self.C-1:
                        plt.ylabel(labels[row])
                        plt.gca().yaxis.set_label_position("right")
                    else:
                        plt.gca().axes.yaxis.set_ticklabels([])
                        
                    if row == 0:
                        plt.xlabel(labels[col])
                        plt.gca().xaxis.set_label_position("top")
                    else:
                        plt.gca().axes.xaxis.set_ticklabels([])
                        
                    # -------------------------------------------------------------
                    # Plot the lower triangular plot
                    # -------------------------------------------------------------
                    
                    plt.subplot(self.C,self.C,self.C*col + row + 1)
                    
                    # plt.pcolormesh(yi.T, xi.T, zi.reshape(xi.shape).T, shading='gouraud', cmap=plt.cm.BuGn_r)
                    plt.pcolormesh(yi.T, xi.T, zi.reshape(xi.shape).T, **kwargs)
                    # plt.contour(yi.T, xi.T, zi.reshape(xi.shape).T,levels=5)
                    
                    # For the lower-triangular plot, row and col must be inverted
                    irow    = col
                    icol    = row
                    if icol == 0:
                        plt.ylabel(labels[irow])
                    else:
                        plt.gca().axes.yaxis.set_ticklabels([])
                        
                    if irow == self.C-1:
                        plt.xlabel(labels[icol])
                    else:
                        plt.gca().axes.xaxis.set_ticklabels([])
    
        return
        
    
#%%


#def weighted_avg_and_std(values, weights):
#    """
#    Return the weighted average and standard deviation.
#
#    values, weights -- Numpy ndarrays with the same shape.
#    """
#    
#    import numpy as np
#    
#    value_collector     = []
#    counter_collector   = []
#    
#    for nw in range(len(values)):
#        
#        for entry in values[nw]:
#            value_collector     += [entry]
#        
#        for entry in weights[nw]:
#            counter_collector   += [entry]
#            
#    values  = np.asarray(value_collector)
#    weights = np.asarray(counter_collector)
#    
#    print(values.shape)
#    print(weights.shape)
#        
#    average = np.average(values, weights=weights,axis = 0)
#    # Fast and numerically precise:
#    variance = np.average((values-average)**2, weights=weights,axis = 0)
#    
#    return (average, np.sqrt(variance))




def weighted_avg_and_std(values, weights):
     """
     Return the weighted average and standard deviation.

     values, weights -- Numpy ndarrays with the same shape.
     """
    
     import numpy as np
     
#     print(values.shape)
#     print(weights.shape)
    
     average = np.average(values, weights=weights,axis = 0)
     # Fast and numerically precise:
     variance = np.average((values-average)**2, weights=weights,axis = 0)
     return (average, np.sqrt(variance))
