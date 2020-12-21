"""
This tutorial serves to guide the user through the construction of the synthetic
test case.
"""

# First, we load a number a libraries we require in the process
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
import pickle
import copy
import os
import shapely
from matplotlib.colors import to_rgb
from matplotlib.gridspec import GridSpec 
from toolbox_AEM import *
from toolbox_MCMC import *

# Set a random seed
np.random.seed(0)

# We will use a custom colormap for plotting; this is merely a cosmetic solution
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", 
    ["xkcd:midnight blue",
     "xkcd:cerulean",
     "xkcd:kermit green",
     "xkcd:goldenrod",
     "xkcd:orangish red",
     "xkcd:brick red"])

# Close all figures we may have created before
plt.close('all')

# Define the model domain
domain_center   = [0,0]
domain_radius   = 250

# We will require the geometry of the elements we define later on; we could just
# define them manually - as a numpy array with two columns - but we have prepared
# some of the necessary files in advance and just load them in
river           = pickle.load(open('synth_case_river.p','rb'))
noflow          = pickle.load(open('synth_case_noflow.p','rb'))
noflow_cs       = pickle.load(open('synth_case_noflow_cellselect.p','rb'))
inhomogeneity   = pickle.load(open('synth_case_inhomogeneity.p','rb'))
inhomogeneity2  = pickle.load(open('synth_case_inhomogeneity2.p','rb'))

#%%

# =============================================================================
# MODEL CREATION
# =============================================================================

# We start by initializing the model domain, specifying its bottom elevation
# (525m) and the aquifer type ('unconfined'). The center and extent of the 
# circular model domain are also specified here.
ml = Model(
    head_offset     = 525,
    aquifer_type    = 'unconfined',
    domain_center   = domain_center,
    domain_radius   = domain_radius)

# For the Möbius base flow, we require a bit of preparation. Since we don't
# just wish to evaluate a single model, but do Bayesian inference with MCMC, 
# we need to provide all elements which have random variables with prior and 
# proposal probability distributions.

# We start by defining the correlation matrix for the proposal of the Möbius
# control point rotations (r). We consider positive correlation between them.
cov         = np.zeros((3,3))
cov[0,2]    = 0.25
cov[2,0]    = 0.25
cov[1,0]    = 0.5
cov[0,1]    = 0.5
cov[1,2]    = 0.5
cov[2,1]    = 0.5
np.fill_diagonal(cov,1)

# The hydraulic conductivity is defined and used in canonical/normal space, but
# we may wish to define its pdf and updates in logarithmic space. To facilitate 
# this conversion we must define a 'converter' and a 'deconverter' function.
def converter_log10(x):
    from numpy import log10
    return log10(x)
def deconverter_log10(x):
    return 10**x

# We set all four of the Möbius base flow element's variables  (r, head_min, 
# head_max, k) as random. We define this by adding their strings to a list.
variables_mb    = ['r','head_min','head_max','k']

# For every entry in 'variables', we must provide a dictionary specifying its
# prior probability density function. This is done as a list of dictionaries.
# We are also free to specify 'limits', for which three options are available:
#   - strict numerical values for strict numerical limits
#   - None to set an upper or lower limit as undefined
#   - a string of another random variable of this element to set it as a 
#     dynamic limit
# Note that the probability density and updates of 'k' are defined in logarithmic
# space, so we specify its 'convert' and 'deconverter' in its prior dictionary.
priors_mb       = [
    
    {'name'             : 'Moebius r',
      'distribution'    : 'vonmises',
      'kappa'           : 1E-4,
      'loc'             : np.asarray([-0.4,0.1,0.8])*np.pi,
      'scale'           : 1},
    
    {'name'             : 'Moebius head_min',
      'distribution'    : 'norm',
      'loc'             : 537,
      'scale'           : 2.5,
      'limits'          : [530, 'head_max']},
    
    {'name'             : 'Moebius head_max',
      'distribution'    : 'norm',
      'loc'             : 540,
      'scale'           : 2.5,
      'limits'          : ['head_min', None]},
    
    {'name'             : 'hydraulic conductivity',
      'distribution'    : 'norm',
      'converter'       : converter_log10,
      'deconverter'     : deconverter_log10,
      'loc'             : -4,
      'scale'           : 0.5,
      'limits'          : [-7,-1]}]

# Equivalently, we create a list of dictionaries for the MCMC proposal 
# distributions.
proposals_mb    = [
    
        {'distribution'     : 'multivariate_normal',
          'mean'            : np.zeros(3),
          'cov'             : cov*(0.5)**2,
          'circular'        : True},
        
        {'distribution'     : 'norm',
          'loc'             : 0,
          'scale'           : 0.5},
        
        {'distribution'     : 'norm',
          'loc'             : 0,
          'scale'           : 0.5},
        
        {'distribution'     : 'norm',
          'loc'             : 0,
          'scale'           : 0.25,
          'converter'       : converter_log10,
          'deconverter'     : deconverter_log10}]

# Now we can define the base flow element.
ElementMoebiusBase(
    model           = ml,
    r               = np.asarray([0,0.5,1])*np.pi,
    head_min        = 537,
    head_max        = 540,
    k               = 1E-4,
    variables       = variables_mb,
    priors          = priors_mb,
    proposals       = proposals_mb)

# -----------------------------------------------------------------------------

# Now we get to define the river, implemented as a prescribed head boundary of
# uncertain connectivity. As a consequence, we start off by marking this 
# variable as uncertain.
variables_hb    = ['connectivity']

# Similarly to the Möbius base, we define it's prior...
priors_hb       = [
    
       {'name'              : 'river connectivity',
        'distribution'      : 'beta',
        'a'                 : 2,
        'b'                 : 4}]

# ...and proposal distribution.
proposals_hb    = [
        {'distribution'     : 'multivariate_normal',
          'mean'            : np.zeros(4),
          'cov'             : np.identity(4)*(0.1)**2}]

# Note that we can define connectivity in two ways: spatially uniform or 
# spatially heterogeneous. In the former case, 'connectivity' is simply a 
# scalar value. In the latter case, 'connectivity' is a vector of length N > 2,
# where N is the number of connectivity nodes. These are placed by normalized
# distances 0 <= d <= 1 along the element in the 'connectivity_normdist' 
# variable. We subsequently assign the connectivity to each segment through
# linear interpolation between the nodes.

# Mind also that the variable 'segments' allows us to subdivide the element 
# more finely than the specified geometry ('river' and 'river_ht). This process
# tries to create segments of equal length. The more segments, the better the
# resolution of the element (no-flow boundaries, for example, may not work with
# an insufficient resolution), but increases computational cost.

# Now we can define the Head Boundary
ElementHeadBoundary(
    model           = ml, 
    line            = river,
    line_ht         = np.linspace(535,545,river.shape[0]),
    segments        = 100,
    connectivity_normdist = np.asarray([0.,0.33,0.66,1.]),
    connectivity    = np.asarray([0.75,0.75,0.75,0.75]),
    variables       = variables_hb,
    priors          = priors_hb,
    proposals       = proposals_hb)

# -----------------------------------------------------------------------------

# For the inhomogeneity in the northern part of the model domain, we only have
# a single uncertain variable: the hydraulic conductivity. It is defined the
# same way as we did in the Möbius base flow element.

# First marking it as uncertain...
variables_in    = ['k']

# ...then defining its prior...
priors_in       = [
    
        {'name'             : 'hydraulic conductivity 1',
          'distribution'    : 'norm',
          'converter'       : converter_log10,
          'deconverter'     : deconverter_log10,
          'loc'             : -5,
          'scale'           : 0.5,
          'limits'          : [-7,-1]}]

# ...and definings its proposal distribution.
proposals_in    = [
    
        {'distribution'     : 'norm',
          'loc'             : 0,
          'scale'           : 0.25,
          'converter'       : converter_log10,
          'deconverter'     : deconverter_log10}]

# With this, we can define the inhomogeneity. 
ElementInhomogeneity(
    model           = ml,
    k               = 1E-5,
    polygon         = inhomogeneity,
    segments        = 100,
    variables       = variables_in,
    priors          = priors_in,
    proposals       = proposals_in)

# -----------------------------------------------------------------------------

# Finally, we can add the extraction well and the no-flow boundary in the 
# south. Since we set no variable of either element as random, their definition
# is rather straightforward.

# First add the extraction well...
ElementWell(
    model           = ml,
    zc              = [-50,30],
    rw              = 0.1,
    strength        = -0.0025)

# ...then the no-flow boundary.
ElementNoFlowBoundary(
    model           = ml,
    line            = noflow,
    segments        = 100)

#%%

# With the model defined, we can now prepare its evaluation. One of the big 
# advantages of AEM is that we can define the evaluation points (or grid) after
# we have defined the flow-relevant elements. Since we focus on 
 
# Define the points at which we wish to evaluate the model. Since our model 
# domain is circular, we are interested in points within the disk. We could, 
# for example, simply create a meshgrid and mask all points inside. For this
# study, however, we have included a function which automatically creates a 
# collection of equi-distant points within a disk.
XY = equidistant_points_in_circle(
    rings           = 101,
    radius          = domain_radius)

# We may not want to evaluate the points inside the (closed) no-flow element,
# so we mask those from our collection of evaluation points.
linepath    = matplotlib.path.Path(noflow_cs)
indices     = linepath.contains_points(XY)
XY          = XY[~indices,:]

# With this, we are ready to evaluate the AEM. Simply call the model.evaluate
# function with the specified evaluation points (XY) and the desired mode. 
# Possible modes are 'head' (hydraulic head), 'potential' (complex potential),
# and 'gradient' (which returns the gradient, with the x coordinate being the
# real and the y coordinate the imaginary part; the variable 'derivatives' 
# defines whether the user wishes to obtain the derivatives of the hydrualic
# potential ('phi') or the stream function ('psi')).
z = ml.evaluate(XY,mode='head') 

# We have also added a rudimentary function to plot the elements the user 
# defined before.
ml.plot()

# Now let's add the contours of our model results.
plt.tricontour(XY[:,0],XY[:,1],np.real(z),levels=21)
plt.colorbar(label='hydraulic head in $m$')

# And make sure the axes are equally spaced.
plt.axis('equal')

# Tricontour will also interpolate over the no-flow boundary; mask this area
plt.fill(noflow_cs[:,0],noflow_cs[:,1],facecolor='w',zorder=2)


#%%

# =============================================================================
# BAYESIAN INFERENCE WITH MCMC
# =============================================================================

# Now we may can start the Bayesian inference with MCMC. Towards this end, we 
# must first define the observations we use for the determination of the 
# likelihood. The observations are defined similarly to the priors and 
# proposals, as a list of dictionaries. Each dictionary contains a (optional)
# name, a location in terms of complex coordinates, and an associated head 
# measurement.
observations    = [
    
    {'name'         : 'Observation 1',
     'location'     : -50+30j,
     'head'         : 535.36},
    
    {'name'         : 'Observation 2',
     'location'     : 10-100j,
     'head'         : 538.71},
    
    {'name'         : 'Observation 3',
     'location'     : 100+115j,
     'head'         : 539.80},
    
    {'name'         : 'Observation 4',
     'location'     : -175+85j,
     'head'         : 536.93}]

# We must also define a dictionary for the likelihood distribution. In this
# implementation, we assume that the same likelihood function is used for all
# observations. If we use a Gaussian pdf ('norm'), we need only define the 
# scale. The mean ('loc') will be defined by the model prediction, and the 
# evaluation point ('x') by the observation ('head' in the 'observations'
# dictionary).
likelihood_dictionary   = {
        'distribution'      : 'norm',
        'scale'             : 0.15}

# With this, we can start defining the MCMC class. The chain length defines for
# how many entries the MCMC will run - the more entries, the more samples are
# drawn from the posterior, but the higher the computational effort.
# The last three variables ('adapt_frequency', 'adapt_proposal', and 
# 'acceptance_target') are only used when 'adapt_proposal' is set to True, and
# define how the proposal is dynamically adapted. The chain will adapt the 
# proposal every 'adapt_frequency' step, and increase or decrease the proposal
# covariance in order to reach 'acceptance_target' - a user-specified ratio of
# accepted proposals.
mcmc = MCMC(
    model                   = ml,
    chain_length            = 10000,
    observations            = observations,
    likelihood_dictionary   = likelihood_dictionary,
    adapt_frequency         = 100,
    adapt_proposal          = True,
    acceptance_target       = 0.2)

# # With the MCMC base defined, let us start the algorithm and store the results.
# mcmc.start()
# pickle.dump(mcmc,open('mcmc.p','wb'))

# Alternatively, we can simply load the results we obtained previously:
mcmc = pickle.load(open('mcmc.p','rb'))

# Let us plot the unique entries in the logposterior chain
plt.figure()
plt.plot(mcmc.chain_logposterior)


# # As of now, MCMC has only evaluated the model at the observation wells (recall
# # that AEM does not require you to evaluate it everywhere). If we wish to plot
# # the uncertainty in the flow fields, we should evaluate the model everywhere.
# # This can be done with the function 'evaluate_chain'.
# # For the evaluation, we can discard the "burn-in" ('cutoff'). Looking at the 
# # logposterior density, a good point would be at around 250 entries. One should
# # select a location for the burn-in after which the chain appears 'normal'.
# # Since it can be computationally very expensive to evaluate every single entry
# # in the chain, we can further subsample the ensemble. 'subsampling' of 0.1 
# # means that only 10% of the entries in the chain are being evaluated, and the
# # 'subsampling_type' can be either 'incremental' (evaluate every 10th entry) or
# # 'random' (10% randomly selected entries are evaluated)
# evaluation_dictionary = mcmc.evaluate_chain(
#     z                       = XY, 
#     cutoff                  = 250, 
#     subsampling             = 0.1, 
#     subsampling_type        = 'incremental')
# pickle.dump(evaluation_dictionary,open('evaluation_dictionary.p','wb'))

# Alternatively, we can simply load the results we obtained previously:
evaluation_dictionary = pickle.load(open('evaluation_dictionary.p','rb'))


#%%

# =============================================================================
# PLOTTING HYDRAULIC HEAD UNCERTAINTY
# =============================================================================

# To plot the hydraulic head uncertainty, let us first calculate the mean and
# standard deviation of the hydraulic head entries. Mind that the chain we 
# we evaluated might still contain duplicates. However, we stored it only with
# the unique samples, and kept count of how often each sample is retained in
# 'evaluation_dictionary['chain_duplicates']'.
mean, std = weighted_avg_and_std(
    values      = np.real(evaluation_dictionary['results']), 
    weights     = evaluation_dictionary['chain_duplicates'])

# Now start plotting the average hydraulic head distribution
plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
for element in ml.elementlist[1:]:
    element.plot()
plt.tricontour(
    np.real(evaluation_dictionary['z']),
    np.imag(evaluation_dictionary['z']),
    mean,
    levels  = 25,
    cmap    = cmap,
    zorder  = 1)
plt.axis('equal')
plt.colorbar(label='hydraulic head $\mu$ in $m$')
plt.xlabel('$x$ coordinate in $m$')
plt.ylabel('$y$ coordinate in $m$')

# Again, mask the are behind the no-flow boundary
plt.fill(noflow_cs[:,0],noflow_cs[:,1],facecolor='w',zorder=2)

# And scatter the observation well locations
for entry in observations:
    
    plt.scatter(
        np.real(entry['location']),
        np.imag(entry['location']))
    plt.text(
        np.real(entry['location']),
        np.imag(entry['location']),
        entry['name'])

# Plot the hydraulic head standard deviation in the second subplot
plt.subplot(1,2,2)
plt.tricontourf(
    np.real(evaluation_dictionary['z']),
    np.imag(evaluation_dictionary['z']),
    std,
    levels  = 25,
    cmap    = cmap,
    zorder  = 1)
plt.axis('equal')
plt.colorbar(label='hydraulic head $\sigma$ in $m$')
plt.xlabel('$x$ coordinate in $m$')
plt.ylabel('$y$ coordinate in $m$')

# Again, mask the are behind the no-flow boundary
plt.fill(noflow_cs[:,0],noflow_cs[:,1],facecolor='w',zorder=2)

# Plot the elements
for element in ml.elementlist[1:]:
    element.plot()

# And scatter the observation well locations
for entry in observations:

    plt.scatter(
        np.real(entry['location']),
        np.imag(entry['location']))
    plt.text(
        np.real(entry['location']),
        np.imag(entry['location']),
        entry['name'])


#%%

# raise Exception

# =============================================================================
# PATHLINE TRACING
# =============================================================================

# # This is a rudimentary code snippet for pathline tracing of the MCMC results.
# # It is quite computationally inefficient (and could be readily parallelized),
# # so use it with caution.

# # Create a list for the flow tracing pathways around the wells
# traces  = []

# # Count the unique entries in the MCMC
# C   = len(mcmc.chain)

# # Then go through all of them
# for c in range(C):
    
#     # Every 10th step, save the current traces evaluated and reset the 'traces'
#     # list to prevent it from becoming to large for the RAM.
#     if (c+1)%10 == 0:
#         pickle.dump(traces,open('traces_'+str(c).zfill(4)+'.p','wb'))
#         traces  = []
    
#     # Print the current progress
#     print(str(c).zfill(4)+' of '+str(C).zfill(4)+' possibilities explored.')
    
#     # Load the parameters of this entry in the MCMC again, and update the model
#     ml.params = copy.deepcopy(mcmc.chain[c])
#     ml.update()
    
#     # Create a random rotation offset for the ring start positions
#     rotation_offset = np.random.uniform(-np.pi,np.pi,1)

#     # Generate 36 points around the well
#     ring    = np.column_stack((
#         np.cos(np.linspace(0,2*np.pi,13)+rotation_offset)[:-1],
#         np.sin(np.linspace(0,2*np.pi,13)+rotation_offset)[:-1] ))
    
#     # Create a dummy to store all traces
#     dummy   = []
    
#     # And trace from each starting point in sequence; 'p' is the offset from the
#     # extraction well
#     for p in ring:
        
#         # Trace the gradients
#         points = ml.trace_gradient(
#             p                   = np.asarray([-50,30])+p,
#             stepsize            = 10,
#             well_snap_distance  = 0.1)
        
#         # And store the trace
#         dummy.append(copy.copy(points))
    
#     # Then save all traces for this result
#     traces.append(copy.deepcopy(dummy))

# # Save the final part of the traces as well
# pickle.dump(traces,open('traces_'+str(c).zfill(4)+'.p','wb'))

# # Now concatenate the traces into a more concise shape
# root_directory = os.path.dirname(os.path.realpath(__file__))
# results = []
# results += [each for each in os.listdir(root_directory) if each.startswith('traces_')]

# # Load the traces, combine them into one
# traces  = []
# for res in results:
#     temp    = pickle.load(open(res,'rb'))
#     for entry in temp:
#         traces  += copy.copy(entry)

# # Then save the results
# pickle.dump(traces,open('combined_traces.p','wb'))

# # And plot the traces
# for idx,entry in enumerate(traces):
#     if idx%15 == 0:
#         plt.plot(entry[:,0],entry[:,1],color='xkcd:silver',alpha=0.01)
        
# If the model has already simulated the traces (caution: this can take one and
# a half days with my current, unoptimized, non-parallel code), load the
# results instead
traces     = pickle.load(open('combined_traces.p','rb'))

# Plot the results
plt.figure()
domain  = np.column_stack((
    np.cos(np.linspace(0,2*np.pi,361))*domain_radius+domain_center[0],
    np.sin(np.linspace(0,2*np.pi,361))*domain_radius+domain_center[1] ))
plt.fill(noflow_cs[:,0],noflow_cs[:,1],facecolor='w',zorder=2)
plt.plot(domain[:,0], domain[:,1], color = 'xkcd:grey', linewidth = 2, zorder = -2)
plt.plot(ml.elementlist[4].line[:,0],ml.elementlist[4].line[:,1], color = 'xkcd:medium grey', linewidth = 2)
mcmc.model.elementlist[3].plot()
mcmc.model.elementlist[1].plot()
plt.axis('equal')
maxdist = 0
for idx,tr in enumerate(traces):
    if len(tr) > maxdist:
        maxdist = len(tr)
    if idx%24 == 0:
        plt.plot(tr[:,0],tr[:,1],color='xkcd:silver',alpha=0.01*mcmc.chain_duplicates[int(idx/12)],linewidth=0.5)