# A simple Analytic Element Method (AEM) toolbox

This repository contains the two Python3 toolboxes and a simple tutorial for the accompanying paper in [Water Resources Research](https://agupubs.onlinelibrary.wiley.com/journal/19447973) (current status: under review).

## Installation

To use these toolboxes, simply download the Python files [`toolbox_AEM.py`](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_AEM.py) and [`toolbox_MCMC.py`](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_MCMC.py) (if desired) and copy them into your working directory. The MCMC toolbox is optional and only required if you also wish to use my MCMC implementation, the AEM toolbox can work as a standalone. Both toolboxes can be imported into your Python script by including the following code snippet:

```
from toolbox_AEM import *
from toolbox_MCMC import *
```

And that's it! 

## Tutorials

I have provided a few simple tutorials to help you familiarize yourself with the toolbox. The [**basic tutorial**](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/Tutorial%2001%20Basic%20AEM) covers the fundamentals of constructing a deterministic flow model with this toolbox, and is available as both a [Python file](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2001%20Basic%20AEM/basic_tutorial.py) and a [Jupyter Notebook](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2001%20Basic%20AEM/basic_tutorial.ipynb). 

The [**advanced tutorial**](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/Tutorial%2002%20Uncertainty%20Estimation) expands the scope of the basic tutorial by showing how to prepare the Analytic Element Model for uncertainty quantification. This tutorial is also available as both a [Python file](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2002%20Uncertainty%20Estimation/uncertainty_estimation_example.py) and a [Jupyter Notebook](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2002%20Uncertainty%20Estimation/uncertainty_estimation_example.ipynb).

For users with an interest in reproducing the some or all of the results in accompanying manuscript, I have also [uploaded the Python files](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Manuscript%20files) required to reproduce the figures in the main manuscript and the supporting information.

## Elements


<img align="left" src="/images/01_uniform.png" width="15%">

### Uniform flow
Basic uniform flow 










## Troubleshooting

**Q: The model seems to create singularities, predicting very high or low water tables at certain isolated locations. What did I do wrong?**

A: This usually happens if the model attempts to evaluate the complex potential Ω directly on an element. This can happen because because two elements share a line segment or because one of the evaluation points lies on an a line segment. Make sure that the elements do not share direct borders, for example by offsetting them by a minuscule amount (e.g., 1E-10). I have implemented protections against this for some but not all elements: inhomogeneity elements, for example, are automatically shrunk by a negligible amount. Also, you should make sure that no inhomogeneities or no-flow boundaries intersect.

**Q: There are still strange artefacts along my no-flow boundaries or inhomogeneities. What happened?**

A: If you have tried the solutions in the answer above and the issue persists, try increasing the resolution of the element by increasing the element's `segments`. Most of the constant-strength line elements I used here require sufficient resolution to induce the desired effect. It is difficult to predict how large this resolution should be in advance.
