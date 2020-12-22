# A simple Analytic Element Method (AEM) toolbox

This repository contains the two Python3 toolboxes and a simple tutorial for the accompanying paper in [Water Resources Research](https://agupubs.onlinelibrary.wiley.com/journal/19447973) (current status: submitted).

## Installation

To use these toolboxes, simply download the Python files [`toolbox_AEM.py`](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_AEM.py) and [`toolbox_MCMC.py`](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_MCMC.py) (if desired) and copy them into your working directory. The MCMC toolbox is optional and only required if you also wish to use my MCMC implementation, the AEM toolbox can work as a standalone. Both toolboxes can be imported into your Python script by including the following code snippet:

```
from toolbox_AEM import *
from toolbox_MCMC import *
```

And that's it! 

## Tutorials

We have provided a few simple tutorials to help you familiarize yourself with the toolbox. The **basic tutorial** is available either as a [Python3 file](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/basic_tutorial) or as a [Jupyter Notebook](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Jupyter%20Notebooks). This tutorial covers some of the most basic elements for a simple, deterministic groundwater simulation.

For more advanced purposes, we have also added the **synthetic test case** of the accompanying paper as a [Python3 file](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/synthetic_reference). This scenario is a bit more complicated and does not only include a deterministic simulation, but also demonstrates the use of the (optional) [MCMC toolbox](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_MCMC.py) for Bayesian posterior inference. While we attempted to also make the synthetic test case approachable to beginners, we recommend starting with the basic tutorial first.

## Troubleshooting

**Q: The model seems to create singularities, predicting very high or low water tables at certain isolated locations. What did I do wrong?**

A: This usually happens if the model attempts to evaluate the complex potential Ω directly on an element. This can happen because because two elements share a line segment or because one of the evaluation points lies on an a line segment. Make sure that the elements do not share direct borders, for example by offsetting them by a minuscule amount (e.g., 1E-10). I have implemented protections against this for some but not all elements: inhomogeneity elements, for example, are automatically shrunk by a negligible amount. Also make sure that no elements intersect.
