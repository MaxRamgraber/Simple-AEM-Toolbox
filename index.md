# A simple Analytic Element Method (AEM) toolbox

This repository contains the two Python3 toolboxes and a simple tutorial for the accompanying paper in [Water Resources Research](https://agupubs.onlinelibrary.wiley.com/journal/19447973) (current status: submitted).

## Installation

To use these toolboxes, simply download the Python files `toolbox_AEM.py` and `toolbox_MCMC.py'` (if desired) and copy them into your working directory. The MCMC toolbox is optional and only required if you also wish to use my MCMC implementation, the AEM toolbox can work as a standalone. Both toolboxes can be imported into your Python script by including the following code snippet:

```
from toolbox_AEM import *
from toolbox_MCMC import *
```

And that's it! 

## Troubleshooting

**Q: The model seems to create singularities, predicting very high or low water tables at certain isolated locations. What did I do wrong?**

_A_: This usually happens if the model attempts to evaluate the complex potential **Ω** directly on an element. This can happen because because two elements share a line segment or because one of the evaluation points lies on an a line segment. Make sure that the elements do not share direct borders, for example by offsetting them by a minuscule amount (e.g., 1E-10). I have implemented protections against this for some but not all elements: inhomogeneity elements, for example, are automatically shrunk by a negligible amount. Also make sure that no elements intersect.
