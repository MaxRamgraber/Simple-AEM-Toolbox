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


<img align="left" src="https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/images/01_uniform.png" width="15%">

### Uniform base flow
To provide a background potential or unidirectional regional flow, the simplest option is to use uniform flow, specified by the AEM toolbox' object `ElementUniformBase`. This element requires the specification of a direction in radians, a minimum and maximum hydraulic head, and a background hydraulic conductivity.
<br /><br />

<img align="right" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/02_moebius.png" width="15%">

### Möbius base flow
Möbius base flow provides a way to implement more intricate regional flow, allowing for curvature, divergence, and convergence. This type of flow is specified by the AEM toolbox' object `ElementMoebiusBase`. This element requires the specification of three control points' direction in radians, a minimum and maximum hydraulic head, and a background hydraulic conductivity.
<br /><br />

<img align="left" src="https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/images/03_extraction_well.png" width="15%">

### Extraction or injection well
Injection or extraction wells - or any other type of pointwise flow - can be implemented using the `ElementWell` object. This element requires the specification of a position, a positive or negative strength value, and a well radius. Alternatively to a strength value, this element can also adjust its strength to induce a desired drawdown on the flow field.
<br /><br />

<img align="right" src="https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/images/04_inhomogeneity.png" width="15%">

### Polygonal inhomogeneity
Zonal inhomogeneities in the aquifer's hydraulic conductivity can be represented using the `ElementInhomogeneity` object. This element requires the specification of a hydraulic conductivity value, as well as a closed or open polygon defining its extent.
<br /><br />

<img align="left" src="https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/images/05_fixed_head_boundary.png" width="15%">

### Prescribed head boundary / River
Prescribed head boundary conditions or rivers can be created using the `ElementHeadBoundary` object. This line element enforces a specified hydraulic head along its path. It requires the specification of its vertices and corresponding head values. It also allows for the implementation of a uniform or spatially varying connectivity value which can limit its influence on the water table.
<br /><br />

<img align="right" src="https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/images/06_no_flow_boundary.png" width="15%">

### No-flow boundary
No-flow boundaries from sheet pile walls or impermeable formations can be created using the `ElementNoFlowBoundary` object. This line element requires only the specification of the vertices along its path and can be either closed or open.
<br /><br />

<img align="left" src="https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/images/07_area_sink.png" width="15%">

### Area sources or sinks
Areal recharge or water extraction can be represented using the `ElementAreaSink` object. This element adds or removes water according to a specified range inside its polygon. It requires the specification of its polygons and a positive or negative strength value. Not that the stream function component of the complex potential is not valid inside an area source or sink.
<br /><br />

## Troubleshooting

**Q: The model seems to create singularities, predicting very high or low water tables at certain isolated locations. What did I do wrong?**

A: This usually happens if the model attempts to evaluate the complex potential Ω directly on an element. This can happen because because two elements share a line segment or because one of the evaluation points lies on an a line segment. Make sure that the elements do not share direct borders, for example by offsetting them by a minuscule amount (e.g., 1E-10). I have implemented protections against this for some but not all elements: inhomogeneity elements, for example, are automatically shrunk by a negligible amount. Also, you should make sure that no inhomogeneities or no-flow boundaries intersect.

**Q: There are still strange artefacts along my no-flow boundaries or inhomogeneities. What happened?**

A: If you have tried the solutions in the answer above and the issue persists, try increasing the resolution of the element by increasing the element's `segments`. Most of the constant-strength line elements I used here require sufficient resolution to induce the desired effect. It is difficult to predict how large this resolution should be in advance.
