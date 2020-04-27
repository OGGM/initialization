.. image:: _figures/logo.png


Initialization of past glacier states
------------------------------------

Reconstructing past glacier mass change is of interest for different applications, e.g. for quantifying their contribution to sea-level change.
One approach is to use a glacier model, forced by reconstructions of climate, to estimate past glacier states. However, glaciers respond to climate variability and
change with time lags between a few years and many centuries, and the backwards reconstruction is impeded by the non-linear interaction between glacier geometry and mass balance.

This repository host a new method, based on OGGM that estimates past glacier states.
The method was published in "The Cryosphere": `Initialization of a global glacier model based on present-day glacier geometry and past climate information: an ensemble approach <https://www.the-cryosphere.net/13/3317/2019/tc-13-3317-2019.html>`_.


Our method consists of 3 steps and the workflow is shown in the following figure:

- generation of glacier states,
- identification of glacier candidates, and
- their evaluation based on the misfit between the modelled and the observed surface elevation in the year of the observation.

.. image:: _figures/workflow.png

The current version was tested on glaciers located in the Alps and we reconstruct their state in 1850 by using synthetic experiments. Thus the reconstructed glacier states did
not represent the reality and need to be handled with caution. A further developement that will also allow the application to real-world cases is under developement and a new manuscript will the submitted soon.

In most cases, the resulting reconstruction is non-unique, as multiple initial states
converge towards the observed state in the year of observation.

Installation
------------

You can install this custom package with::

     pip install git+https://github.com/OGGM/reconstruction.git

However, what you'd probably like to do is to `fork <https://help.github.com/articles/fork-a-repo/>`_ this repository and use
it as a template for your own project. You can install it locally with::

    pip install -e .

Usage
-----
How to use this method is shown in this `example`_.
The usage of our method takes some time (hundreds of OGGM simulations are required to reconstruct the state of ONE year), and we recommend to use a super-computer environment for the reconstruction of
multiple glaciers.


Get in touch
------------

- View the source code `on GitHub`_.
- Report bugs or share your ideas on the `issue tracker`_.
- Improve the model by submitting a `pull request`_.
- Or you can always send us an `e-mail`_ the good old way.

.. _e-mail: jeis@uni-bremen.de
.. _on GitHub: https://github.com/OGGM/reconstruction
.. _issue tracker: https://github.com/OGGM/reconstruction/issues
.. _pull request: https://github.com/OGGM/reconstruction/pulls
.. _example: https://github.com/OGGM/reconstruction/tree/master/examples
