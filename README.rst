.. image:: _figures/logo.png


Reconstruct estimated glacier states
------------------------------------

Reconstructing past glacier change is of interest for different applications, e.g. for quantifying their contribution to sea-level change.
One approach is to use a glacier model, forced by reconstructions of climate, to estimate past glacier states. However, glaciers respond to climate variability and
change with time lags between a few years and many centuries, and the backwards reconstruction is impeded by the non-linear interaction between glacier geometry and mass balance.

This repository host a new method, based and developed for the usage of OGGM, that estimates past glacier states

Our method consists of 3 steps:

- generation of glacier states
- identification of glacier candidates, and
- their evaluation based on the misfit between the modelled and the observed surface elevation at the year of the observation.

We tested our method on glaciers located in the Alps and reconstruct their state in 1850. In most cases, the resulting reconstruction is non-unique, as multiple initial states
converge towards the observed state in the year of observation.


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
.. _example: https://github.com/OGGM/reconstruction/example
