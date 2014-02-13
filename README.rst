====================================
Welcome to the PyTrap package
====================================

This package is destined to physicists using or studying the
(`Electron Ion Beam Trap <http://pra.aps.org/abstract/PRA/v55/i3/pR1577_1>`_).
The main features are:

    - a fast (analytic) calculation of the **potential inside the trap** depending on the set of potentials
    - a **stability map** indicating which potentials lead to stable trapping
    - trajectories simulations.
    
Most of the proofs can be found in this `article`_ and for details see my `thesis`_.

Examples are the best way to dive into PyTrap.

.. index:: quickstart

Quick start
===========

Basic calculations
------------------

Calculations involving **PyTrap** can be performed
even without knowing anything about the Python_ programming language.

After installing this package and invoking the Python interpreter, calculations can be performed directly:
.. code-block:: python
   
   >>> from pytrap import *
   # Imagine you want to study Oxygen 4+ at 5.2 keV with V1=8000 V, V2=5850 V ...
   >>> setup = Setup(V = [8000, 5850, 4150, 1650, 3300], ma = 16, ne = 4, ener = 5200)
   # the conversion to SI units is automatic for mass and charge
   >>> print setup.m
   2.6761968e-26
   # creates an object Trap with the given setup
   >>> trap = Trap(setup)
   # gives the potential in Volts inside the trap at coordinates z = 0.15 m and r = 1 mm
   >>> trap.potential(0.15, 1e-3)
   1240.31320095 
   # calculates the oscillation frequency (in Hz) of the ions in the given setup (returns -1 if ions escaped)
   >>> print trap.fz()
   1187832.86276
   # Note that the ion is supposed to stay on the axis.
   # Note that each fuction is commented: trap.delta_r?
    
Stability of the trajectories
-----------------------------

The main difficulty of the EIBT is finding the good set of potentials in
order to have a stable trajectory of the ions. 
As described in the `article`_, the stability of the radial motion is described by the coefficient
delta:

    >>> trap.delta_r()
    1.919245181
    
Theoretically, the radial motion is stable if delta is lower than 1. 
So with this setup, there will be no trapping.
In order to find a good configuration of the potential, we fix V2, V3 and V4 and vary V1 and Vz.
This can be represented in a stability map.
The **Zgraph** class enables you to plot the stability map
of your setup. It is parallelized using `multiprocessing 
<http://docs.python.org/library/multiprocessing.html>`_.
For example, with the same ion at a lower energy, 4.2keV:

    >>> setup = Setup(V=[6000, 4850, 3150, 1650, 2000], ener = 4200)
    # the first range is V1 and the second is Vz and we want 20 points on each side (see Zgraph?)
    >>> graph = Zgraph(setup, (6000,8000), (3000, 5000), 20)
    # to launch the calculation:
    >>> graph.calc()
    # When it is finished, we can plot delta-1:
    >>> graph.plot()
    
should give you, where only the black areas are stable (delta-1<0):

.. figure:: _static/stab42.png
   :align: center
   :width: 500

   The radial stability map of oxygen 4+ at 4.2 kV. The black area show where trapping is possible.
     
Space charge and friction
-------------------------
All the previous functions can take space charge into account with the parameter **alpha**
and friction with **b**. See references for the definitions of these parameters.

Animation
---------
To see the influence of the space charge, for example, you can produce an animation using the
class **movie**:

    >>> mov = movie(setup, (6000,8000), (3000, 5000), 8, (0, 3), (0, 0), 5, path="")
    # with only 8 points along each side, 5 frames, we can vary alpha from 0 to 3:
    >>> mov.calc()
    # you can produce a .mov sequence (requires ffmpeg):
    >>> mov.mk_movie()
    # and don't forget to clean the data:
    >>> mov.clean()

.. figure:: _static/movie.gif
   :align: center
   :width: 500
   
    Here you see how space charge destroyes the stability areas.

Studying synchronization
------------------------

Since synchronization is beleived to obey an Hill's equation, we can use the same tools:

    >>> setup = Setup(V=[6000, 4850, 3150, 1650, 2000], ener = 4200, alpha=0.6)
    >>> Trap(setup).delta_s()
    1.45973959868
    # synchronization seems not possible with this setup
    # To make a map, just have to specify 'sync' in the mode:
    >>> graph = Zgraph(setup, (6000,8000), (3000, 5000), 20, mode='sync')
    # to launch the calculation:
    >>> graph.calc()
    # When it is finished, we can plot delta-1:
    >>> graph.plot()

.. figure:: _static/sync42.png
   :align: center
   :width: 500

   We see that synchronization can only occur in the lower part of the diagram.
   
Using the C-based version (cytrap)
=================================

The previous functions are purely developped in python. When a lot of particles need to be simulated and with high
precision, it is more suitable to implement the code in C with `GSL`_ library. Note that you should have GSL installed
and you should compile yourself cytrap before you can use the folowing functionalities (go in the cytrap directory and type Make).
Since it is very convenient to use only python, pytrap comes with wrapping functions:

Plotting Poincaré sections
--------------------------

Plotting Poincaré sections is simple and fast using pytrap:

    >>> setup = Setup(V = [6000, 5850, 4150, 1650, 4800], ma = 16, ne = 4, ener = 5200)
    >>> trap = Trap(setup)
    >>> trap.poincare(50,0.0003, 0.01, 60)
    #here we simulate 60 particles, initially located at z=0 and with r in [0.0003, 0.01] oscillating 50 times

.. figure:: _static/poincare.png
   :align: center
   :width: 500

   A Poincaré section in the EIBT.
   
Trajectory calculation
----------------------

To compute trajectories of ions in the EIBT use **trajectory**:

    >>> trap.trajectory(10, 0.0003, 0.003, 5)
    #here we simulate 5 particles, initially located at z=0 and with r in [0.0003, 0.003] oscillating for 10 units of time
    
.. figure:: _static/trajectory.png
   :align: center
   :width: 500
   
Installation instructions
=========================

The python part should not be difficult to install. Extract the pytrap directory from the tar.gz file,
place it where it is convenient, move to it and type:

    >>> sudo python setup.py install
    
From then you shoould be able to type:

    >>> import pytrap
    
without any error and use all the python based commands.

In order to user cython's commands, you should first get `GSL`_ installed. Then go into pytrap/cytrap directory
and type **make**. If the compilation succeeds, you can use cytrap.

Miscealenous
============

Source code
-----------

The `code <http://pypi.python.org/pypi/pytrap/>`_ is written in
python and uses: Matplotlib, Numpy and Scipy.
    
Available documentation
-----------------------

The :doc:`user_guide` details many of the features of this package.

The part :doc:`numpy_guide` describes how arrays of numbers with
uncertainties can be created and used.

The :doc:`tech_guide` gives advanced technical details.

In addition to this web documentation, the pydoc_ gives access to many
of the documentation strings included in the code.

.. index:: license
    
License
-------

This software is released under a **dual license**; one of the
following options can be chosen:

1. The `BSD license`_.
2. Any other license, as long as it is obtained from the creator of
   this package.

.. index:: support

Contact
-------

Please send feature requests, bug reports, or feedback to the creator
of :mod:`PyTrap`, `Alexandre Vallette`_.


Acknowledgments
---------------

The author wishes to thank Eric O. LEBIGOT (EOL) (author of the famous `uncertainties <http://packages.python.org/uncertainties/>`_ package) for is support on python and numpy.
The C-based part of this program was first developped by Khanh-Dang Nguyen Thu Lam, whom I thank a lot.


.. toctree::
   :hidden:
   :maxdepth: 1

   Overview <self>


.. _Python: http://python.org/
.. _article: http://hal.archives-ouvertes.fr/hal-00505782/fr/
.. _thesis: http://tel.archives-ouvertes.fr/tel-00605746/en/
.. _GSL: http://www.gnu.org/s/gsl/
.. _Alexandre Vallette: mailto:alexandre.vallette@spectro.jussieu.fr
.. _PyTrap package: http://

