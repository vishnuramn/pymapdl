"""
.. _ref_how_to_add_an_example_reference_key:

Statically Indeterminate Reaction Force Analysis
------------------------------------------------
Problem Descrption:
 - A prismatical bar with built-in ends is loaded axially at two
   intermediate cross sections.  Determine the reactions :math:`R_1`
   and :math:`R_2`.

Reference:
 - S. Timoshenko, Strength of Materials, Part I, Elementary Theory and
   Problems, 3rd Edition, D. Van Nostrand Co., Inc., New York, NY, 1955,
   pg. 26, problem 10.

Analysis Type(s):
 - Static Analysis ``ANTYPE=0``

Element Type(s):
 - 3-D Spar (or Truss) Elements (LINK180)

.. image:: ../../_static/vm1_setup.png
   :width: 400
   :alt: VM1 Problem Sketch

Material Properties
 - :math:`E = 30 \cdot 10^6 psi`

Geometric Properties:
 - :math:`a = b = 0.3`
 - :math:`l = 10 in`

Loading:
 - :math:`F_1 = 2*F_2 = 1000 lb`

Analytical Equations:
 - :math:`P = R_1 + R_2` where :math:`P` is load.
 - :math:`\frac{R_2}{R_1} = \frac{a}{b}`
   Where :math:`a` and :math:`b` are the ratios of distances between
   the load and the wall.
"""

from ansys.mapdl.core import launch_mapdl

# start MAPDL and enter the pre-processing routine
mapdl = launch_mapdl()
mapdl.verify('VM47')

# Torsional Frequency of a Suspended Disk
mapdl.prep7()
mapdl.antype('MODAL')
mapdl.modopt('LANB', 1)

###############################################################################
# Animations
# ~~~~~~~~~~
# You can even create animations.  See :ref:`ref_pyvista_mesh` for an example.
# Incidentally that is also how you link to another example (via `ref_pyvista_mesh`).
#
#
# Making a Pull Request
# ~~~~~~~~~~~~~~~~~~~~~
# Once your example is complete and you've verified builds locally, you can make a pull request (PR).
# Branches containing examples should be prefixed with `doc/` as per the branch
# naming conventions found here: :ref:`contributing`.
#
# Note that you only need to create the python source example (*.py).  The jupyter
# notebook, the example html and the demo script will all be auto-generated via ``sphinx-gallery``.


mapdl.et(1, 'COMBIN14', '', '', 1)
mapdl.et(2, 'mASS21', '', '', 3)
mapdl.r(1, 4.8)
mapdl.r(2, 1, 0.30312)
mapdl.n(1)
mapdl.n(2, '', '', -1)
mapdl.e(1, 2)
mapdl.type(2)
mapdl.real(2)
mapdl.e(2)
mapdl.d(1, 'ALL')
mapdl.d(2, 'UX', '', '', '', '', 'UY', 'UX')
mapdl.finish()
mapdl.slashsolu()
mapdl.solve()
frequency = mapdl.get('FREQ', 'MODE', 1, 'FREQ')


results = f"""
------------------- VM47 RESULTS COMPARISON ---------------

     |   TARGET   |   Mechanical APDL   |   RATIO

         0.6333            {frequency:1.4f}          {frequency/0.6333:1.4f}

-----------------------------------------------------------

"""

print(results)

