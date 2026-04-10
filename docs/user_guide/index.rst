User guide
==========

.. contents::
   :local:
   :depth: 2

----

Two-step nesting workflow
-------------------------

Nesting couples a coarse *overall* model with a fine *detail* model.  The
detail model obtains its boundary conditions from the output of the overall
model.  Because the overall model must first be run before its output is
available, the workflow is split into two steps.

Step 1 — observation points
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nest1(overall, detail)`` inspects the boundary-point locations of the detail
model and adds matching observation points to the overall model.  After this
call the overall model can be saved and executed; it will write time-series
output at the observation-point locations.

Step 2 — boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nest2(overall, detail)`` reads the time-series output produced by the overall
model and sets it as boundary conditions on the detail model.  Depending on
the model combination, the boundary conditions are water levels, wave
parameters, or spectral data.

Supported model combinations
-----------------------------

The table below shows which overall/detail pairs are supported and what type
of boundary condition is transferred.

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 40

   * - Overall
     - Detail
     - BC type
     - Notes
   * - SFINCS
     - SFINCS
     - Water level
     - Supports ``filter_incoming`` and ``return_maximum``
   * - SFINCS
     - XBeach
     - Water level
     -
   * - SFINCS
     - BEWARE
     - Water level
     -
   * - HurryWave
     - HurryWave
     - Wave spectra
     -
   * - HurryWave
     - SFINCS
     - Wave parameters
     - Sets SnapWave boundary conditions (Hs, Tp, Wd, Ds)
   * - HurryWave
     - XBeach
     - Wave parameters
     - Supports deshoaling via ``option``
   * - HurryWave
     - BEWARE
     - Wave parameters
     -
   * - Delft3D-FM
     - Delft3D-FM
     - Water level
     -
   * - Delft3D-FM
     - SFINCS
     - Water level
     -
   * - Delft3D-FM
     - BEWARE
     - Water level
     -
   * - BEWARE
     - SFINCS
     - Water level + waves
     - Reads WL and R2_setup from BEWARE output

Model type resolution
---------------------

``nest1`` and ``nest2`` automatically detect the model type from the Python
class name of the object passed in.  The mapping is:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Class name
     - Resolved type
     - Package
   * - ``SfincsModel``
     - ``sfincs``
     - hydromt_sfincs
   * - ``HurrywaveModel``
     - ``hurrywave``
     - hydromt_hurrywave
   * - ``Delft3DFM``
     - ``delft3dfm``
     - cht_delft3dfm
   * - ``XBeach``
     - ``xbeach``
     - cht_xbeach
   * - ``BEWARE``
     - ``beware``
     - cht_beware

Any other class is resolved by lowercasing its name.

Observation-point naming
------------------------

Observation points created by ``nest1`` follow a consistent naming convention:

.. code-block:: text

   <prefix>_0001
   <prefix>_0002
   ...

The prefix defaults to the ``detail.name`` attribute, or ``"nest"`` if no name
is set.  You can override it with the ``obs_point_prefix`` parameter.

In ``nest2``, the same prefix is used to look up the corresponding stations in
the overall model output file.

Coordinate reference systems
-----------------------------

Both the overall and detail models carry a ``crs`` attribute (a ``pyproj.CRS``
object).  Coordinate transformations between models are handled automatically.

If a model object does not have its CRS set, you can pass it explicitly:

.. code-block:: python

   nest2(
       overall_model,
       detail_model,
       overall_crs="EPSG:4326",
       detail_crs="EPSG:32631",
   )

Water-level correction
----------------------

A constant offset can be added to all boundary water levels using the
``boundary_water_level_correction`` parameter (default ``0.0``).  This is
useful for correcting datum differences between models:

.. code-block:: python

   nest2(
       overall_model,
       detail_model,
       boundary_water_level_correction=0.15,  # add 15 cm
   )

Filtering and maximum values
-----------------------------

Two optional flags modify the behaviour of water-level nesting:

``filter_incoming``
   When ``True``, an incoming-wave filter is applied to the water-level
   signal.  This removes outgoing waves at the boundary and keeps only the
   incoming component.

``return_maximum``
   When ``True``, ``nest2`` returns a ``pandas.Series`` of peak water levels
   per time step instead of setting boundary conditions.

Writing boundary-condition files
---------------------------------

If ``bc_path`` is provided, ``nest2`` writes the boundary conditions to disk
after setting them on the detail model.  This is useful when the detail model
is run separately from the Python session:

.. code-block:: python

   nest2(
       overall_model,
       detail_model,
       output_path="/path/to/overall/output",
       bc_path="/path/to/detail/input",
   )
