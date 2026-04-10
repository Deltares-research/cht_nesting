API reference
=============

.. contents::
   :local:
   :depth: 2

----

Top-level functions
-------------------

The two entry points are re-exported from the package root:

.. code-block:: python

   from cht_nesting import nest1, nest2

nest1
~~~~~

.. py:function:: cht_nesting.nest1(overall, detail, option=None, obs_point_prefix=None)

   Add observation points to the overall model at the boundary-point locations
   of the detail model.

   The correct implementation is selected automatically based on the class
   names of *overall* and *detail*.

   :param overall: The coarse model that receives the new observation points.
   :type overall: Any
   :param detail: The fine model whose boundary points are used as observation
       point locations.
   :type detail: Any
   :param option: Reserved for future use.
   :type option: str or None
   :param obs_point_prefix: Prefix for observation-point names.  Defaults to
       ``detail.name`` or ``"nest"``.
   :type obs_point_prefix: str or None
   :returns: ``True`` if observation points were added, ``False`` if the
       combination is not supported.
   :rtype: bool

nest2
~~~~~

.. py:function:: cht_nesting.nest2(overall, detail, obs_point_prefix=None, output_path=None, output_file=None, bc_path=None, bc_file=None, overall_crs=None, detail_crs=None, option=None, boundary_water_level_correction=0.0, filter_incoming=False, bctype="waterlevel", return_maximum=False)

   Read output from the overall model and set boundary conditions on the
   detail model.

   The correct implementation is selected automatically based on the class
   names of *overall* and *detail*.

   *overall* may be a string type name (``"sfincs"``, ``"hurrywave"``,
   ``"delft3dfm"``, ``"beware"``, ``"xbeach"``) instead of a model object.

   :param overall: The overall model, or a string type name.
   :type overall: Any or str
   :param detail: The detail model.
   :type detail: Any
   :param obs_point_prefix: Prefix for observation-point names.  Defaults to
       ``detail.name`` or ``"nest"``.
   :type obs_point_prefix: str or None
   :param output_path: Path to the overall model output directory.
   :type output_path: str or None
   :param output_file: Name of the output file (e.g. ``"sfincs_his.nc"``).
   :type output_file: str or None
   :param bc_path: If given, boundary-condition files are written to this
       directory.
   :type bc_path: str or None
   :param bc_file: Name of the boundary-condition file.
   :type bc_file: str or None
   :param overall_crs: CRS of the overall model (e.g. ``"EPSG:4326"``).
   :type overall_crs: str or None
   :param detail_crs: CRS of the detail model.
   :type detail_crs: str or None
   :param option: Implementation-specific option string.
   :type option: str or None
   :param boundary_water_level_correction: Constant offset added to all
       boundary water levels (metres).
   :type boundary_water_level_correction: float
   :param filter_incoming: Apply an incoming-wave filter to water levels.
   :type filter_incoming: bool
   :param bctype: Boundary-condition type (``"waterlevel"`` or ``"wave"``).
   :type bctype: str
   :param return_maximum: Return peak values instead of setting BCs.
   :type return_maximum: bool
   :returns: Implementation-specific result.  ``True`` on success, ``False``
       if the combination is not supported.
   :rtype: Any

----

Nest1 implementations
---------------------

Each function adds observation points to the *overall* model at the
boundary-point locations of the *detail* model.  They are called automatically
by ``nest1()``; direct use is rarely needed.

SFINCS as overall model
~~~~~~~~~~~~~~~~~~~~~~~

``nest1_sfincs_in_sfincs(overall, detail)``
   Adds observation points at SFINCS water-level boundary locations.

``nest1_xbeach_in_sfincs(overall, detail)``
   Adds observation points at XBeach flow-boundary locations.

``nest1_beware_in_sfincs(overall, detail)``
   Adds observation points at BEWARE flow-boundary locations.

HurryWave as overall model
~~~~~~~~~~~~~~~~~~~~~~~~~~

``nest1_hurrywave_in_hurrywave(overall, detail)``
   Adds spectral observation points at HurryWave boundary-condition locations.

``nest1_sfincs_in_hurrywave(overall, detail)``
   Adds observation points at SFINCS SnapWave boundary locations.

``nest1_xbeach_in_hurrywave(overall, detail)``
   Adds observation points at XBeach wave-boundary locations.

``nest1_beware_in_hurrywave(overall, detail)``
   Adds observation points at BEWARE wave-boundary locations.

Delft3D-FM as overall model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nest1_delft3dfm_in_delft3dfm(overall, detail)``
   Adds observation points at Delft3D-FM boundary locations.

``nest1_sfincs_in_delft3dfm(overall, detail)``
   Adds observation points at SFINCS water-level boundary locations.

``nest1_xbeach_in_delft3dfm(overall, detail)``
   Adds observation points at XBeach flow-boundary locations.

``nest1_beware_in_delft3dfm(overall, detail)``
   Adds observation points at BEWARE flow-boundary locations.

BEWARE as overall model
~~~~~~~~~~~~~~~~~~~~~~~~

BEWARE output points are fixed.  ``nest1`` returns ``True`` without adding
any observation points.

----

Nest2 implementations
---------------------

Each function reads output from the *overall* model and sets boundary
conditions on the *detail* model.  They are called automatically by
``nest2()``; direct use is rarely needed.

All functions accept the keyword arguments listed under ``nest2()`` above.

SFINCS as overall model
~~~~~~~~~~~~~~~~~~~~~~~

``nest2_sfincs_in_sfincs``
   Reads water-level time series from ``sfincs_his.nc`` and sets water-level
   BCs on the detail SFINCS model.  Supports ``filter_incoming`` (incoming-wave
   filter) and ``return_maximum`` (return peak values).

``nest2_xbeach_in_sfincs``
   Reads water-level time series and sets flow BCs on the detail XBeach model.

``nest2_beware_in_sfincs``
   Reads water-level time series and sets water-level BCs on the detail
   BEWARE model.

HurryWave as overall model
~~~~~~~~~~~~~~~~~~~~~~~~~~

``nest2_hurrywave_in_hurrywave``
   Reads spectral output and sets spectral BCs on the detail HurryWave model.

``nest2_sfincs_in_hurrywave``
   Reads wave-parameter time series (Hs, Tp, Wd, Ds) from
   ``hurrywave_his.nc`` and sets SnapWave BCs on the detail SFINCS model.

``nest2_xbeach_in_hurrywave``
   Reads wave-parameter time series and sets wave BCs on the detail XBeach
   model.  Supports deshoaling via the ``option`` parameter.

``nest2_beware_in_hurrywave``
   Reads wave-parameter time series and sets wave BCs on the detail BEWARE
   model.

Delft3D-FM as overall model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nest2_delft3dfm_in_delft3dfm``
   Reads water-level time series from Delft3D-FM output and sets water-level
   BCs on the detail Delft3D-FM model.

``nest2_sfincs_in_delft3dfm``
   Reads water-level time series via ``overall.read_timeseries_output()`` and
   sets water-level BCs on the detail SFINCS model.

``nest2_beware_in_delft3dfm``
   Reads water-level time series and sets water-level BCs on the detail
   BEWARE model.

BEWARE as overall model
~~~~~~~~~~~~~~~~~~~~~~~~

``nest2_sfincs_in_beware``
   Reads water level (WL) and wave setup (R2_setup) from BEWARE output and
   sets boundary conditions on the detail SFINCS model.
