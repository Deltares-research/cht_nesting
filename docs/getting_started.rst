Getting started
===============

Installation
------------

Install from PyPI:

.. code-block:: bash

   pip install cht_nesting

For development:

.. code-block:: bash

   git clone https://github.com/Deltares-research/cht_nesting.git
   cd cht_nesting
   pip install -e ".[tests]"

Quick example
-------------

The typical workflow has two steps: first add observation points to the overall
model, then — after the overall model has been run — read its output and set
boundary conditions on the detail model.

.. code-block:: python

   from cht_nesting import nest1, nest2

   # Step 1 — add observation points
   nest1(overall_model, detail_model, obs_point_prefix="nested")

   # ... run the overall model ...

   # Step 2 — transfer boundary conditions
   nest2(
       overall=overall_model,
       detail=detail_model,
       obs_point_prefix="nested",
       output_path="/path/to/overall/output",
       boundary_water_level_correction=0.1,
   )

The ``overall`` argument of ``nest2`` can also be a string instead of a model
object.  This is useful when you only have the output files and do not need to
instantiate the full model:

.. code-block:: python

   nest2(
       overall="sfincs",
       detail=detail_model,
       obs_point_prefix="nested",
       output_path="/path/to/sfincs/output",
   )

Accepted string values: ``"sfincs"``, ``"hurrywave"``, ``"delft3dfm"``,
``"beware"``, ``"xbeach"``.
