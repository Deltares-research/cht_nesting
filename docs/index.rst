cht_nesting documentation
#########################

**cht_nesting** is a Python toolkit for two-step nesting of coastal hydrodynamic
and wave models.  It couples a coarse *overall* model with a fine *detail* model
by transferring boundary conditions through observation-point output.

Workflow
--------

1. **Nesting step 1** — ``nest1()`` adds observation points to the
   overall model at every boundary-point location of the detail model.
2. **Run the overall model** (outside cht_nesting).
3. **Nesting step 2** — ``nest2()`` reads the overall model output
   at those observation points and sets boundary conditions on the detail model.

Supported models
----------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Model
     - Package / class
   * - SFINCS
     - ``hydromt_sfincs.SfincsModel``
   * - HurryWave
     - ``hydromt_hurrywave.HurrywaveModel``
   * - Delft3D-FM
     - ``cht_delft3dfm.Delft3DFM``
   * - XBeach
     - ``cht_xbeach.XBeach``
   * - BEWARE
     - ``cht_beware.BEWARE``

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   user_guide/index
   api/index
   changelog
