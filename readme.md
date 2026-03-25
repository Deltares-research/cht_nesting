# Coastal Hazards Toolkit - Nesting

Two-step nesting workflow for coupling coastal hydrodynamic and wave models.
Reads output from a coarse (overall) model and applies it as boundary conditions
to a fine (detail) model.

## Supported model combinations

| Overall | Detail |
|---|---|
| SFINCS / SfincsModel | SFINCS / SfincsModel |
| SFINCS / SfincsModel | XBeach |
| SFINCS / SfincsModel | BEWARE |
| HurryWave | HurryWave |
| HurryWave | SFINCS / SfincsModel |
| HurryWave | XBeach |
| HurryWave | BEWARE |
| Delft3D-FM | Delft3D-FM |
| Delft3D-FM | SFINCS / SfincsModel |
| Delft3D-FM | XBeach |
| Delft3D-FM | BEWARE |
| BEWARE | SFINCS / SfincsModel |

## Usage

```python
from cht_nesting import nest1, nest2

# Step 1: add observation points to overall model at detail boundary locations
nest1(overall_model, detail_model)

# Run overall model ...

# Step 2: extract time series and set as boundary conditions on detail model
nest2(overall_model, detail_model)
```

## Installation

```bash
pip install -e .
```

Or directly from GitHub:

```bash
pip install git+https://github.com/deltares-research/cht_nesting.git
```
