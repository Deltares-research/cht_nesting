"""Coastal Hazards Toolkit - Nesting.

Two-step nesting workflow for coupling coastal hydrodynamic and wave models.
"""

from .nest1.nest1 import nest1 as nest1
from .nest2.nest2 import nest2 as nest2

__version__ = "1.0.0"
