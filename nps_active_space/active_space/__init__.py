"""Active space generation modules.

Heavy imports (GDAL, NMSim propagation) are lazy so AAM-only code paths and unit
tests can import ``aam_propagation_model`` without loading the full geospatial stack.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nps_active_space.active_space.active_space_generator import (
        ActiveSpaceGenerator,
        human_hearing_threshold,
    )
    from nps_active_space.active_space.layered_active_space import LayeredActiveSpace

__all__ = ["ActiveSpaceGenerator", "human_hearing_threshold", "LayeredActiveSpace"]


def __getattr__(name: str):
    if name in ("ActiveSpaceGenerator", "human_hearing_threshold"):
        from nps_active_space.active_space.active_space_generator import (
            ActiveSpaceGenerator,
            human_hearing_threshold,
        )
        if name == "ActiveSpaceGenerator":
            return ActiveSpaceGenerator
        return human_hearing_threshold
    if name == "LayeredActiveSpace":
        from nps_active_space.active_space.layered_active_space import LayeredActiveSpace
        return LayeredActiveSpace
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
