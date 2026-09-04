"""Active space generation modules.

``ActiveSpaceGenerator`` is imported lazily so NMSim installs do not load the
optional AAM adapter (``aam_translator`` / ``.[aam]`` extra).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
    from nps_active_space.active_space.audibility import human_hearing_threshold
    from nps_active_space.active_space.layered_active_space import LayeredActiveSpace

__all__ = ["ActiveSpaceGenerator", "human_hearing_threshold", "LayeredActiveSpace"]


def __getattr__(name: str):
    if name == "ActiveSpaceGenerator":
        from nps_active_space.active_space.active_space_generator import ActiveSpaceGenerator
        return ActiveSpaceGenerator
    if name == "human_hearing_threshold":
        from nps_active_space.active_space.audibility import human_hearing_threshold
        return human_hearing_threshold
    if name == "LayeredActiveSpace":
        from nps_active_space.active_space.layered_active_space import LayeredActiveSpace
        return LayeredActiveSpace
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
