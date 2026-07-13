"""Production continuous-preprocessing support for the DEMI EEG reanalysis.

The package implements only continuous preprocessing. It deliberately exposes
no epoch, AutoReject, CSD, event-repair, or participant-inclusion API.
"""

from .contracts import (
    ALL_RETAINED_CHANNELS,
    EEG_TARGET_CHANNELS,
    MASTOID_CHANNELS,
    SCALP_SOURCE_CHANNELS,
    load_config,
)

__all__ = [
    "ALL_RETAINED_CHANNELS",
    "EEG_TARGET_CHANNELS",
    "MASTOID_CHANNELS",
    "SCALP_SOURCE_CHANNELS",
    "load_config",
]
