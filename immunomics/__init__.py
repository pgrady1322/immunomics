#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Package initialization and version metadata.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

__version__ = "0.1.0"

__all__ = ["analysis", "data", "integration", "visualization"]


def __getattr__(name: str):
    """Lazy-import subpackages to avoid circular-import errors in CI."""
    if name in __all__:
        import importlib

        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
