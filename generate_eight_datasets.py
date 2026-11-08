"""User-facing entry point for the eight-source dataset workflow.

The implementation remains in :mod:`generate_all_datasets` so the established
generic entry point and imports remain backward compatible.
"""
from __future__ import annotations

from generate_all_datasets import *  # noqa: F401,F403


if __name__ == "__main__":
    raise SystemExit(main())
