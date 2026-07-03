"""Repo-root pytest bootstrap.

Tests live in tests/ but import the model (palustris_me) from the repo root, so make the
root importable and register the `slow` marker. This lets a bare `pytest` run at the repo
root discover and run the whole suite (the config in tests/pytest.ini is used when pytest is
invoked from inside tests/ instead). No test logic here.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: LP-solve-heavy tests (~1.4 s per solve); deselect with -m 'not slow'"
    )
