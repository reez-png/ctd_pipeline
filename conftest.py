"""pytest configuration for the ctd_pipeline project.

The mere presence of this file at the project root fixes imports: pytest adds the
directory that holds the root conftest.py to sys.path, so tests can `import ctd_lib`
(and any other top-level module) regardless of how pytest is launched, bare `pytest`,
`python -m pytest`, PyCharm's run arrow, or CI. Without it, bare `pytest` only puts the
tests/ folder on the path and the import fails.

Intentionally empty otherwise.
"""
