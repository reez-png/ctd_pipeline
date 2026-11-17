"""pytest configuration for the ctd_pipeline project.

The mere presence of this file at the project root fixes imports: pytest adds the
directory that holds the root conftest.py to sys.path, so tests can `import ctd_lib`
(and any other top-level module) regardless of how pytest is launched, bare `pytest`,
`python -m pytest`, PyCharm's run arrow, or CI. Without it, bare `pytest` only puts the
tests/ folder on the path and the import fails.

Otherwise it only pins BLAS to one thread (below).
"""

# One BLAS thread for the test run. The lag-model tests solve many small least-squares systems; a multi-threaded
# BLAS spends longer coordinating threads than computing on matrices this size, badly so on many-core Windows
# machines. Set here, before any test imports numpy; a value the user has set is left alone.
import os as _os

for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")
