"""suncat: SUNA + CTD processing toolkit.

Two layers, kept deliberately separate:
  - a PORTABLE science core (io, calc, qc, nitrate, productivity, sections, viz) that works on
    already-produced .cnv and L2/L3 products, installable anywhere;
  - an OPTIONAL, Windows-only SBE backend (suncat.sbe) that drives the proprietary
    SBEDataProcessing executables via the .psa files, used on the acquisition machine only.
"""
__version__ = "0.0.0"
