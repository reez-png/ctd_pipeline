from __future__ import annotations
import argparse
from . import __version__


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="suncat", description="SUNA + CTD toolkit")
    p.add_argument("--version", action="version", version=f"suncat {__version__}")
    sub = p.add_subparsers(dest="command")
    for name, help_ in [("process", "run the processing chain for a cruise"),
                        ("qc", "run the QC/QA diagnostics for a cruise")]:
        sp = sub.add_parser(name, help=help_)
        sp.add_argument("cruise", help="cruise id, e.g. P45_06")
    args = p.parse_args(argv)
    if not args.command:
        p.print_help(); return 0
    print(f"suncat {args.command}: not implemented yet (scaffold). Cruise: {args.cruise}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
