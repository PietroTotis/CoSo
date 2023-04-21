import os
import argparse

from src.launcher import run_coso, run_viscoso

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Run solver on CoLa file")
    parser.add_argument(
        "-v",
        nargs='?',
        const=os.path.join("src", "VisCoSo", "viscoso.html"),
        help="Generate a visual representation of CoSo reasoning to an html file",
    )
    parser.add_argument("--debug", action="store_true", help="Print log")
    args = parser.parse_args()
    if args.filename:
        print("Running solver...")
        if args.v:
            sol = run_viscoso(file=args.filename, visual=args.v)
        else:
            sol = run_coso(file=args.filename, debug=args.debug)
        print(f"Solution: {sol}")
    else:
        parser.print_help()
