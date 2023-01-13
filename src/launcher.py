import os
import argparse
from pathlib import Path
from parser import EmptyException, Parser

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Run solver on file F")
    parser.add_argument(
        "-v",
        help="Generate a visual representation of CoSo reasoning to an html file V",
    )
    parser.add_argument("--debug", action="store_true", help="Print log")
    args = parser.parse_args()
    if args.f:
        parser = Parser(args.f)
        try:
            parser.parse()
            print("Running solver...")
            sol = parser.problem.solve(debug=args.debug)
            print(f"Solution: {sol}")
            if args.v:
                with open(args.v, "w") as out:
                    out.write(sol.log.to_viscoso())
        except EmptyException:
            print("Could not find a problem :(")
    else:
        parser.print_help()
