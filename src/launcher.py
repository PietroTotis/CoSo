import os
import argparse
from pathlib import Path
from parser import EmptyException, Parser

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Run solver on file")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--visualize", help="Generate VisCoSo log to given file")
    args = parser.parse_args()
    if args.f:
        parser = Parser(args.f)
        try:
            parser.parse()
            print("Running solver...")
            count = parser.problem.solve(debug=args.debug)
            print(f"Solution: {count}")
            if args.visualize:
                with open(args.visualize, "w") as out:
                    out.write(count.log.to_viscoso())
        except EmptyException:
            print("Could not find a problem :(")
    else:
        parser.print_help()
