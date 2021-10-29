import os
import argparse
from pathlib import Path
from parser import EmptyException, Parser

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', action='store_true')
    parser.add_argument('-f', help='Run solver on file')
    args = parser.parse_args()
    if args.f:
        parser = Parser(args.f)
        try:
            parser.parse()
            print("Running solver...")
            count = parser.problem.solve(log=args.v)
            print(f"Solution: {count}")
        except EmptyException:
            print("Could not find a problem :(")
    else:
        parser.print_help()
