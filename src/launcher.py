import os
import argparse
from pathlib import Path
from parser import Parser
from formulas import PosFormula, InFormula

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Run solver on file')
    args = parser.parse_args()
    if args.f:
        parser = Parser(args.f)
        problem = parser.parsed
        print("Running solver...")
        count = problem.solve(log=True)
    else:
        parser.print_help()
