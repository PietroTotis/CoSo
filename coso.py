import argparse

from src.cola_parser import Parser
from src.problem import EmptyException


def coso(debug=False, file=None, cola=None, visual=None):
    parser = Parser(file, cola)
    try:
        parser.parse()
        print("Running solver...")
        sol = parser.problem.solve(debug=debug)
        print(f"Solution: {sol}")
        if visual is not None:
            with open(visual, "w") as out:
                out.write(sol.log.to_viscoso())
    except EmptyException:
        print("Could not find a problem :(")


def viscoso(file=None, cola=None):
    parser = Parser(file, cola)
    try:
        parser.parse()
        sol = parser.problem.solve(debug=False)
        return sol.log.to_viscoso()
    except EmptyException:
        print("Could not find a problem :(")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Run solver on CoLa file")
    parser.add_argument(
        "-v",
        help="Generate a visual representation of CoSo reasoning to an html file",
    )
    parser.add_argument("--debug", action="store_true", help="Print log")
    args = parser.parse_args()
    if args.filename:
        if args.v:
            coso(file=args.filename, debug=args.debug, visual=args.v)
        else:
            coso(file=args.filename, debug=args.debug)
    else:
        parser.print_help()
