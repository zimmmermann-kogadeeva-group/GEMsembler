import argparse
import sys
import yaml

from .pathsfinding import pathsfinding


def pathsfinding_cli(args):
    parser = argparse.ArgumentParser("pathsfinding")
    parser.add_argument("path_to_model")
    parser.add_argument("path_to_config")
    parser.add_argument("output_dir")

    args = parser.parse_args(args)

    with open(args.path_to_config) as fh:
        config = yaml.safe_load(fh)

    pathsfinding(args.path_to_model, args.output_dir, **config)


def main():
    parser = argparse.ArgumentParser(prog="gemsembler")

    parser.add_argument("command", help="Subcommand to run: pathsfinding")
    # parse_args defaults to [1:] for args, but you need to
    # exclude the rest of the args too, or validation will fail
    args = parser.parse_args(sys.argv[1:2])

    if args.command == "pathsfinding":
        pathsfinding_cli(sys.argv[2:])


if __name__ == "__main__":
    main()
