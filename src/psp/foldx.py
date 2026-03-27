"""FoldX integration (WIP).

This is a placeholder for migrating FoldX5 workflows from the legacy
`phd_generic_analysis` repo. FoldX is a local binary; no web submission is
involved. When wiring this up, the submit step should:
 - ensure FoldX binary is available (e.g., `foldx` or `foldx5` on PATH),
 - run RepairPDB / BuildModel for the provided mutations,
 - write job/result metadata similarly to other tools.

Currently, all commands raise a RuntimeError to signal that migration is
pending.
"""

import argparse
from typing import Iterable, Optional


def _not_implemented() -> None:
    raise RuntimeError(
        "FoldX integration is WIP. Install FoldX5 locally and migrate workflows from the legacy repo."
    )


def submit(input_csv: str, output_csv: str, pdb_dir: str) -> None:  # pragma: no cover - placeholder
    _not_implemented()


def fetch_jobs(output_csv: str, results_dir: str, force: bool = False) -> None:  # pragma: no cover - placeholder
    _not_implemented()


def format_results(results_dir: str, formatted_csv: str) -> None:  # pragma: no cover - placeholder
    _not_implemented()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="FoldX wrapper (WIP)")
    sub = parser.add_subparsers(dest="command", required=True)

    submit_p = sub.add_parser("submit", help="Submit mutations to FoldX (WIP)")
    submit_p.add_argument("--input", required=True, help="Input CSV path")
    submit_p.add_argument("--output", required=True, help="Output CSV for job info")
    submit_p.add_argument("--pdb-dir", default="input/pdb", help="Directory to cache PDB files")

    fetch_p = sub.add_parser("fetch", help="Fetch job results (WIP)")
    fetch_p.add_argument("--jobs", required=True, help="Job CSV from submit step")
    fetch_p.add_argument("--results-dir", default="results/foldx", help="Directory to store results")
    fetch_p.add_argument("--force", action="store_true", help="Overwrite existing result files")

    format_p = sub.add_parser("format", help="Format downloaded results (WIP)")
    format_p.add_argument("--results-dir", default="results/foldx", help="Directory containing result files")
    format_p.add_argument("--output", required=True, help="Formatted CSV output path")

    return parser


def main(argv: Optional[Iterable[str]] = None) -> None:  # pragma: no cover - placeholder
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command == "submit":
        submit(args.input, args.output, args.pdb_dir)
    elif args.command == "fetch":
        fetch_jobs(args.jobs, args.results_dir, force=args.force)
    elif args.command == "format":
        format_results(args.results_dir, args.output)
    else:
        parser.error("Unknown command")


if __name__ == "__main__":
    main()
