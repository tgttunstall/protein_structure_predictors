import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List
from urllib.parse import urljoin

import pandas as pd

from .common import (
    INPUT_FIELDS,
    download_pdb,
    find_first_link,
    get_session,
    load_mutations,
    save_binary,
    unique_path,
    write_rows,
)

SUBMIT_URL = "https://biosig.lab.uq.edu.au/mcsm_lig/prediction"
REQUIRED_FIELDS = INPUT_FIELDS + ["ligand_id", "wt_affinity"]


def _load_rows(path: str) -> List[Dict[str, str]]:
    rows = load_mutations(path)
    missing = [field for field in ["ligand_id", "wt_affinity"] if field not in rows[0]] if rows else []
    if missing:
        raise ValueError(f"Missing required columns for mCSM-LIG: {', '.join(missing)}")
    return rows


def submit(input_csv: str, output_csv: str, pdb_dir: str) -> None:
    rows = _load_rows(input_csv)
    session = get_session()
    records: List[Dict[str, str]] = []

    for row in rows:
        pdb_path = download_pdb(row["pdb_id"], pdb_dir)
        mutation = f"{row['wt_aa'].strip()}{row['residue_number'].strip()}{row['mut_aa'].strip()}"
        data = {
            "mutation": mutation,
            "chain": row["chain"].strip(),
            "lig_id": row["ligand_id"].strip(),
            "affin_wt": row["wt_affinity"].strip(),
            "run": "single",
        }

        with open(pdb_path, "rb") as pdb_file:
            files = {"wild": pdb_file}
            response = session.post(SUBMIT_URL, data=data, files=files, timeout=60)
            response.raise_for_status()
        job_url = find_first_link(
            response.text,
            keywords=["mcsm_lig/output", "mcsm_lig", "prediction"],
        )
        if not job_url:
            raise RuntimeError("Could not locate job URL in mCSM-LIG response")
        job_url_abs = urljoin(SUBMIT_URL, job_url)
        job_id = job_url_abs.rstrip("/").split("/")[-1]

        record = {**{field: row[field] for field in REQUIRED_FIELDS}, "job_id": job_id, "job_url": job_url_abs}
        records.append(record)

    write_rows(output_csv, REQUIRED_FIELDS + ["job_id", "job_url"], records)


def fetch_jobs(output_csv: str, results_dir: str, force: bool = False) -> None:
    with open(output_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "job_url" not in reader.fieldnames:
            raise ValueError("output CSV must contain job_url column")
        rows = list(reader)

    session = get_session()

    for row in rows:
        job_id = row.get("job_id") or Path(row["job_url"]).name
        target = Path(results_dir) / f"mcsm_lig_{job_id}.html"
        if target.exists() and not force:
            continue
        response = session.get(row["job_url"], timeout=60)
        response.raise_for_status()
        path = unique_path(results_dir, f"mcsm_lig_{job_id}", ".html") if target.exists() and not force else target
        save_binary(response.content, path)


def format_results(results_dir: str, formatted_csv: str) -> None:
    frames: List[pd.DataFrame] = []
    for path in Path(results_dir).glob("mcsm_lig_*.html"):
        html = path.read_text(encoding="utf-8", errors="ignore")
        try:
            tables = pd.read_html(html)
        except ValueError:
            continue
        if not tables:
            continue
        df = tables[0]
        df.insert(0, "source_file", str(path))
        frames.append(df)
    if not frames:
        raise RuntimeError("No mCSM-LIG result tables found to format")
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(formatted_csv, index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="mCSM-LIG wrapper")
    sub = parser.add_subparsers(dest="command", required=True)

    submit_p = sub.add_parser("submit", help="Submit single mutations to mCSM-LIG")
    submit_p.add_argument("--input", required=True, help="Input CSV path (must include ligand_id, wt_affinity)")
    submit_p.add_argument("--output", required=True, help="Output CSV for job info")
    submit_p.add_argument("--pdb-dir", default="input/pdb", help="Directory to cache PDB files")

    fetch_p = sub.add_parser("fetch", help="Fetch job results")
    fetch_p.add_argument("--jobs", required=True, help="Job CSV from submit step")
    fetch_p.add_argument("--results-dir", default="results/mcsm_lig", help="Directory to store HTML results")
    fetch_p.add_argument("--force", action="store_true", help="Overwrite existing result files")

    format_p = sub.add_parser("format", help="Format downloaded results to CSV")
    format_p.add_argument("--results-dir", default="results/mcsm_lig", help="Directory containing result HTML files")
    format_p.add_argument("--output", required=True, help="Formatted CSV output path")

    return parser


def main(argv: Iterable[str] = None) -> None:
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
