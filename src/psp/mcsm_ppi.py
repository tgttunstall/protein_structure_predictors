import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import urljoin
from io import StringIO

import pandas as pd

from .common import (
    INPUT_FIELDS,
    download_pdb,
    extract_job_url_and_id,
    get_session,
    load_mutations,
    log,
    save_binary,
    unique_path,
    write_rows,
)

SUBMIT_URL = "https://biosig.lab.uq.edu.au/mcsm/prot_prot_prediction_list"


def build_mutation_line(row: Dict[str, str]) -> str:
    return f"{row['chain'].strip()} {row['wt_aa'].strip()}{row['residue_number'].strip()}{row['mut_aa'].strip()}"


def submit(input_csv: str, output_csv: str, pdb_dir: str) -> None:
    rows = load_mutations(input_csv)
    by_pdb: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        by_pdb.setdefault(row["pdb_id"].strip(), []).append(row)

    session = get_session()
    records: List[Dict[str, str]] = []

    for idx, (pdb_id, mutations) in enumerate(by_pdb.items(), start=1):
        log(f"Submitting mCSM-PPI job {idx}/{len(by_pdb)} for PDB {pdb_id} ({len(mutations)} mutations)")
        pdb_path = download_pdb(pdb_id, pdb_dir)
        mutation_lines = "\n".join(build_mutation_line(m) for m in mutations)

        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile("w+", delete=False) as tmp:
            tmp.write(mutation_lines)
            tmp.flush()
            tmp_path = tmp.name

        try:
            with open(pdb_path, "rb") as pdb_file, open(tmp_path, "rb") as mutation_file:
                files = {"wild": pdb_file, "mutation_list": mutation_file}
                response = session.post(SUBMIT_URL, files=files, timeout=60)
                response.raise_for_status()
                job_url_abs, job_id = extract_job_url_and_id(
                    response.text,
                    SUBMIT_URL,
                    keywords=["/mcsm/output", "prot_prot", "results"],
                    path_hint="results_prot_prot",
                    response_url=response.url,
                )
        finally:
            Path(tmp_path).unlink(missing_ok=True)

        for mutation in mutations:
            record = {
                **{field: mutation[field] for field in INPUT_FIELDS},
                "job_id": job_id,
                "job_url": job_url_abs,
            }
            records.append(record)

    output_fields = INPUT_FIELDS + ["job_id", "job_url"]
    write_rows(output_csv, output_fields, records)


def fetch_jobs(output_csv: str, results_dir: str, force: bool = False) -> None:
    import csv

    with open(output_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "job_url" not in reader.fieldnames:
            raise ValueError("output CSV must contain job_url column")
        rows = list(reader)

    session = get_session()

    for row in rows:
        job_id = row.get("job_id") or Path(row["job_url"]).name
        target_html = Path(results_dir) / f"mcsm_ppi_{job_id}.html"
        if target_html.exists() and not force:
            continue
        log(f"Fetching mCSM-PPI job {job_id} (HTML)")
        response = session.get(row["job_url"], timeout=60)
        response.raise_for_status()
        path_html = (
            unique_path(results_dir, f"mcsm_ppi_{job_id}", ".html")
            if target_html.exists() and not force
            else target_html
        )
        save_binary(response.content, path_html)


def format_results(results_dir: str, formatted_csv: str) -> None:
    records: List[pd.DataFrame] = []
    paths = list(Path(results_dir).glob("mcsm_ppi_*.html"))
    for path in paths:
        html = path.read_text(encoding="utf-8", errors="ignore")
        try:
            tables = pd.read_html(StringIO(html))
        except ValueError:
            continue
        if not tables:
            continue
        df = tables[0]
        drop_cols = [
            col
            for col in df.columns
            if (
                (isinstance(col, str) and (col.strip().lower() == "index" or col.startswith("Unnamed")))
                or (not isinstance(col, str) and str(col).startswith("Unnamed"))
            )
        ]
        if drop_cols:
            df = df.drop(columns=drop_cols)
        df.insert(0, "source_file", str(path))
        records.append(df)
    if not records:
        raise RuntimeError("No mCSM-PPI result tables found to format")
    combined = pd.concat(records, ignore_index=True)
    combined.to_csv(formatted_csv, index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="mCSM-PPI wrapper")
    sub = parser.add_subparsers(dest="command", required=True)

    submit_p = sub.add_parser("submit", help="Submit mutation list to mCSM-PPI")
    submit_p.add_argument("--input", required=True, help="Input CSV path")
    submit_p.add_argument("--output", required=True, help="Output CSV for job info")
    submit_p.add_argument("--pdb-dir", default="tests/fixtures/pdb", help="Directory to cache PDB files")

    fetch_p = sub.add_parser("fetch", help="Fetch job results")
    fetch_p.add_argument("--jobs", required=True, help="Job CSV from submit step")
    fetch_p.add_argument("--results-dir", default="results/mcsm_ppi", help="Directory to store HTML results")
    fetch_p.add_argument("--force", action="store_true", help="Overwrite existing result files")

    format_p = sub.add_parser("format", help="Format downloaded results to CSV")
    format_p.add_argument("--results-dir", default="results/mcsm_ppi", help="Directory containing result HTML files")
    format_p.add_argument("--output", required=True, help="Formatted CSV output path")

    return parser


def main(argv: Optional[Iterable[str]] = None) -> None:
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
