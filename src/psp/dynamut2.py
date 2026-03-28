import argparse
from io import StringIO
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from .common import (
    INPUT_FIELDS,
    extract_job_url_and_id,
    get_session,
    load_mutations,
    log,
    save_binary,
    unique_path,
    write_rows,
)

SUBMIT_URL = "https://biosig.lab.uq.edu.au/dynamut2/run_list_prediction"


def build_mutation_line(row: Dict[str, str]) -> str:
    return f"{row['chain'].strip()} {row['wt_aa'].strip()}{row['residue_number'].strip()}{row['mut_aa'].strip()}"


def submit(input_csv: str, output_csv: str, pdb_dir: str) -> None:
    rows = load_mutations(input_csv)
    by_pdb: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        by_pdb.setdefault(row["pdb_id"].strip(), []).append(row)

    session = get_session()
    records: List[Dict[str, str]] = []

    from tempfile import NamedTemporaryFile

    for idx, (pdb_id, mutations) in enumerate(by_pdb.items(), start=1):
        log(f"Submitting DynaMut2 job {idx}/{len(by_pdb)} for PDB {pdb_id} ({len(mutations)} mutations)")
        mutation_lines = "\n".join(build_mutation_line(m) for m in mutations)

        with NamedTemporaryFile("w+", delete=False) as tmp:
            tmp.write(mutation_lines)
            tmp.flush()
            tmp_path = tmp.name

        try:
            with open(tmp_path, "rb") as mutation_file:
                files = {"mutation_list": ("mutation_list.txt", mutation_file, "text/plain")}
                data = {"pdb_accession_list": pdb_id}
                response = session.post(SUBMIT_URL, data=data, files=files, timeout=60)
                response.raise_for_status()
                job_url_abs, job_id = extract_job_url_and_id(
                    response.text,
                    SUBMIT_URL,
                    keywords=["dynamut2", "results_prediction"],
                    path_hint="results_prediction",
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

    write_rows(output_csv, INPUT_FIELDS + ["job_id", "job_url"], records)


def fetch_jobs(
    output_csv: str, results_dir: str, force: bool = False, max_attempts: int = 6, wait_seconds: int = 30
) -> None:
    import csv
    import time
    from bs4 import BeautifulSoup

    with open(output_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "job_url" not in reader.fieldnames:
            raise ValueError("output CSV must contain job_url column")
        rows = list(reader)

    session = get_session()

    seen: set[str] = set()

    for row in rows:
        job_id = row.get("job_id") or Path(row["job_url"]).name
        if job_id in seen and not force:
            continue
        seen.add(job_id)
        target = Path(results_dir) / f"dynamut2_{job_id}.html"
        if target.exists() and not force:
            continue

        content = None
        url = row["job_url"]

        for attempt in range(1, max_attempts + 1):
            log(f"Fetching DynaMut2 job {job_id} (attempt {attempt}/{max_attempts})")
            response = session.get(url, timeout=60)
            response.raise_for_status()
            html = response.text
            lower = html.lower()

            soup = BeautifulSoup(html, "html.parser")
            # Follow meta refresh if present
            refresh = soup.find("meta", attrs={"http-equiv": "refresh"})
            if refresh and "content" in refresh.attrs:
                parts = refresh["content"].split(";")
                for part in parts:
                    if "url=" in part.lower():
                        next_url = part.split("=", 1)[1].strip()
                        url = urljoin(response.url, next_url)
                        break

            processing = ("preloader-wrapper" in lower) or ("processing" in lower)
            if processing and attempt < max_attempts:
                time.sleep(wait_seconds)
                continue

            content = response.content
            break

        if content is None:
            log(f"Failed to retrieve completed results for job {job_id}")
            continue

        path = (
            unique_path(results_dir, f"dynamut2_{job_id}", ".html")
            if target.exists() and not force
            else target
        )
        save_binary(content, path)


def format_results(results_dir: str, formatted_csv: str) -> None:
    frames: List[pd.DataFrame] = []
    for path in Path(results_dir).glob("dynamut2_*.html"):
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
                (isinstance(col, str) and (col.strip().lower() in {"index", "#"} or col.startswith("Unnamed")))
                or (not isinstance(col, str) and str(col).startswith("Unnamed"))
            )
        ]
        if drop_cols:
            df = df.drop(columns=drop_cols)
        df.insert(0, "source_file", str(path))
        frames.append(df)
    if not frames:
        raise RuntimeError("No DynaMut2 result tables found to format")
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(formatted_csv, index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DynaMut2 wrapper")
    sub = parser.add_subparsers(dest="command", required=True)

    submit_p = sub.add_parser("submit", help="Submit mutation lists to DynaMut2")
    submit_p.add_argument("--input", required=True, help="Input CSV path")
    submit_p.add_argument("--output", required=True, help="Output CSV for job info")
    submit_p.add_argument("--pdb-dir", default="tests/fixtures/pdb", help="(Unused) kept for interface consistency")

    fetch_p = sub.add_parser("fetch", help="Fetch job results")
    fetch_p.add_argument("--jobs", required=True, help="Job CSV from submit step")
    fetch_p.add_argument("--results-dir", default="results/dynamut2", help="Directory to store HTML results")
    fetch_p.add_argument("--force", action="store_true", help="Overwrite existing result files")
    fetch_p.add_argument("--max-attempts", type=int, default=6, help="Max attempts while job processes")
    fetch_p.add_argument("--wait-seconds", type=int, default=30, help="Seconds to wait between attempts")

    format_p = sub.add_parser("format", help="Format downloaded results to CSV")
    format_p.add_argument("--results-dir", default="results/dynamut2", help="Directory containing result HTML files")
    format_p.add_argument("--output", required=True, help="Formatted CSV output path")

    return parser


def main(argv: Optional[Iterable[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command == "submit":
        submit(args.input, args.output, args.pdb_dir)
    elif args.command == "fetch":
        fetch_jobs(
            args.jobs,
            args.results_dir,
            force=args.force,
            max_attempts=args.max_attempts,
            wait_seconds=args.wait_seconds,
        )
    elif args.command == "format":
        format_results(args.results_dir, args.output)
    else:
        parser.error("Unknown command")


if __name__ == "__main__":
    main()
