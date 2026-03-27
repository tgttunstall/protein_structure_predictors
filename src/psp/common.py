import csv
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import requests
from requests.adapters import HTTPAdapter, Retry

INPUT_FIELDS = [
    "mutation_id",
    "gene",
    "pdb_id",
    "chain",
    "residue_number",
    "wt_aa",
    "mut_aa",
]


def load_mutations(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [field for field in INPUT_FIELDS if field not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")
        rows = [row for row in reader]
    return rows


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def get_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": "protein-structure-predictors/0.1"})
    return session


def download_pdb(pdb_id: str, pdb_dir: str) -> str:
    ensure_dir(pdb_dir)
    pdb_id_clean = pdb_id.strip().upper()
    path = Path(pdb_dir) / f"{pdb_id_clean}.pdb"
    if path.exists():
        return str(path)
    url = f"https://files.rcsb.org/download/{pdb_id_clean}.pdb"
    session = get_session()
    response = session.get(url, timeout=30)
    response.raise_for_status()
    path.write_bytes(response.content)
    return str(path)


def write_rows(path: str, fieldnames: Iterable[str], rows: Iterable[Dict[str, str]]) -> None:
    path_obj = Path(path)
    ensure_dir(str(path_obj.parent))
    exists = path_obj.exists()
    with open(path_obj, "a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def unique_path(directory: str, stem: str, suffix: str) -> Path:
    ensure_dir(directory)
    base = Path(directory) / f"{stem}{suffix}"
    if not base.exists():
        return base
    counter = 1
    while True:
        candidate = Path(directory) / f"{stem}_{counter}{suffix}"
        if not candidate.exists():
            return candidate
        counter += 1


def find_first_link(html: str, keywords: Optional[List[str]] = None) -> Optional[str]:
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(html, "html.parser")
    anchors = soup.find_all("a", href=True)
    if not keywords:
        return anchors[0].get("href") if anchors else None
    for anchor in anchors:
        href = anchor.get("href")
        if any(keyword in href for keyword in keywords):
            return href
    return anchors[0].get("href") if anchors else None


def save_binary(content: bytes, path: Path) -> str:
    ensure_dir(str(path.parent))
    path.write_bytes(content)
    return str(path)


def extract_job_url_and_id(
    html: str,
    base_url: str,
    keywords: Optional[List[str]] = None,
    path_hint: Optional[str] = None,
) -> Tuple[str, str]:
    def _to_abs(href: str) -> str:
        return requests.compat.urljoin(base_url, href)

    # 1) Prefer anchor tag links filtered by keywords/path_hint
    link = find_first_link(html, keywords=keywords)
    if link:
        abs_url = _to_abs(link)
        job_id = Path(requests.utils.urlparse(abs_url).path).name
        if job_id:
            return abs_url, job_id

    # 2) Regex search for URLs in the HTML
    url_regex = re.compile(r"https?://[^\s\"']+", re.IGNORECASE)
    candidates = url_regex.findall(html)
    if path_hint:
        candidates = [c for c in candidates if path_hint in c]
    if keywords and candidates:
        candidates = [c for c in candidates if any(k in c for k in keywords)] or candidates
    if candidates:
        abs_url = candidates[0]
        job_id = Path(requests.utils.urlparse(abs_url).path).name
        if job_id:
            return abs_url, job_id

    # 3) Regex search for likely path with a job identifier
    path_regex = re.compile(r"/(?:[A-Za-z0-9_\-]+/)*([A-Za-z0-9_.-]+)")
    path_candidates = []
    if path_hint:
        hint_regex = re.compile(rf"/(?:[A-Za-z0-9_\-]+/)*{re.escape(path_hint)}[^\s\"']*", re.IGNORECASE)
        path_candidates = hint_regex.findall(html)
    if not path_candidates:
        path_candidates = path_regex.findall(html)
    if path_candidates:
        candidate = path_candidates[0]
        abs_url = _to_abs(candidate)
        job_id = Path(requests.utils.urlparse(abs_url).path).name
        if job_id:
            return abs_url, job_id

    # 4) Fallback: try to scrape a Job ID text pattern
    match = re.search(r"Job ID[:\s]*([A-Za-z0-9_.-]+)", html, flags=re.IGNORECASE)
    if match:
        job_id = match.group(1)
        abs_url = _to_abs(job_id)
        return abs_url, job_id

    raise RuntimeError("Could not locate job URL/ID in response")


def log(msg: str) -> None:
    print(f"[psp] {msg}")
