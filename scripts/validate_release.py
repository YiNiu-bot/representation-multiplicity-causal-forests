#!/usr/bin/env python3
"""Validate public repository structure, metadata, and privacy boundaries."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "replication" / "metadata" / "release_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pdf_pages(path: Path) -> int | None:
    try:
        output = subprocess.check_output(
            ["pdfinfo", str(path)], text=True, stderr=subprocess.STDOUT
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    match = re.search(r"^Pages:\s+(\d+)$", output, flags=re.MULTILINE)
    return int(match.group(1)) if match else None


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for artifact in manifest["artifacts"]:
        path = ROOT / artifact["path"]
        if not path.is_file():
            raise SystemExit(f"Missing manifest artifact: {artifact['path']}")
        observed = sha256(path)
        if observed != artifact["sha256"]:
            raise SystemExit(f"Hash mismatch: {artifact['path']}")

    paper = ROOT / manifest["paper"]["path"]
    pages = pdf_pages(paper)
    if pages is not None and pages != manifest["paper"]["pages"]:
        raise SystemExit(f"PDF page mismatch: expected {manifest['paper']['pages']}, found {pages}")
    if pages is not None and pages > 45:
        raise SystemExit(f"PDF exceeds 45 pages: {pages}")

    forbidden_suffixes = {
        ".dta", ".rdata", ".rds", ".sav", ".sqlite", ".xlsx", ".zip",
        ".bst", ".cfg", ".cls",
    }
    forbidden_name_fragments = {
        "release_workflow": "internal release workflow",
        "finalize_release": "internal release finalizer",
        "referee": "referee material",
        "submission": "journal-submission material",
        "temporary-chat": "temporary-chat material",
    }
    forbidden_patterns = {
        "/Users/": "absolute macOS home path",
        "Submitted to Econometrica": "journal-submission header",
        "\\bibliographystyle{ecta}": "journal bibliography style",
        "\\documentclass[ecta": "journal document class",
        "RI_ENFORCE_REFERENCE_HARDWARE": "private hardware enforcement",
        "Hardware UUID": "hardware identifier",
        "Provisioning UDID": "hardware identifier",
        "Serial Number (system)": "hardware identifier",
    }
    forbidden_regexes = {
        re.compile(r"\bghp_[A-Za-z0-9]{30,}\b"): "GitHub personal access token",
        re.compile(r"\bgithub_pat_[A-Za-z0-9_]{40,}\b"): "GitHub fine-grained token",
        re.compile(r"\bAKIA[0-9A-Z]{16}\b"): "AWS access key",
        re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{20,}\b"): "Slack token",
        re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"): "private key",
    }
    text_suffixes = {
        ".R", ".bib", ".cff", ".csv", ".json", ".js", ".md",
        ".py", ".sh", ".tex", ".txt", "",
    }
    excluded_parts = {".git", ".r-library", "rerun_outputs", "__pycache__"}
    for path in ROOT.rglob("*"):
        if excluded_parts.intersection(path.parts) or not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.is_symlink():
            raise SystemExit(f"Symlink not allowed in release: {relative}")
        lowered_name = relative.as_posix().lower()
        for fragment, label in forbidden_name_fragments.items():
            if fragment in lowered_name:
                raise SystemExit(f"Forbidden {label}: {relative}")
        if path.suffix.lower() in forbidden_suffixes:
            raise SystemExit(f"Forbidden repository file: {relative}")
        if path.stat().st_size > 10 * 1024 * 1024:
            raise SystemExit(f"Release file exceeds 10 MiB: {relative}")
        if path.suffix in text_suffixes and path.name != "SHA256SUMS":
            text = path.read_text(encoding="utf-8", errors="strict")
            for pattern, label in forbidden_patterns.items():
                if pattern in text:
                    raise SystemExit(
                        f"Forbidden {label} in {relative}: {pattern}"
                    )
            for pattern, label in forbidden_regexes.items():
                if pattern.search(text):
                    raise SystemExit(f"Forbidden {label} in {relative}")

    if shutil.which("pdftotext"):
        text = subprocess.check_output(
            ["pdftotext", str(paper), "-"], text=True, stderr=subprocess.STDOUT
        )
        if "Submitted to Econometrica" in text:
            raise SystemExit("Journal-submission header found in public PDF")

    if any(path.name.startswith(".") and path.name == ".DS_Store" for path in ROOT.rglob("*")):
        raise SystemExit(".DS_Store found in release")

    print("PASS: repository manifest, page limit, paths, and privacy boundary.")


if __name__ == "__main__":
    main()
