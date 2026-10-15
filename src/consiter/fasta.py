"""Pure-Python helpers for FASTA headers, sequence dictionaries and coverage intervals."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path


def read_dict_names(dict_path: Path) -> list[str]:
    """Return sequence names, in order, from a SAM-style ``.dict`` file."""
    names = []
    with open(dict_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("@SQ"):
                continue
            for field in line.rstrip("\n").split("\t")[1:]:
                if field.startswith("SN:"):
                    names.append(field[3:])
                    break
    return names


def read_fasta_names(fasta_path: Path) -> list[str]:
    """Return the first word of every header in a FASTA file."""
    with open(fasta_path, encoding="utf-8") as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def rename_headers(src: Path, dest: Path, names: list[str], prefix: str | None = None) -> None:
    """Copy ``src`` to ``dest``, replacing the i-th header with ``names[i]``.

    If ``prefix`` is given, headers become ``prefix|name``.
    """
    with open(src, encoding="utf-8") as fin, open(dest, "w", encoding="utf-8") as fout:
        n = 0
        for line in fin:
            if line.startswith(">"):
                if n >= len(names):
                    raise ValueError(f"{src} has more records than the {len(names)} names given")
                name = names[n] if prefix is None else f"{prefix}|{names[n]}"
                line = f">{name}\n"
                n += 1
            fout.write(line)
    if n != len(names):
        raise ValueError(f"{src} has {n} records but {len(names)} names were given")


def low_coverage_intervals(bedgraph: Iterable[str], min_cov: int) -> Iterator[str]:
    """Yield BED lines for bedgraph intervals with coverage below ``min_cov``."""
    for line in bedgraph:
        chrom, start, end, depth = line.split("\t")[:4]
        if float(depth) < min_cov:
            yield f"{chrom}\t{start}\t{end}\n"
