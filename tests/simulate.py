"""Small deterministic simulator for a mutated genome and paired-end reads."""

from __future__ import annotations

import random
from pathlib import Path

BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def mutate(rng: random.Random, seq: str, n_snps: int, indels: list[tuple[int, str | int]]) -> str:
    """Apply ``n_snps`` random SNPs plus indels ``(pos, inserted_seq | deleted_length)``.

    SNPs stay 400 bp away from the ends (where simulated coverage is low and
    strand-biased) and 30 bp away from each other, so each one is callable.
    """
    s = list(seq)
    taken = {pos for pos, _ in indels}
    snps = []
    while len(snps) < n_snps:
        pos = rng.randrange(400, len(seq) - 400)
        if all(abs(pos - t) > 30 for t in taken):
            taken.add(pos)
            snps.append(pos)
    for pos in snps:
        s[pos] = rng.choice([b for b in BASES if b != s[pos]])
    for pos, change in sorted(indels, reverse=True):
        if isinstance(change, str):
            s[pos + 1 : pos + 1] = list(change)
        else:
            del s[pos + 1 : pos + 1 + change]
    return "".join(s)


def write_fasta(path: Path, records: dict[str, str], width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name} simulated\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def simulate_reads(
    rng: random.Random,
    genome: dict[str, str],
    left: Path,
    right: Path,
    coverage: int = 100,
    read_len: int = 100,
    insert: int = 300,
) -> None:
    """Write error-free paired reads (FR orientation) uniformly over ``genome``."""
    with open(left, "w") as l_fh, open(right, "w") as r_fh:
        n = 0
        for seq in genome.values():
            for _ in range(len(seq) * coverage // (2 * read_len)):
                start = rng.randrange(0, len(seq) - insert + 1)
                frag = seq[start : start + insert]
                r1 = frag[:read_len]
                r2 = frag[-read_len:].translate(_COMPLEMENT)[::-1]
                qual = "I" * read_len
                l_fh.write(f"@r{n}/1\n{r1}\n+\n{qual}\n")
                r_fh.write(f"@r{n}/2\n{r2}\n+\n{qual}\n")
                n += 1
