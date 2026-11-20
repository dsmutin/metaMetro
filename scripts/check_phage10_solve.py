"""Check a MetaMetro phage_10 default build and a metaMalevich solve summary (CI job ``phage10-solve``).

    python scripts/check_phage10_solve.py identity BENCH_DIR
    python scripts/check_phage10_solve.py solve SUMMARY_JSON [SOLVE_WORKDIR]

``identity`` compares ``BENCH_DIR/identity.sha256`` with
``metametro.bench.canonical_identity("phage_10")`` and ``gold_sha256`` in
``BENCH_DIR/manifest.yaml`` with ``canonical_gold("phage_10")``, and checks the
``ground_truth/`` files against ``gold.sha256``. ``solve`` prints the
phage_10 row of the summary and requires ``build_status == "built"``,
``solved is True``, ``kraken2`` applied and no colouring skipped,
``kraken2_input == "metametro_pin"``, and non-empty ``metrics``. It checks the
identity of the build that solve made: ``identity_as_built`` must equal the
canonical identity when the row has it; otherwise ``SOLVE_WORKDIR/<bench>/
identity.sha256`` must equal ``SOLVE_IDENTITIES`` for the applied colourings
(the canonical identity for any other set), and its ``gold_sha256`` must equal
``canonical_gold``. It does not trust the solve exit code or the summary's
``ok`` flag (metaMalevich #275), only the row.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

NAMES = {"phage_10", "phage_species_5_x10"}


def check_identity(bench_dir: Path) -> int:
    from metametro.bench import canonical_gold, canonical_identity
    from metametro.bench.ground_truth import check_gold_digest, gold_digest
    from metametro.errors import ContractError
    from metametro.tables import read_yaml

    expected = canonical_identity("phage_10")
    found = (bench_dir / "identity.sha256").read_text(encoding="utf-8").split()[0]
    print(f"phage_10 identity {found}, canonical {expected}")
    gold = read_yaml(bench_dir / "manifest.yaml").get("gold_sha256")
    print(f"phage_10 gold_sha256 {gold}, canonical {canonical_gold('phage_10')}")
    try:
        check_gold_digest(bench_dir)
    except ContractError as exc:
        print(exc)
        return 1
    ok = found == expected and gold == canonical_gold("phage_10") == gold_digest(bench_dir)
    return 0 if ok else 1


# identity.sha256 of the phage_10 build that ``metamalevich solve`` makes, by the
# colourings it applied. A solve that builds with colourings=None must give the
# canonical identity; older metaMalevich solves build with SOLVE_COLOURINGS.
SOLVE_IDENTITIES = {
    ("kraken2", "composition_kmeans", "decaying"): "91ddcb0f98f1c57ce093f7b53a4667addd5a267abe591acce60945df8e1a5184",
}


def _solve_problems(row: dict, workdir: Path | None) -> list[str]:
    from metametro.bench import canonical_gold, canonical_identity
    from metametro.tables import read_yaml

    problems = []
    if row.get("build_status") != "built" or row.get("solved") is not True:
        problems.append(f"not solved: build_status={row.get('build_status')!r} solved={row.get('solved')!r}")
    if "kraken2" not in row.get("colourings", []):
        problems.append(f"kraken2 was not applied: colourings={row.get('colourings')!r}")
    if row.get("colourings_skipped"):
        problems.append(f"colourings were skipped: {row['colourings_skipped']!r}")
    if row.get("kraken2_input") != "metametro_pin":
        problems.append(f"the MetaMetro Kraken2 pin was not used: kraken2_input={row.get('kraken2_input')!r}")
    if not isinstance(row.get("metrics"), dict) or not row["metrics"]:
        problems.append(f"metrics are missing: {row.get('metrics')!r} ({row.get('metrics_reason', '')})")
    canonical = canonical_identity("phage_10")
    if "identity_as_built" in row:
        if row["identity_as_built"] != canonical:
            problems.append(f"identity_as_built {row['identity_as_built']} != canonical {canonical}")
        if row.get("identity_canonical", canonical) != canonical:
            problems.append(f"identity_canonical {row['identity_canonical']} != canonical {canonical}")
        return problems
    if workdir is None:
        problems.append("the row has no identity_as_built; pass the solve --workdir to check the solved build")
        return problems
    bench = workdir / str(row.get("canonical_name", "phage_species_5_x10"))
    found = (bench / "identity.sha256").read_text(encoding="utf-8").split()[0]
    expected = SOLVE_IDENTITIES.get(tuple(row.get("colourings", [])), canonical)
    print(f"solved build identity {found}, expected {expected}")
    if found != expected:
        problems.append(f"solved build identity {found} != expected {expected}")
    gold = read_yaml(bench / "manifest.yaml").get("gold_sha256")
    if gold != canonical_gold("phage_10"):
        problems.append(f"solved build gold_sha256 {gold} != canonical {canonical_gold('phage_10')}")
    return problems


def check_solve(summary_path: Path, workdir: Path | None = None) -> int:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = [
        row
        for row in summary.get("benches", [])
        if row.get("bench") in NAMES or row.get("canonical_name") in NAMES
    ]
    if len(rows) != 1:
        print(f"expected one phage_10 row in {summary_path}, found {len(rows)}")
        return 1
    row = rows[0]
    print(json.dumps(row, indent=2, sort_keys=True))
    problems = _solve_problems(row, workdir)
    for problem in problems:
        print(f"phage_10: {problem}")
    if not problems:
        print("phage_10 solve OK: built, solved, kraken2 from the pin, metrics present, identity and gold as expected")
    return 1 if problems else 0


def main(argv: list[str]) -> int:
    if len(argv) not in {3, 4} or argv[1] not in {"identity", "solve"} or (argv[1] == "identity" and len(argv) == 4):
        print(__doc__)
        return 2
    path = Path(argv[2])
    if argv[1] == "identity":
        return check_identity(path)
    return check_solve(path, Path(argv[3]) if len(argv) == 4 else None)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
