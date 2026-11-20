"""Check a MetaMetro phage_10 default build and a metaMalevich solve summary (CI job ``phage10-solve``).

    python scripts/check_phage10_solve.py identity BENCH_DIR
    python scripts/check_phage10_solve.py solve SUMMARY_JSON

``identity`` compares ``BENCH_DIR/identity.sha256`` with
``metametro.bench.canonical_identity("phage_10")``. ``solve`` prints the
phage_10 row of the summary and requires ``build_status == "built"`` and
``solved is True``. It does not trust the solve exit code or the summary's
``ok`` flag (metaMalevich #275), only the row.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

NAMES = {"phage_10", "phage_species_5_x10"}


def check_identity(bench_dir: Path) -> int:
    from metametro.bench import canonical_identity

    expected = canonical_identity("phage_10")
    found = (bench_dir / "identity.sha256").read_text(encoding="utf-8").split()[0]
    print(f"phage_10 identity {found}, canonical {expected}")
    return 0 if found == expected else 1


def check_solve(summary_path: Path) -> int:
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
    if row.get("build_status") != "built" or row.get("solved") is not True:
        print(f"phage_10 not solved: build_status={row.get('build_status')!r} solved={row.get('solved')!r}")
        return 1
    return 0


def main(argv: list[str]) -> int:
    if len(argv) != 3 or argv[1] not in {"identity", "solve"}:
        print(__doc__)
        return 2
    path = Path(argv[2])
    return check_identity(path) if argv[1] == "identity" else check_solve(path)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
