from pathlib import Path
from astroquery.vizier import Vizier

ROOT = Path(__file__).resolve().parents[2]
outdir = ROOT / "ngc1866"
outdir.mkdir(exist_ok=True)

Vizier.ROW_LIMIT = -1

catalog = "J/A+A/440/487"

tables = Vizier.get_catalogs(catalog)

for name, table in zip(tables.keys(), tables.values()):
    print(name, len(table), table.colnames)
    safe = name.replace("/", "_")
    table.write(outdir / f"{safe}.csv", format="csv", overwrite=True)

print("Saved NGC 1866 Cepheid tables in ngc1866/")
