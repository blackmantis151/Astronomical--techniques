from astroquery.vizier import Vizier
import os

os.makedirs("tf", exist_ok=True)

Vizier.ROW_LIMIT = -1

catalog = "J/ApJ/902/145"

tables = Vizier.get_catalogs(catalog)

for name, table in zip(tables.keys(), tables.values()):
    print(name, len(table), table.colnames)
    safe_name = name.replace("/", "_")
    table.write(f"tf/{safe_name}.csv", format="csv", overwrite=True)

print("Done. Saved Cosmicflows-4 Tully-Fisher tables in tf/")

