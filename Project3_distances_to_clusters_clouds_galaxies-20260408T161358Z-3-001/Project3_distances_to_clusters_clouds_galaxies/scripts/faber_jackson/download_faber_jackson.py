from astroquery.vizier import Vizier
import os

os.makedirs("fj", exist_ok=True)
Vizier.ROW_LIMIT = -1

# 6dFGS Fundamental Plane / velocity dispersion catalog
catalogs_to_try = [
    "J/MNRAS/443/1231",
    "J/MNRAS/443/1231/catalog",
    "J/MNRAS/443/1231/table1",
]

for catalog in catalogs_to_try:
    print(f"\nTrying catalog: {catalog}")
    try:
        tables = Vizier.get_catalogs(catalog)
        print("Found tables:", tables.keys())

        for name, table in zip(tables.keys(), tables.values()):
            print(name, len(table), table.colnames)
            safe_name = name.replace("/", "_")
            table.write(f"fj/{safe_name}.csv", format="csv", overwrite=True)

        print("Success.")
        break

    except Exception as e:
        print("Failed:", e)

print("Done.")

