from astroquery.gaia import Gaia
import os

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
os.makedirs("gaia", exist_ok=True)

queries = {
    "gaia_nearby_ms.csv": """
        SELECT TOP 5000
        source_id, ra, dec, parallax, parallax_error,
        pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe
        FROM gaiadr3.gaia_source
        WHERE parallax_over_error > 20
        AND parallax > 10
        AND phot_g_mean_mag < 14
        AND bp_rp IS NOT NULL
        AND ruwe < 1.4
    """,

    "gaia_hyades.csv": """
        SELECT TOP 20000
        source_id, ra, dec, parallax, parallax_error,
        pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', 66.75, 15.87, 8.0)
        )
        AND phot_g_mean_mag < 17
        AND bp_rp IS NOT NULL
        AND parallax > 10
        AND ruwe < 1.4
    """,

    "gaia_pleiades.csv": """
        SELECT TOP 20000
        source_id, ra, dec, parallax, parallax_error,
        pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', 56.75, 24.12, 3.0)
        )
        AND phot_g_mean_mag < 18
        AND bp_rp IS NOT NULL
        AND parallax > 4
        AND parallax < 10
        AND ruwe < 1.4
    """
}

for filename, query in queries.items():
    print(f"Downloading {filename}...")
    job = Gaia.launch_job(query)
    table = job.get_results()
    outpath = os.path.join("gaia", filename)
    table.write(outpath, format="csv", overwrite=True)
    print(f"Saved {outpath} with {len(table)} rows")

print("Done.")
