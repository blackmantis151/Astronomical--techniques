from config import SCIENCE_DIR, LOG_DIR, OUTPUT_SUBDIRS
from utils import ensure_directories, build_science_dataframe


def group_science_files():
    ensure_directories(OUTPUT_SUBDIRS)

    science_files = sorted(SCIENCE_DIR.glob("*.fits"))
    if len(science_files) == 0:
        raise FileNotFoundError(f"No science files found in {SCIENCE_DIR}")

    df = build_science_dataframe(science_files)

    csv_path = LOG_DIR / "science_file_table.csv"
    df.to_csv(csv_path, index=False)

    print(df)
    print("\nCounts by group:")
    print(df.groupby("group_id").size())
    print(f"\nSaved table: {csv_path}")

    return df


if __name__ == "__main__":
    group_science_files()