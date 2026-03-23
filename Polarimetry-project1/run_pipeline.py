from step01_master_bias import make_master_bias
from step02_group_science import group_science_files
from step03_bias_subtract import bias_subtract_science
from step04_stack_groups import stack_groups
from step05_detect_sources import detect_sources_in_reference_group
from step06_pair_beams import pair_beams
from step07_aperture_photometry import run_photometry
from step08_compute_polarimetry import run_polarimetry
from step09_check_alignment_star import mark_alignment_star


def main():
    print("\n[1] Creating master bias")
    make_master_bias()

    print("\n[2] Grouping science files")
    group_science_files()

    print("\n[3] Bias subtracting science")
    bias_subtract_science()

    print("\n[4] Stacking groups")
    stack_groups()
    
    print("\n[4.5] Checking alignment star")
    mark_alignment_star()

    print("\n[5] Detecting sources")
    detect_sources_in_reference_group(group_id=1)

    print("\n[6] Pairing O/E beams")
    pair_beams()

    print("\n[7] Aperture photometry")
    run_photometry()

    print("\n[8] Computing polarimetry")
    run_polarimetry(chosen_aperture=5)

    print("\nPipeline completed")
    


if __name__ == "__main__":
    main()