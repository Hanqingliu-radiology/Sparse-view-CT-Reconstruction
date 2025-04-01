import argparse
import cupy as cp
import numpy as np
import logging
from pathlib import Path
from read_data import read_dicom
from helper import save_to_tiff_stack_with_metadata
from rebinning_functions import rebin_curved_to_flat_detector, rebin_helical_to_fan_beam_trajectory

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)

def run(parser):
    args = parser.parse_args()
    logging.info(f'Processing scan {args.scan_id}.')
    raw_projections, args = read_dicom(args)
    logging.info('Finish reading DICOM data.')

    # Convert to float64 before loading to GPU for better precision
    raw_projections = raw_projections.astype(np.float64)

    try:
        # Process curved to flat detector conversion
        proj_flat_detector = rebin_curved_to_flat_detector(args, raw_projections)
        # if args.save_all:
        #     save_path = Path(args.path_out) / f'{args.scan_id}_flat_helix_projections.tif'
        #     save_to_tiff_stack_with_metadata(
        #         proj_flat_detector.get(),  # Explicitly convert to NumPy array
        #         save_path,
        #         metadata=vars(args))
        logging.info('Finished flattening projections.')
        # Process helical to fan beam conversion
        proj_fan_geometry = rebin_helical_to_fan_beam_trajectory(args, proj_flat_detector)
        logging.info('Finished fan beam conversion.')

        # Save final results
        save_path = Path(args.path_out) / f'{args.scan_id}_flat_fan_projections.tif'
        if isinstance(proj_fan_geometry, cp.ndarray):  # Check if it's a CuPy array
            proj_fan_geometry = proj_fan_geometry.get()  # Convert to NumPy array if necessary

        save_to_tiff_stack_with_metadata(
            proj_fan_geometry,  # Explicitly convert to NumPy array
            save_path,
            metadata=vars(args)
        )
        logging.info(f'Finished. Results saved at {save_path.resolve()}')

    except Exception as e:
        print(f"Error during processing: {e}")
        raise
    finally:
        # Final cleanup
        cp.get_default_memory_pool().free_all_blocks()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path_dicom', type=str, required=False,
                        default=r'F:\example-ptr-CR-raw-data\AB_NA.h5',  help='Local path of helical projection data.')
    parser.add_argument('--path_out', type=str,
                        default=r'F:\example-ptr-CR-raw-data', help='Output path of rebinned data.')
    parser.add_argument('--scan_id', type=str,
                        default='scan_001_patient', help='Custom scan ID.')
    parser.add_argument('--idx_proj_start', type=int, default=0,
                        help='First index of helical projections that are processed.')
    parser.add_argument('--idx_proj_stop', type=int, default=None,
                        help='Last index of helical projections that are processed.')
    parser.add_argument('--save_all', dest='save_all', action='store_true', default=False,
                        help='Save all intermediate results.')
    parser.add_argument('--nu', type=int, required=False, default=None,
                        help='Number of detector elements in x direction.')
    parser.add_argument('--nv', type=int, required=False, default=None,
                        help='Number of detector elements in z direction.')
    parser.add_argument('--du', type=float, required=False, default=None,
                        help='Detector element size in x direction.')
    parser.add_argument('--dv', type=float, required=False, default=None,
                        help='Detector element size in z direction.')
    parser.add_argument('--dv_rebinned', type=float, required=False, default=None,
                        help='Detector element size for rebinned data.')
    parser.add_argument('--dso', type=float, required=False,
                        default=None, help='Source-to-object distance.')
    parser.add_argument('--dsd', type=float, required=False,
                        default=None, help='Source-to-detector distance.')
    parser.add_argument('--rotview', type=int, required=False,
                        default=None, help='Number of rotation views.')
    parser.add_argument('--det_central_element', type=float, nargs=2, default=[369.625, 32.5],
                        help='Detector central element [nu, nv].')
    parser.add_argument('--hu_factor', type=float,
                        default=0.019, help='HU factor for conversion.')
    run(parser)