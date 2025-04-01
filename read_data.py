import numpy as np
import h5py

def read_dicom(args):
    with h5py.File(args.path_dicom, "r", swmr=True) as h5_file:
        actual_proj_count = h5_file['projections'].shape[0]
        print(f"H5 file contains {actual_proj_count} projections")

        if not hasattr(args, 'idx_proj_stop') or args.idx_proj_stop is None or args.idx_proj_stop > actual_proj_count:
            args.idx_proj_stop = actual_proj_count
            print(f"Setting idx_proj_stop to {actual_proj_count}")

        if not hasattr(args, 'idx_proj_start') or args.idx_proj_start is None:
            args.idx_proj_start = 0
        print(f"Using projection range: {args.idx_proj_start} to {args.idx_proj_stop}")

        indices = slice(args.idx_proj_start, args.idx_proj_stop)
        projections = np.array(h5_file['projections'][indices, :, :])
        z_positions = np.array(h5_file['TablePosition'][indices])
        angles_raw = np.array(h5_file['angles'][indices])
        angles = -np.unwrap(angles_raw * 2 * np.pi / 360) - np.pi

    param_defaults = {
        'nu': 736,
        'nv': 64,
        'du': 1.2858393,
        'dv': 1.0947227,
        'dv_rebinned': 0.6,
        'det_central_element': [369.625, 32.5],
        'dso': 595,
        'dsd': 1085.6,
        'hu_factor': 0.019,
        'rotview': 2304,
    }

    for param, default in param_defaults.items():
        if not getattr(args, param, None):
            setattr(args, param, default)

    args.ddo = args.dsd - args.dso
    args.pitch = (z_positions[2304] - z_positions[0]).item() if len(z_positions) > 2304 else 0
    args.nz_rebinned = int((z_positions[-1] - z_positions[0]) / args.dv_rebinned)

    args.angles = angles.tolist()
    args.z_positions = z_positions.tolist()

    return projections, args