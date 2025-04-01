import numpy as np
import struct
import h5py

def unpack_tag(data, tag):
    return struct.unpack('f', data[tag].value)[0]


def read_projections(folder, indices):
    with h5py.File(folder, "r") as f:
        TablePosition = f['TablePosition'][()][indices]
        angles = f['angles'][()][indices]
        projections = f['projections'][()][indices,:,:]
        return projections, TablePosition, angles

def read_dicom(args):
    """ 读取 DICOM-CT-PD 投影数据 """
    indices = slice(args.idx_proj_start, args.idx_proj_stop)
    # Add this to see what's actually in the file before slicing
    with h5py.File(args.path_dicom, "r") as h5_file:
        actual_proj_count = h5_file['projections'].shape[0]
        print(f"H5 file contains {actual_proj_count} projections")
        print(f"Requested range: {args.idx_proj_start} to {args.idx_proj_stop}")

    projections, z_positions, angles = read_projections(args.path_dicom, indices)

    # Add this to see the actual returned data size
    print(f"Loaded projections shape: {projections.shape}")

    # 角度转换（确保在 with 作用域内计算）
    angles = angles*2*np.pi/360
    angles = - np.unwrap(angles) - np.pi  # 计算 angles 确保文件未关闭

    # 传递从命令行传递的参数，并赋值
    args.nu = args.nu or 736
    args.nv = args.nv or 64
    args.du = args.du or 1.2858393
    args.dv = args.dv or 1.0947227
    args.dv_rebinned = args.dv_rebinned or 0.6
    args.det_central_element = args.det_central_element or [369.625, 32.5]
    args.dso = args.dso or 595
    args.dsd = args.dsd or 1085.6
    args.ddo = args.dsd - args.dso
    args.pitch = (z_positions[2304] - z_positions[0]).item()
    args.nz_rebinned = int((z_positions[-1] - z_positions[0]) / args.dv_rebinned)
    args.hu_factor = args.hu_factor or 0.019
    args.rotview = args.rotview or 2304

    # 将 angles 和 z_positions 转为 list，方便存储
    args.angles = angles.tolist()
    args.z_positions = z_positions.tolist()

    return projections, args