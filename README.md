
# Sparse-view-CT-Reconstruction
Welcome to the organ-based Sparse-view-CT-Reconstruction, this projection mainly focuses on the radiation-sensitive organs such as the eye lens, thyroid, and breast, and so on.
## Convert raw data to h5 files. 
## Convert h5 files to tif files
1. Generate a bat file in the format:
```
@echo off 
python main.py --path_dicom patient_0001 --path_out save_path --scan_id patient_0001 --save_all & ^
python main.py --path_dicom patient_0002 --path_out save_path --scan_id patient_0002 --save_all & ^
python main.py --path_dicom patient_0003 --path_out save_path --scan_id patient_0003 --save_all
```
2. Run the batch file to batch process the data.
## 

