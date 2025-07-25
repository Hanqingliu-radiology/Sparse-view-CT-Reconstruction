function [reconstructed_image, metaData, angles] = ogran_based_sparse_view_reconstruction(sinogram_tiff_path, recostructed_url, ...
    full_body_mask, eye_lens_mask, thyroid_mask, breast_mask, ...
    reconstruction_method, number_of_iterations, sparse_ratio, z_start, z_end, image_dimension, tigre_affine_z)

%% 1. Read file information
[metaData, angles] = ReadFanInfo(sinogram_tiff_path);
if isequal(z_end, "all")
    z_end = metaData.nz_rebinned;
end

Volume_tiff = tiffreadVolume(sinogram_tiff_path); 
Volume_tiff = flip(Volume_tiff, 2);
Volume_tiff = Volume_tiff / 2294.5;

%% 2. Calculate the difference in the number of layers between Siemens Image and Tigre Image
full_body_info = niftiinfo(full_body_mask);
siemens_affine_z = full_body_info.Transform.T(4,3);
Decimal = siemens_affine_z - floor(siemens_affine_z);
tigre_affine_z = fix(tigre_affine_z) - Decimal;
bias_slice = siemens_affine_z - tigre_affine_z;

%% 3. Calculate angle ranges and Z-axis ranges for each tissue
fprintf('=====================================================================\n');
fprintf('Calculating angle ranges and Z-axis ranges for each tissue\n');

% Initialize tissue information structure
organ_info = struct();

% Define organ list and corresponding mask files
organs = {'eye_lens', 'thyroid', 'breast'};
organ_names = {'Eye Lens', 'Thyroid', 'Breast'};
mask_files = {eye_lens_mask, thyroid_mask, breast_mask};

% Use auxiliary function to calculate each organ's range
for i = 1:length(organs)
    organ_info.(organs{i}) = calculateOrganRange(mask_files{i}, organs{i}, organ_names{i}, bias_slice);
end

fprintf('=====================================================================\n');

%% 4. Reconstruction process - Key fix: following full_dose pattern
fprintf('Starting organ-specific sparse view reconstruction\n');
fprintf('Unified sparse ratio: %.1f%% (keeping %.1f%% of angles)\n', sparse_ratio*100, (1-sparse_ratio)*100);
% reconstructed_z_slice = size(niftiread(niftiinfo(full_body_mask)), 3);
% Fix 1: Use the same data type as full_dose
reconstructed_image = zeros(image_dimension, image_dimension, z_end - z_start + 1, 'int16');
% Create progress bar
h = waitbar(0, 'Reconstruction in progress...');

for axial_slice = 1:size(reconstructed_image, 3)
    % Update progress bar
    if isvalid(h)
        waitbar(axial_slice/size(reconstructed_image, 3), h, ...
            sprintf('Reconstructing slice %d/%d', axial_slice, size(reconstructed_image, 3)));
    end
        % prog_bar_yazdan(axial_slice, size(reconstructed_image,3), 'Reconstruction    ', 'sdfdsf')
    
    % Get the actual Z coordinate of current slice
    current_z = z_start + axial_slice - 1;
    
    % Determine which organ region the current slice belongs to
    organ_type = 'background';
    angle_range = 360; % Default to full angle
    
    % Check if in eye lens region
    if organ_info.eye_lens.valid && ...
       current_z >= organ_info.eye_lens.z_range(1) && ...
       current_z <= organ_info.eye_lens.z_range(2)
        organ_type = 'eye_lens';
        angle_range = organ_info.eye_lens.angle_range;
    % Check if in thyroid region
    elseif organ_info.thyroid.valid && ...
           current_z >= organ_info.thyroid.z_range(1) && ...
           current_z <= organ_info.thyroid.z_range(2)
        organ_type = 'thyroid';
        angle_range = organ_info.thyroid.angle_range;
    % Check if in breast region
    elseif organ_info.breast.valid && ...
           current_z >= organ_info.breast.z_range(1) && ...
           current_z <= organ_info.breast.z_range(2)
        organ_type = 'breast';
        angle_range = organ_info.breast.angle_range;
    end
    
    % Calculate angles used for current slice
    if strcmp(organ_type, 'background')
        % Background region uses all angles
        selected_indices = 1:length(angles);
        selected_angles = angles;
    else
        % Organ region uses sparse sampling
        [selected_angles, selected_indices] = selectSparseAngles(angles, angle_range, sparse_ratio, organ_type);
    end
    
    % Fix 2: Ensure sufficient angles for reconstruction
    if length(selected_angles) < 10  % If too few angles, use all angles
        fprintf('Warning: Slice #%d has too few selected angles (%d), using full angle reconstruction\n', axial_slice, length(selected_angles));
        selected_indices = 1:length(angles);
        selected_angles = angles;
    end

    % Extract projection data for current slice - Fix 3: Following full_dose approach
    projections = Volume_tiff(:, axial_slice, selected_indices);
    projections = squeeze(projections);
    
    % Display debug information
    if mod(axial_slice, 20) == 0
        fprintf('Slice #%d (Z=%d, %s): Using %d/%d projection angles\n', ...
            axial_slice, current_z, organ_type, length(selected_angles), length(angles));
        fprintf('Projection data range: [%.3f, %.3f]\n', min(projections(:)), max(projections(:)));
    end
    
    % Perform reconstruction - Fix 4: Following the full_dose complete reconstruction method list
    try
        if isequal(reconstruction_method, "MLEM")
            img = MLEM_patient(projections, selected_angles, number_of_iterations);
        elseif isequal(reconstruction_method, "FDK")
            img = FDK_patient(projections, selected_angles);
        elseif isequal(reconstruction_method, "FBP")
            img = FDK_patient(projections, selected_angles);
        elseif isequal(reconstruction_method, "OSSART")
            img = OSSART_patient(projections, selected_angles, number_of_iterations);
        elseif isequal(reconstruction_method, "SIRT")
            img = SIRT_patient(projections, selected_angles, number_of_iterations);   
        else
            error('Unsupported reconstruction method: %s', reconstruction_method);
        end
        
        % Fix 5: Completely following the full_dose processing approach
        img_rotated = rot90(img, 1);
        % Convert image values to standard CT value range and round to integers
        img_ct = round(img_rotated);  % Round to nearest integer
        % Ensure CT values are within a reasonable range (-1000 to 3000 HU) - Fix 6: Use the same range as full_dose
        img_ct = max(-1000, min(3000, img_ct));
        reconstructed_image(:, :, axial_slice) = int16(img_ct);
        
        % Debug: Check reconstruction results
        if mod(axial_slice, 20) == 0
            fprintf('Reconstruction result range: [%.1f, %.1f]\n', min(img_ct(:)), max(img_ct(:)));
        end
        
    catch ME
        warning('Failed to reconstruct slice #%d: %s', axial_slice, ME.message);
        reconstructed_image(:, :, axial_slice) = zeros(image_dimension, image_dimension, 'int32');
    end
end
% Close progress bar
if isvalid(h)
    close(h);
end
% Fix 7: Following full_dose order processing
% reconstructed_image = flip(reconstructed_image, 3);

%% 5. Post-processing and saving - Completely following the full_dose approach
reference_spacing = [0.9766, 0.9766, 1];
reconstructed_info = raw_nifti_info(3);
reconstructed_info.PixelDimensions = reference_spacing;
reconstructed_info.ImageSize = nifti_info_size(size(reconstructed_image));

original_size = size(reconstructed_image);
original_z_slices = original_size(3);

% Use same Z-axis resampling as full_dose
target_z_slices = z_end;
try
    % [X, Y, Z] = ndgrid(1:original_size(1), 1:original_size(2), linspace(1, original_z_slices, reconstructed_z_slice));
    [X, Y, Z] = ndgrid(1:original_size(1), 1:original_size(2), linspace(1, original_z_slices, target_z_slices));
    % Fix 8: Following the full-dose resampling approach
    reconstructed_double = double(reconstructed_image);
    corrected_image_double = interpn(reconstructed_double, X, Y, Z, 'linear');
    corrected_image = int16(round(corrected_image_double));  % Convert back to integer
catch err
    disp(['Z-axis resampling error: ', err.message]);
    corrected_image = zeros(original_size(1), original_size(2), target_z_slices, 'int32');
    for z = 1:target_z_slices
        orig_z = min(round((z * original_z_slices) / target_z_slices), original_z_slices);
        if orig_z < 1
            orig_z = 1;
        end
        corrected_image(:,:,z) = reconstructed_image(:,:,orig_z);
    end
end

% Replace original image with resampled image
reconstructed_image = corrected_image;
clear corrected_image corrected_image_double reconstructed_double;

% Fix 9: Completely following full_dose circular mask processing
[rows, cols, ~] = size(reconstructed_image);
center_x = cols / 2;
center_y = rows / 2;
radius = min(rows, cols) / 2;

[X, Y] = meshgrid(1:cols, 1:rows);
distance_from_center = sqrt((X - center_x).^2 + (Y - center_y).^2);
circular_mask = distance_from_center <= radius;

% Apply circular mask to all slices
for slice_idx = 1:size(reconstructed_image, 3)
    current_slice = reconstructed_image(:, :, slice_idx);
    % Fix 10: Use the same mask value as full_dose
    current_slice(~circular_mask) = int16(-999);
    reconstructed_image(:, :, slice_idx) = current_slice;
end

% Fix 11: Completely following the full-dose saving approach
% Read reference file information (use first valid mask file)
segmentation_file = full_body_mask;
if organ_info.breast.valid && exist(breast_mask, 'file')
    segmentation_file = breast_mask;
elseif organ_info.thyroid.valid && exist(thyroid_mask, 'file')
    segmentation_file = thyroid_mask;
elseif organ_info.eye_lens.valid && exist(eye_lens_mask, 'file')
    segmentation_file = eye_lens_mask;
end

% ===== Save NIfTI file - Using exactly the same method as full_dose =====
try
    if ~isempty(segmentation_file) && exist(segmentation_file, 'file')
        fprintf('Using reference file as template for saving...\n');
        
        template_info = niftiinfo(segmentation_file);
        template_info.Filename = [recostructed_url, '.nii.gz'];
        template_info.ImageSize = size(reconstructed_image);
        template_info.Datatype = 'int16';
        template_info.Description = sprintf('Organ-based reconstruction using %s method', reconstruction_method);
        
        img_size = size(reconstructed_image);
        template_info.raw.dim(2:4) = img_size;
        template_info.Transform.T(4,3) = tigre_affine_z;
        
        niftiwrite(reconstructed_image, recostructed_url, template_info, 'Compressed', true);
        
        fprintf('Successfully saved file: %s.nii.gz\n', recostructed_url);
        
        if exist([recostructed_url, '.nii.gz'], 'file')
            saved_info = niftiinfo([recostructed_url, '.nii.gz']);
            fprintf('Verified saved Origin: [%.1f, %.1f, %.1f]\n', ...
                saved_info.raw.qoffset_x, saved_info.raw.qoffset_y, saved_info.raw.qoffset_z);
            fprintf('Verified saved dimensions: %s\n', mat2str(saved_info.ImageSize));
            fprintf('Verified saved data type: %s\n', saved_info.Datatype);
        end
        
    else
        fprintf('No reference file, saving with default template...\n');
        niftiwrite(reconstructed_image, recostructed_url, 'Compressed', true);
        fprintf('Successfully saved file (default template): %s.nii.gz\n', recostructed_url);
    end
    
catch err
    fprintf('Failed to save NIfTI file: %s\n', err.message);
    fprintf('Trying backup saving method...\n');
    
    try
        niftiwrite(reconstructed_image, recostructed_url);
        fprintf('Successfully saved file (uncompressed): %s.nii\n', recostructed_url);
    catch err2
        fprintf('All saving methods failed: %s\n', err2.message);
    end
end

% Clean up memory
try
    astra_clear;
catch
end
clear Volume_tiff proj_geom vol_geom projections

fprintf('\n=====================================================================\n');
fprintf('Organ-specific sparse view reconstruction summary:\n');
fprintf('Unified sparse ratio: %.1f%% (keeping %.1f%% of angles)\n', sparse_ratio*100, (1-sparse_ratio)*100);
if organ_info.eye_lens.valid
    fprintf('- Eye lens: Angle range %.1f°, Z-axis [%d-%d]\n', ...
        organ_info.eye_lens.angle_range, organ_info.eye_lens.z_range(1), ...
        organ_info.eye_lens.z_range(2));
end
if organ_info.thyroid.valid
    fprintf('- Thyroid: Angle range %.1f°, Z-axis [%d-%d]\n', ...
        organ_info.thyroid.angle_range, organ_info.thyroid.z_range(1), ...
        organ_info.thyroid.z_range(2));
end
if organ_info.breast.valid
    fprintf('- Breast: Angle range %.1f°, Z-axis [%d-%d]\n', ...
        organ_info.breast.angle_range, organ_info.breast.z_range(1), ...
        organ_info.breast.z_range(2));
end
fprintf('Reconstructed image dimensions: [%s]\n', num2str(size(reconstructed_image)));
fprintf('NIfTI file saved to: %s\n', recostructed_url);
fprintf('=====================================================================\n');

end

%% Auxiliary function 1: Calculate organ range
function organ_data = calculateOrganRange(mask_file, ~, organ_name_en, bias_slice)
    organ_data = struct();  
    if nargin < 4
        bias_slice = 0;
    end
    
    if ~isempty(mask_file) && exist(mask_file, 'file')
        fprintf('\n--- Calculating %s angle range (bias: %d) ---\n', organ_name_en, bias_slice);
        try
            [angle_range, z_range] = calculateOgransTangentWithRadialRays(mask_file);    
            z_range_biased = z_range + bias_slice;
            organ_data.angle_range = angle_range;
            organ_data.z_range = z_range_biased;  
            organ_data.z_range_original = z_range;  
            organ_data.bias_slice = bias_slice;  
            organ_data.valid = true;
            
            fprintf('%s angle range: %.1f degrees, Z-axis range: [%d, %d] (original: [%d, %d], bias: %d)\n', ...
                organ_name_en, angle_range, z_range_biased(1), z_range_biased(2), ...
                z_range(1), z_range(2), bias_slice);
                
        catch ME
            warning('Error calculating %s angle: %s', organ_name_en, ME.message);

            organ_data.valid = false;
        end
    else
        organ_data.valid = false;
    end
end

%% Auxiliary function 2: Improved sparse angle selection
function [selected_angles, selected_indices] = selectSparseAngles(all_angles, angle_range, sparse_ratio, organ_type)
    if strcmp(organ_type, 'background') || angle_range >= 360
        selected_indices = 1:length(all_angles);
        selected_angles = all_angles;
        return;
    end
    
    total_angles = length(all_angles);
    keep_ratio = 1 - sparse_ratio;  
    min_angles = max(10, round(total_angles * 0.1));  
    target_angles = max(min_angles, round(total_angles * keep_ratio));
    
    if target_angles >= total_angles
        selected_indices = 1:total_angles;
    else
        step = max(1, floor(total_angles / target_angles));
        selected_indices = 1:step:total_angles;
        if length(selected_indices) > target_angles
            selected_indices = selected_indices(1:target_angles);
        end
    end
    
    selected_angles = all_angles(selected_indices);
    
    fprintf('Organ %s: Total angles %d, Target angles %d, Actual select %d\n', ...
        organ_type, total_angles, target_angles, length(selected_indices));
end
