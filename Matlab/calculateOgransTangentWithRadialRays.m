function [organ_range, organ_z_range] = calculateOgransTangentWithRadialRays(mask_file)
% Calculate organ tangent angle range using radial ray method
% Input parameters:
%   mask_file - Organ segmentation mask file path (NIfTI format)
% Output parameters:
%   organ_range - Calculated organ angle range
%   organ_z_range - Organ range along Z-axis

    % Read mask file
    mask_nii = niftiread(mask_file);
    
    % Calculate Z-axis range
    slice_sums = squeeze(sum(sum(mask_nii > 0, 1), 2));
    non_zero_slices = find(slice_sums > 0);
    
    if isempty(non_zero_slices)
        organ_z_range = [1, size(mask_nii, 3)];
    else
        organ_z_range = [non_zero_slices(1), non_zero_slices(end)];
    end
    
    % Find slice with maximum area
    [~, slice_index] = max(slice_sums);
    
    % Rotate mask clockwise by 270 degrees (counterclockwise by 90 degrees)
    mask_nii = rot90(mask_nii, 3);
    
    % Extract slice with maximum area
    mask_slice = mask_nii(:, :, slice_index);
    
    % Ensure slice is binary image
    if ~islogical(mask_slice)
        organ_mask = mask_slice > 0;
    else
        organ_mask = mask_slice;
    end
    
    % Get image dimensions
    [rows, cols] = size(organ_mask);
    
    % Determine image center point
    center_x = cols / 2;
    center_y = rows / 2;
    
    % Radial ray detection
    [ray_mask_intersections, ray_angles, boundary_points] = performRadialRayDetection(organ_mask, center_x, center_y);
    
    % Analyze ray-mask intersection results to find tangent points
    [left_tangent_angle, right_tangent_angle] = ...
        findBottomTangentPoints(ray_mask_intersections, ray_angles, boundary_points, center_x, center_y, rows, cols);
    
    % Calculate angle difference and final range
    angle_diff = abs(right_tangent_angle - left_tangent_angle);
    if angle_diff > 180
        angle_diff = 360 - angle_diff;
    end
    
    organ_range = angle_diff * 1.2; % Add 20% safety margin
    organ_range = min(180, organ_range);
end

function [ray_mask_intersections, ray_angles, boundary_points] = performRadialRayDetection(organ_mask, center_x, center_y)
% Perform radial ray detection
    [rows, cols] = size(organ_mask);
    
    % Define angle step size
    angle_step = 1;
    ray_angles = 0:angle_step:(360-angle_step);
    num_rays = length(ray_angles);
    
    % Calculate maximum ray length
    max_radius = max([center_x, center_y, cols-center_x, rows-center_y]) * 1.5;
    
    % Initialize results
    ray_mask_intersections = zeros(1, num_rays);
    boundary_points = cell(1, num_rays);
    
    % Cast rays for each angle
    for i = 1:num_rays
        angle_deg = ray_angles(i);
        angle_rad = deg2rad(angle_deg);
        
        % Calculate ray direction
        dx = cos(angle_rad);
        dy = sin(angle_rad);
        
        % Sample points along the ray
        intersection_points = [];
        
        for r = 1:max_radius
            x = center_x + r * dx;
            y = center_y + r * dy;
            
            if x < 1 || x > cols || y < 1 || y > rows
                break;
            end
            
            xi = round(x);
            yi = round(y);
            
            if xi >= 1 && xi <= cols && yi >= 1 && yi <= rows
                if organ_mask(yi, xi) > 0
                    intersection_points(end+1, :) = [yi, xi];
                end
            end
        end
        
        if ~isempty(intersection_points)
            ray_mask_intersections(i) = 1;
            boundary_points{i} = intersection_points;
        else
            ray_mask_intersections(i) = 0;
            boundary_points{i} = [];
        end
    end
end

function [left_tangent_angle, right_tangent_angle] = ...
    findBottomTangentPoints(ray_mask_intersections, ray_angles, boundary_points, center_x, center_y, rows, cols)
% Find tangent points based on ray detection results

    % Find continuous intervals of all intersecting rays
    intersecting_indices = find(ray_mask_intersections == 1);
    if isempty(intersecting_indices)
        error('No rays intersect with the mask');
    end
    
    % Analyze boundaries of continuous intervals
    consecutive_groups = findConsecutiveGroups(intersecting_indices);
    
    if length(consecutive_groups) == 1
        % Single continuous region
        group = consecutive_groups{1};
        left_boundary_idx = group(1);
        right_boundary_idx = group(end);
        
        left_tangent_angle = ray_angles(left_boundary_idx);
        right_tangent_angle = ray_angles(right_boundary_idx);
        
    elseif length(consecutive_groups) >= 2
        % Multiple regions - select outer boundaries of the two largest regions
        group_sizes = cellfun(@length, consecutive_groups);
        [~, sorted_idx] = sort(group_sizes, 'descend');
        
        group1 = consecutive_groups{sorted_idx(1)};
        group2 = consecutive_groups{sorted_idx(2)};
        
        % Determine left and right sides
        center1 = mean(ray_angles(group1));
        center2 = mean(ray_angles(group2));
        
        if center1 < center2
            left_group = group1;
            right_group = group2;
        else
            left_group = group2;
            right_group = group1;
        end
        
        left_tangent_angle = ray_angles(left_group(1));
        right_tangent_angle = ray_angles(right_group(end));
    end
end

function consecutive_groups = findConsecutiveGroups(indices)
% Find consecutive index groups
    consecutive_groups = {};
    if isempty(indices)
        return;
    end
    
    current_group = indices(1);
    
    for i = 2:length(indices)
        if indices(i) == indices(i-1) + 1
            current_group(end+1) = indices(i);
        else
            consecutive_groups{end+1} = current_group;
            current_group = indices(i);
        end
    end
    
    consecutive_groups{end+1} = current_group;
    
    % Handle cases crossing the 360-degree boundary
    if length(consecutive_groups) > 1
        first_group = consecutive_groups{1};
        last_group = consecutive_groups{end};
        
        if first_group(1) == 1 && last_group(end) == 360
            merged_group = [last_group, first_group];
            consecutive_groups = consecutive_groups(2:end-1);
            consecutive_groups{end+1} = merged_group;
        end
    end
end
