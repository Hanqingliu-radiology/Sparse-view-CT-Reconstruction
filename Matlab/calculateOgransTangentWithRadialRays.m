function [organ_range, organ_z_range] = calculateOgransTangentWithRadialRays(mask_file)
% 使用径向射线法计算器官切线角度范围
% 输入参数:
%   mask_file - 器官分割掩膜文件路径 (NIfTI格式)
% 输出参数:
%   organ_range - 计算得到的器官角度范围
%   organ_z_range - 器官在Z轴上的范围

    % 读取掩膜文件
    mask_nii = niftiread(mask_file);
    
    % 计算Z轴范围
    slice_sums = squeeze(sum(sum(mask_nii > 0, 1), 2));
    non_zero_slices = find(slice_sums > 0);
    
    if isempty(non_zero_slices)
        organ_z_range = [1, size(mask_nii, 3)];
    else
        organ_z_range = [non_zero_slices(1), non_zero_slices(end)];
    end
    
    % 找到最大面积的切片
    [~, slice_index] = max(slice_sums);
    
    % 顺时针旋转掩膜270度（逆时针旋转90度）
    mask_nii = rot90(mask_nii, 3);
    
    % 提取最大面积切片
    mask_slice = mask_nii(:, :, slice_index);
    
    % 确保切片是二值图像
    if ~islogical(mask_slice)
        organ_mask = mask_slice > 0;
    else
        organ_mask = mask_slice;
    end
    
    % 获取图像尺寸
    [rows, cols] = size(organ_mask);
    
    % 确定图像中心点
    center_x = cols / 2;
    center_y = rows / 2;
    
    % 径向射线检测
    [ray_mask_intersections, ray_angles, boundary_points] = performRadialRayDetection(organ_mask, center_x, center_y);
    
    % 分析射线-掩膜交叉结果，找到切线点
    [left_tangent_angle, right_tangent_angle] = ...
        findBottomTangentPoints(ray_mask_intersections, ray_angles, boundary_points, center_x, center_y, rows, cols);
    
    % 计算角度差和最终范围
    angle_diff = abs(right_tangent_angle - left_tangent_angle);
    if angle_diff > 180
        angle_diff = 360 - angle_diff;
    end
    
    organ_range = angle_diff * 1.2; % 增加20%安全边界
    organ_range = min(180, organ_range);
end

function [ray_mask_intersections, ray_angles, boundary_points] = performRadialRayDetection(organ_mask, center_x, center_y)
% 执行径向射线检测
    [rows, cols] = size(organ_mask);
    
    % 定义角度步长
    angle_step = 1;
    ray_angles = 0:angle_step:(360-angle_step);
    num_rays = length(ray_angles);
    
    % 计算最大射线长度
    max_radius = max([center_x, center_y, cols-center_x, rows-center_y]) * 1.5;
    
    % 初始化结果
    ray_mask_intersections = zeros(1, num_rays);
    boundary_points = cell(1, num_rays);
    
    % 对每个角度发射射线
    for i = 1:num_rays
        angle_deg = ray_angles(i);
        angle_rad = deg2rad(angle_deg);
        
        % 计算射线方向
        dx = cos(angle_rad);
        dy = sin(angle_rad);
        
        % 沿射线采样点
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
% 基于射线检测结果找到切线点

    % 找到所有相交射线的连续区间
    intersecting_indices = find(ray_mask_intersections == 1);
    if isempty(intersecting_indices)
        error('没有射线与掩膜相交');
    end
    
    % 分析连续区间的边界
    consecutive_groups = findConsecutiveGroups(intersecting_indices);
    
    if length(consecutive_groups) == 1
        % 单一连续区域
        group = consecutive_groups{1};
        left_boundary_idx = group(1);
        right_boundary_idx = group(end);
        
        left_tangent_angle = ray_angles(left_boundary_idx);
        right_tangent_angle = ray_angles(right_boundary_idx);
        
    elseif length(consecutive_groups) >= 2
        % 多个区域 - 选择最大的两个区域的外侧边界
        group_sizes = cellfun(@length, consecutive_groups);
        [~, sorted_idx] = sort(group_sizes, 'descend');
        
        group1 = consecutive_groups{sorted_idx(1)};
        group2 = consecutive_groups{sorted_idx(2)};
        
        % 确定左右侧
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
% 找到连续的索引组
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
    
    % 处理跨越360度边界的情况
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