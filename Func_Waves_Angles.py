import numpy as np
from scipy.ndimage import binary_erosion,  binary_dilation
from math import atan2, degrees
from scipy.spatial.distance import cdist
from scipy.spatial import distance_matrix
from skimage.measure import regionprops
from skimage.segmentation import flood
from scipy.ndimage import label as label
from skimage.measure import label as lab
from skimage.morphology import skeletonize
from skimage.draw import line

def is_island(island_mask, lake):
    """
    Check if a land blob is a real island (surrounded by water).
    """
    on_border = (
        np.any(island_mask[0, :]) or
        np.any(island_mask[-1, :]) or
        np.any(island_mask[:, 0]) or
        np.any(island_mask[:, -1])
    )
    return not on_border and np.all(lake[island_mask] == 0)


def find_nearest_non_zero(area_maskTemp, middle_point, radius=2):
    y, x = middle_point  # Middle point (row, col)

    # Start with an empty list to track all the neighboring coordinates to check
    # Define the offsets for the square search radius (from -radius to +radius)
    for r in range(-radius, radius + 1):
        for c in range(-radius, radius + 1):
            # Calculate new coordinates
            new_y, new_x = y + r, x + c

            # Check if new coordinates are within the bounds of the array
            if 0 <= new_y < area_maskTemp.shape[0] and 0 <= new_x < area_maskTemp.shape[1]:
                if area_maskTemp[new_y, new_x] == 1:
                    return new_y, new_x  # Return first matching point

    return None  # Return None if no point was found
def split_and_label_area_mask_by_islands(lake, area_mask, propagated_lake, start_pt):
    """
    Split area_mask by cutting at island locations and return a labeled mask of subregions.

    Parameters:
    - lake: 2D array (1 = water, 0 = land)
    - area_mask: 2D boolean array to be split
    - start_point: (x, y) = (col, row)

    Returns:
    - labeled_mask: 2D integer array with separate labels for each area chunk after island cuts
    """
    # 1. Intersection of the two rasters
    # Binary masks
    maskA = area_mask == 1
    maskB = propagated_lake == 1
    print('size area', np.sum(area_mask))
    # Dilate each mask slightly
    dilatedA = binary_dilation(maskA)
    dilatedB = binary_dilation(maskB)

    # Touching zone = where A is adjacent to B (but not overlapping)
    touch_zone = dilatedA & dilatedB & ~(maskA & maskB)  # just in case of overlap (safety)    print(intersection)

    # 2. Optional cleaning: remove small objects or noise
    # e.g., from skimage.morphology import remove_small_objects
    # intersection = remove_small_objects(intersection, min_size=20)

    # 3. Skeletonize to get the centerline
    skeleton = skeletonize(touch_zone)

    # 4. Label connected components (lines)
    labeled = lab(skeleton)

    # 5. For each line, extract endpoints
    endpoints = []
    for region in regionprops(labeled):
        coords = region.coords
        if len(coords) < 2:
            continue

        # Use a distance matrix to find the two furthest points
        dists = np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=-1)
        i, j = np.unravel_index(np.argmax(dists), dists.shape)

        pointA = coords[i]
        pointB = coords[j]
        endpoints.append((tuple(pointA[::-1]), tuple(pointB[::-1])))  # reverse (row, col) to (x, y)

    At = []
    Bt = []
    i = 0
    start_point = [start_pt[1], start_pt[0]]
    for a, b in endpoints:

        # Compute distances
        dist_a = np.linalg.norm(np.array(a) - np.array(start_point))
        dist_b = np.linalg.norm(np.array(b) - np.array(start_point))
        # Assign A and B
        if dist_a <= dist_b:
            At.append(np.array(a))
            Bt.append(np.array(b))
        else:
            At.append(np.array(b))
            Bt.append(np.array(a))
        # aT = At[i]
        # bT = At[i]
        i += 1

    if At == []:
        return At, Bt, area_mask, propagated_lake
    AB_pairs = list(zip(At, Bt))
    # Sort pairs by distance from A to start_point
    AB_pairs.sort(key=lambda pair: np.linalg.norm(pair[0] - np.array(start_point)))

    # Unzip back to A and B
    AL, BL = zip(*AB_pairs)
    AL = list(AL)
    BL = list(BL)

    # New filtered lists
    AL_filtered = []
    BL_filtered = []

    area_maskN = area_mask.copy()
    area_maskTemp = area_mask.copy()

    area_maskN = area_maskN[np.newaxis, :, :]
    j = 0
    for A, B in zip(AL, BL):
        MidllePointX = (A[0] + B[0]) / 2
        MidllePointY = (A[1] + B[1]) / 2
        MiddlePoint = [int(MidllePointX), int(MidllePointY)]

        if np.sum(area_maskTemp) == 0:
            print('skip Mask')

            continue
        if area_maskTemp[MiddlePoint[1], MiddlePoint[0]] == 0:

            try:
                MiddlePointProp = [MiddlePoint[1], MiddlePoint[0]]
                new_y, new_x = find_nearest_non_zero(area_maskTemp, MiddlePointProp, radius=2)
                MiddlePointProp = [new_y, new_x]

            except:
                print('other side of island')
                continue
        else:
            MiddlePointProp = [MiddlePoint[1], MiddlePoint[0]]


        ####   change area mask_temp to area mask to then check overlap ####
        propagated_lakeT = propagate_rays(area_mask, MiddlePointProp, max_fill_distance=1)
        # best_point, best_area_size, propagated_lakeT = propagate_rays_from_line(area_maskTemp, A, B)

        AL_filtered.append(A)
        BL_filtered.append(B)


        if propagated_lake.ndim >= 3:
            propagated_lakeT = propagated_lakeT.squeeze()
        ###################################################################################################
        area_maskTemp = area_mask  #  np.logical_xor(area_maskTemp, propagated_lakeT)
        ###################################################################################################

        if j ==0:
            propagated_lakes = propagated_lakeT
            propagated_lakes = propagated_lakes[np.newaxis, :, :]
        else:
            propagated_lakes = np.concatenate((propagated_lakes,  propagated_lakeT[np.newaxis, :, :]), axis=0)
        j += 1
    AL = AL_filtered
    BL = BL_filtered

    structure = np.ones((3, 3))  # 8-connectivity
    remaining_areas, num_features = label((area_maskTemp == 1) & (propagated_lake == 0), structure)

    for i in range(len(propagated_lakes)):
        PropTemp = propagated_lakes[i,:,:]
        dilated = binary_dilation(PropTemp == 1, structure=structure)

        # Check which labels in remaining_areas touch the propagated lake
        touching_labels = np.unique(remaining_areas[dilated & (remaining_areas > 0)])

        # Create a mask of all touching regions
        touching_mask = np.isin(remaining_areas, touching_labels)
        remaining_areas[np.isin(remaining_areas, touching_labels)] = 0
        # Merge with the current propagated lake layer
        merged = ((PropTemp == 1) | touching_mask).astype(np.uint8)

        if i == 0:
            area_maskN[0,:,:] = merged
        else:
            area_maskN = np.concatenate((area_maskN,  merged[np.newaxis, :, :]), axis=0)

    return AL, BL, area_maskN, propagated_lakes

def erase_Islands(lake, size_lim = 100000):
    lake = lake.astype(np.uint8)
    land = (lake == 0)

    labeled_islands, num_islands = label(land)

    merged_mask = lake

    for island_label in range(1, num_islands + 1):
        island_mask = (labeled_islands == island_label)

        if np.sum(island_mask) > size_lim:
            print(island_label, np.sum(island_mask))
            continue

        if is_island(island_mask, lake):
            merged_mask += island_mask
    return merged_mask



def move_point_inside(area_mask, point):
    """
    Ensures the point is inside the area mask. If the point is outside,
    moves it to the closest point inside the area.

    Parameters:
    - area_mask: The binary mask (1 for area, 0 for non-area).
    - point: The (x, y) coordinates of the point.

    Returns:
    - New point coordinates if moved inside, else the original point.
    """
    # Round the point to integer coordinates
    x, y = np.round(point).astype(int)

    # Check if the point is inside the area
    if 0 <= x < area_mask.shape[0] and 0 <= y < area_mask.shape[1] and area_mask[x, y] == 1:
        return x, y  # Point is already inside, do nothing

    # If the point is outside the area, find the closest point inside the area
    area_coords = np.column_stack(np.where(area_mask == 1))  # Get all (x, y) coordinates inside the area
    point_coord = np.array([x, y])

    # Compute distances from the point to all points inside the area
    distances = distance_matrix([point_coord], area_coords)

    # Find the closest point inside the area
    closest_point_idx = np.argmin(distances)
    closest_point = area_coords[closest_point_idx]

    return closest_point


# Helper function: Bresenham's line algorithm
def bresenham(x1, y1, x2, y2):
    """Returns the points on a line between two coordinates using Bresenham's algorithm."""
    points = []
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx - dy
    while (x1 != x2 or y1 != y2):
        points.append((x1, y1))
        e2 = err * 2
        if e2 > -dy:
            err -= dy
            x1 += sx
        if e2 < dx:
            err += dx
            y1 += sy
    points.append((x2, y2))  # Include last point
    return points


# Find the borders of the lake
def find_borders(lake):
    structure = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
    lake = lake.astype(bool)
    lake_eroded = binary_erosion(lake, structure=structure).astype(lake.dtype)
    borders = lake & ~lake_eroded
    # borders = lake - lake_eroded
    return np.argwhere(borders == 1)


# Cast rays from the start point to all border points
def propagate_rays(lake, start_point, max_fill_distance=1):
    # print(start_point)
    propagated_lake = np.zeros_like(lake)
    border_points = find_borders(lake)
    x, y = zip(*border_points)
    # print('here')
    # plt.figure()
    # plt.scatter(x,y)
    # plt.scatter(start_point[0], start_point[1])
    # # plt.imshow(propagated_lake_filled)
    # plt.show()

    for border_point in border_points:
        x2, y2 = border_point
        x1, y1 = start_point
        rr, cc = line(y1, x1, y2, x2)  # note: (y, x) order
        line_points= list(zip(cc, rr))
        # line_points = bresenham(x1, y1, x2, y2)
        for (x, y) in line_points:
            if lake[x, y] == 1:  # Stop at the border
                propagated_lake[x, y] = 1
            else:
                break
    propagated_lake_filled = fill_isolated_pixels(lake, propagated_lake, max_fill_distance=max_fill_distance)

    return propagated_lake_filled


def fill_isolated_pixels(lake, propagated_lake, max_fill_distance=1):
    """
    Fills isolated water pixels adjacent to both the 100% area and the ground with 100% fill.
    Allows tuning the number of pixels to fill by increasing max_fill_distance.

    Parameters:
    - lake: The original binary map of the lake (1 for water, 0 for land).
    - propagated_lake: The current state of the propagated lake (1 for 100% filled, 0 for unfilled).
    - max_fill_distance: The number of pixel layers to consider for filling. Defaults to 1 layer (direct adjacency).

    Returns:
    - propagated_lake: Updated lake with isolated pixels filled.
    """
    # Define the neighborhood for checking adjacency
    structure = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])

    # Find the pixels adjacent to the 100% propagated lake area
    adjacency_to_100_percent = binary_dilation(propagated_lake, structure=structure,
                                               iterations=max_fill_distance).astype(lake.dtype)

    # Find the pixels adjacent to the ground (land) area
    adjacency_to_ground = binary_dilation(lake == 0, structure=structure, iterations=max_fill_distance).astype(
        lake.dtype)

    # Isolate pixels that are water, not part of the 100% area, adjacent to both the 100% area and the ground
    isolated_pixels = (lake == 1) & (propagated_lake == 0) & (adjacency_to_100_percent == 1) & (
            adjacency_to_ground == 1)

    # Set these isolated pixels to 100% area
    propagated_lake[isolated_pixels == 1] = 1

    return propagated_lake


# Compute the angles for each cell in the area
def compute_angle_for_cell(A, B, C):
    """Compute the angle ABC."""
    BA = np.array(A) - np.array(B)
    BC = np.array(C) - np.array(B)
    angle = atan2(BA[0] * BC[1] - BA[1] * BC[0], BA[0] * BC[0] + BA[1] * BC[1])
    return degrees(angle)


def propagate_rays_from_line(lake, start_point, end_point):
    """Propagate waves from each point along the line defined by start_point and end_point."""
    propagated_lake = np.zeros_like(lake)  # Create an empty array to store the propagated area
    # start_point = [start_point[1], start_point[0]]
    rr, cc = line(start_point[1], start_point[0], end_point[1], end_point[0])  # note: (y, x) order
    line_points = list(zip(cc, rr))
    # line_points = bresenham(start_point[1], start_point[0], end_point[1], end_point[0])  # Get the line of points
    # line_points = bresenham(start_point[0], start_point[1], end_point[0], end_point[1])  # Get the line of points
    line_points = line_points[0::]
    best_area_size = 0  # Variable to track the largest area
    best_point = None  # Variable to store the point that gives the largest area
    # plt.figure()
    # plt.imshow(lake)
    # plt.scatter(start_point[0], start_point[1])
    # plt.show()
    for point in line_points:
        point1 = [point[1], point[0]]
        shifted_point = point1
        x1, y1 = shifted_point  # The current point on the line
        # Ensure indices are within bounds
        if x1 < 0 or x1 >= lake.shape[0] or y1 < 0 or y1 >= lake.shape[1]:
            continue  # Skip out-of-bounds points

        # Propagate the wave from the current point along the line
        for border_point in find_borders(lake):  # Iterate over all boundary points
            x2, y2 = border_point
            # Trace the line from (x1, y1) to (x2, y2) using Bresenham's algorithm
            rr, cc = line(x1, y1, x2, y2)  # note: (y, x) order
            line_path = list(zip(cc, rr))
            # line_path = bresenham(x1, y1, x2, y2)
            for (x, y) in line_path:
                # Ensure indices are within bounds
                if x < 0 or x >= lake.shape[0] or y < 0 or y >= lake.shape[1]:
                    continue  # Skip out-of-bounds points

                if lake[x, y] == 1:  # Continue propagation only in water area
                    propagated_lake[x, y] = 1
                else:
                    break  # Stop if we hit land

        # Count the number of 1s in the propagated_lake (which is the area covered)
        area_size = np.sum(propagated_lake)

        # If this is the largest area, update the best point and best area size
        if area_size > best_area_size:
            best_area_size = area_size
            best_point = shifted_point

    return best_point, best_area_size, propagated_lake


def compute_lake_angles_recursive(lake1, A, B, Angles_matrices, island_mask, threshold=20,
                                  recursion_depth=0):
    """
    Recursively computes angles for wave propagation, handling branches and sub-branches.
    """

    max_angles = []
    stack = [(lake1, A, B, Angles_matrices, threshold, recursion_depth, max_angles)]  # Initial state
    while stack:
        lakeCall, Acall, Bcall, Angles_matrices, threshold, recursion_depth, max_angles = stack.pop()
        MAX_DEPTH = 10
        print('this is ACall, bcall:  ', Acall, Bcall)
        print('recursion :', recursion_depth)
        if recursion_depth > MAX_DEPTH:  # To avoid infinite recursion, limit the recursion depth
            return Angles_matrices

        if recursion_depth > 0:
            previous_angles = Angles_matrices[recursion_depth - 1, :, :]
    ### for each recursion level compute each branches first find the branches from starting point or lines A or AB
        angles = np.zeros_like(lakeCall[0,:,:], dtype=float)
        new_lake = np.zeros((1, *lakeCall[0,:,:].shape), dtype=float)

        first_lake = 1
        An, Bn = [], []
        for Call in range(len(Acall)):
            if recursion_depth > 0:
                A = Acall[Call]
                # print(A)
                A = A[::-1]
                B = Bcall[Call]
                B = B[::-1]
            else:
                A = Acall[Call]
                B = Bcall[Call]

            lake = lakeCall[Call, :,:]

            structure = np.ones((3, 3))  # 3x3 structure for 8-connected neighborhood

            if recursion_depth == 0:
                propagated_lake = propagate_rays(lake, A, max_fill_distance=1)
            else:
                MidllePointX = (A[0] + B[0]) / 2
                MidllePointY = (A[1] + B[1]) / 2
                MiddlePoint = [int(MidllePointX), int(MidllePointY)]
                print('Mid: ', MiddlePoint, np.sum(lake))
                MiddlePoint = move_point_inside(lake, MiddlePoint)
                    #                     #

                propagated_lake = propagate_rays(lake, MiddlePoint, max_fill_distance=1)

                # Label unfilled areas
            if propagated_lake.ndim >= 3:
                print('squeeze')
                propagated_lake = propagated_lake.squeeze()


            remaining_areas, num_features = label((lake == 1) & (propagated_lake == 0), structure)
            # Find the coordinates of the boundary of the 100% region
            # boundary_points = np.argwhere(propagated_lake == 1)


            # if recursion_depth > MAX_DEPTH-1:  # To avoid infinite recursion, limit the recursion depth
            #     return Angles_matrices
            for area_label in range(1, num_features + 1):
                # Get mask for the current unfilled area
                area_mask = remaining_areas == area_label
                if np.sum(area_mask) < threshold:
                    propagated_lake[area_mask] = 1
                    continue


                AnTT, BnTT, area_masks, propagated_lakes = split_and_label_area_mask_by_islands(lake, area_mask, propagated_lake, A)


                for i in range(len(AnTT)):

                    AnT = AnTT[i]
                    BnT = BnTT[i]
                    area_maskT = area_masks[i]
                    area_points = np.argwhere(area_maskT)
                    area_points = area_points[:, [1, 0]]
                    if np.sum(area_maskT) == 0:
                        continue


                    area_border = binary_dilation(area_maskT) & (~area_maskT)
                    ## check if touches_island


                    if first_lake == 1:
                        new_lake[0,:,:] = area_maskT
                        first_lake += 1
                    else:
                        new_lake = np.append(new_lake, area_maskT[np.newaxis, :, :], axis=0)

                    if AnT is None or AnT.size == 0:  # For NumPy arrays
                        continue
                    An.append(AnT)
                    Bn.append(BnT)
                    # Compute the angle for each cell in this area
                    for C in area_points:
                        angle = compute_angle_for_cell(BnT, AnT, C)
                        angles[C[1], C[0]] = abs(angle)

                    if np.any(area_border & island_mask):
                        for C in area_points:
                            angle = compute_angle_for_cell(BnT, AnT, C)
                            if recursion_depth > 0:
                                if angles[C[1], C[0]] == 0 and previous_angles[C[1], C[0]] == 0:
                                    angles[C[1], C[0]] = abs(angle)
                                else:
                                    previous_ang = angles[C[1], C[0]]
                                    previous_angIt = previous_angles[C[1], C[0]]

                                    angles[C[1], C[0]] = abs(min(abs(angle), previous_ang) - previous_angIt)
                            else:
                                if angles[C[1], C[0]] == 0:
                                    angles[C[1], C[0]] = abs(angle)
                                else:
                                    previous_ang = angles[C[1], C[0]]
                                    angles[C[1], C[0]] = abs(abs(angle) - previous_ang)
                    else:
                        for C in area_points:
                            angle = compute_angle_for_cell(BnT, AnT, C)
                            angles[C[1], C[0]] = abs(angle)

                    if recursion_depth > 0:
                        rr, cc = line(AnT[0], AnT[1], BnT[0], BnT[1])  # note: (y, x) order
                        line_points = list(zip(rr, cc))

                        # line_points = bresenham(AnT[0], AnT[1], BnT[0], BnT[1])
                        ################### append matrices 1 by one so we can compute cumulative angles afterwards... ######
                        max_angle = float('-inf')  # Start with the smallest possible value

                        # Iterate over the line points and print the angles
                        for (x, y) in line_points:
                            # Ensure that the indices are within bounds
                            # print(x, y, angles.shape)
                            if 0 <= x < angles.shape[1] and 0 <= y < angles.shape[0]:
                                current_angle = previous_angles[y, x]
                                # print('here')

                                # Update the maximum angle if the current angle is larger
                                if current_angle > max_angle:
                                    max_angle = current_angle
                                    max_angle_point = (x, y)


                        area_maskT = area_maskT.astype(bool)
                        previous_angles[area_maskT] = max_angle

                        max_angles.append(max_angle)






        if recursion_depth == 0:
            Angles_matrices[0, :, :] = angles

        else:

            Angles_matrices[recursion_depth - 1, :, :] = previous_angles
            Angles_matrices = np.append(Angles_matrices, angles[np.newaxis, :, :], axis=0)



        if Call == len(Acall)-1:
            stack.append((new_lake, An, Bn, Angles_matrices, threshold, recursion_depth + 1, max_angles))


    return Angles_matrices
def extract_islands_from_lake(lake_mask):
    from scipy.ndimage import binary_fill_holes

    # Fill holes in lake -> gives you lake + islands
    filled = binary_fill_holes(lake_mask)

    # Islands are what was added to fill the holes
    island_mask = filled & ~lake_mask

    return island_mask





def compute_lake_angles(lake1, start_point, threshold=20):
    ## compute the angles from the starting point as soon as the wave front needs to turn.
    ## input : shape of the lake as a binary 1 for water, 0 for ground
    ##         coordinate of the starting point
    ## return a matrix with the lake binary shape, with 1 where there is no angles, 0 for grounds and the angles for each pixels.

    ## if min size == -1 filter all islands
    lake = erase_Islands(lake1, size_lim = 100000)

    # Propagate the rays and fill the lake area
    result = propagate_rays(lake, start_point, max_fill_distance=1)
    # # Step to fill isolated pixels before computing the angles
    # previous_angles = np.zeros_like(lake)
    Angles_matrices = np.zeros((1, *lake.shape), dtype=float)
    lakeN = np.expand_dims(lake, axis=0)
    island_mask = extract_islands_from_lake(lake)

    Angles_matrices = compute_lake_angles_recursive(lakeN, [start_point], [start_point], Angles_matrices, island_mask, threshold,
                                                    recursion_depth=0)

    stacked_angles = np.stack(Angles_matrices, axis=0)
    final_raster_stack = np.zeros_like(stacked_angles, dtype=float)
    if final_raster_stack.ndim > 2:
        for i in range(final_raster_stack.shape[0]):
            # Set NaN for non-water areas
            final_raster_stack[i, lake == 0] = np.nan

            # Set remaining angles for non-ground areas where result == 0 and lake == 1
            mask = (lake == 1) & (result == 0)
            final_raster_stack[i, mask] = stacked_angles[i, mask]
            # Set to 0 for areas fully propagated (where result == 1)
            final_raster_stack[i, result == 1] = 0

    return final_raster_stack

