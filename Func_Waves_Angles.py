import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import binary_erosion, distance_transform_edt, binary_dilation, center_of_mass
from math import atan2, degrees
from scipy.spatial.distance import cdist
from scipy.spatial import distance_matrix
import matplotlib
from skimage.measure import regionprops
from skimage.segmentation import flood
from scipy.ndimage import label as label
from skimage.measure import label as lab
from skimage.morphology import skeletonize
from scipy.ndimage import distance_transform_edt
from skimage.graph import route_through_array
from skimage.morphology import thin
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
    # print('touch:', touch_zone[touch_zone])
    # plt.figure(figsize=(8, 8))
    # plt.imshow(touch_zone)
    # # Create masks for where values are 1
    # maskA = area_mask == 1
    # maskB = propagated_lake == 1
    # print(maskA)
    # # Plot rasterA in red
    # plt.imshow(maskA, cmap='Reds', alpha=0.5)
    #
    # # Plot rasterB in blue
    # plt.imshow(maskB, cmap='Blues', alpha=0.5)
    #
    # # Optional: add title and grid
    # plt.title("Overlay of rasterA (red) and rasterB (blue)")
    # plt.axis('off')  # Hide axes

    # plt.show()
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
    # plt.figure()
    # # 6. Optional: visualize
    # plt.imshow(skeleton, cmap='gray')
    # plt.scatter(start_pt[1], start_pt[0])
    At = []
    Bt = []
    i = 0
    start_point = [start_pt[1], start_pt[0]]
    for a, b in endpoints:
        # plt.plot([a[0], b[0]], [a[1], b[1]], 'r-o')
        # plt.scatter(a[0], a[1],  color='r')
        # plt.scatter(b[0], b[1],  color='b')
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
        # plt.scatter(aT[0], aT[1], color='r')
        # plt.scatter(bT[0], bT[1], color='b')
    if At == []:
        return At, Bt, area_mask, propagated_lake
    AB_pairs = list(zip(At, Bt))
    # Sort pairs by distance from A to start_point
    AB_pairs.sort(key=lambda pair: np.linalg.norm(pair[0] - np.array(start_point)))

    # Unzip back to A and B
    AL, BL = zip(*AB_pairs)
    AL = list(AL)
    BL = list(BL)

    # plt.figure()
    # plt.imshow(area_mask)
    # # plt.title("Shared Lines and Endpoints")
    # plt.show()
    # New filtered lists
    AL_filtered = []
    BL_filtered = []

    area_maskN = area_mask.copy()
    area_maskTemp = area_mask.copy()
    # propagated_lakes = propagated_lake.copy()
    # propagated_lakes = propagated_lakes[np.newaxis, :, :]
    area_maskN = area_maskN[np.newaxis, :, :]
    j = 0
    # print(At, Bt)
    for A, B in zip(AL, BL):
        # print(A, B)
        MidllePointX = (A[0] + B[0]) / 2
        MidllePointY = (A[1] + B[1]) / 2
        MiddlePoint = [int(MidllePointX), int(MidllePointY)]
        # plt.figure()
        # plt.imshow(area_maskTemp)
        # plt.scatter(A[0], A[1])
        # plt.scatter(B[0], B[1])
        #
        # plt.scatter(MiddlePoint[0], MiddlePoint[1])
        # plt.show()
        if np.sum(area_maskTemp) == 0:
            print('skip Mask')

            continue
        if area_maskTemp[MiddlePoint[1], MiddlePoint[0]] == 0:
            # plt.figure()
            # plt.imshow(area_maskTemp)
            # plt.scatter(A[0], A[1])
            # plt.scatter(B[0], B[1])
            #
            # plt.scatter(MiddlePoint[0], MiddlePoint[1])
            # plt.show()
            try:
                MiddlePointProp = [MiddlePoint[1], MiddlePoint[0]]
                print(MiddlePointProp)
                new_y, new_x = find_nearest_non_zero(area_maskTemp, MiddlePointProp, radius=2)
                MiddlePointProp = [new_y, new_x]
                print(MiddlePointProp)

            except:
                print('other side of island')
                continue
        else:
            MiddlePointProp = [MiddlePoint[1], MiddlePoint[0]]
        print(MiddlePoint, MiddlePointProp)

        # Return None if no point was found
        # MiddlePoint = move_point_inside(area_maskTemp, MiddlePoint1)

        # plt.figure()
        # plt.imshow(area_maskTemp)
        # plt.scatter(MiddlePoint1[0], MiddlePoint1[1])
        # plt.scatter(MiddlePoint[0], MiddlePoint[1])
        ###################################################################################################

        ####   change area mask_temp to area mask to then check overlap ####
        propagated_lakeT = propagate_rays(area_mask, MiddlePointProp, max_fill_distance=1)
        # best_point, best_area_size, propagated_lakeT = propagate_rays_from_line(area_maskTemp, A, B)
        ###################################################################################################
        #
        # plt.figure()
        # plt.imshow(propagated_lakeT, alpha = 0.5, cmap='Greens')
        # plt.imshow(area_mask, alpha = 0.2, cmap='Blues')
        # if j > 0:
        #     plt.imshow(propagated_lakes[j-1,:,:], alpha = 0.5)
        # # plt.scatter(MiddlePoint1[0], MiddlePoint1[1])
        # plt.scatter(MiddlePoint[0], MiddlePoint[1])
        # plt.show()
        # print(np.sum(propagated_lakeT))
        # print(np.sum(area_maskTemp))
        # if np.sum(propagated_lakeT) < 10:
        #     continue
        AL_filtered.append(A)
        BL_filtered.append(B)


        if propagated_lake.ndim >= 3:
            print('squeeze')
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
    # if num_features == 0:
    #     print('num0', np.sum(propagated_lakes[0,:,:]))
    #     plt.figure()
    #     plt.imshow(area_maskTemp)
    #     plt.show()
    #     plt.figure()
    #     plt.imshow(propagated_lake)
    #     plt.show()
    #     print(len(propagated_lakes))
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
        # plt.figure()
        # plt.imshow(merged)
        # plt.title('mask_merged '+str(i))
        # plt.show()
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

def Merge_islands_erase_small_Islands(lake, min_size=100, is_island=None, max_islands=None):
    lake = lake.astype(np.uint8)
    land = (lake == 0)

    labeled_islands, num_islands = label(land)

    merged_mask = np.zeros_like(lake, dtype=bool)
    i = 0

    for island_label in range(1, num_islands + 1):
        island_mask = (labeled_islands == island_label)

        # 🔎 Use your custom island check
        if is_island is not None and not is_island(island_mask, lake):
            continue

        i += 1
        if max_islands is not None and i > max_islands:
            continue
        #
        # dilated = binary_dilation(island_mask, iterations=proximity)
        # merged_mask |= dilated
        merged_mask |=island_mask
    # Re-label after merging
    merged_labeled, merged_num = label(merged_mask)

    # Filter out small components
    final_mask = np.zeros_like(lake, dtype=bool)
    for label_id in range(1, merged_num + 1):
        component = (merged_labeled == label_id)
        if np.sum(component) >= min_size:
            final_mask |= component

    return final_mask.astype(np.uint8)


def erase_all_islands(lake, is_island=None):
    lake = lake.astype(np.uint8)

    # Identify land (0s)
    land = (lake == 0)

    # Label all connected land patches
    labeled_land, num_patches = label(land)

    # Copy lake to modify
    cleaned_lake = lake.copy()

    for label_id in range(1, num_patches + 1):
        island_mask = (labeled_land == label_id)

        # Use your custom is_island function
        if is_island is not None and is_island(island_mask, lake):
            cleaned_lake[island_mask] = 1  # Turn island into water

    return cleaned_lake

def split_area_mask_around_islands(lake, area_mask, start_point):
    """
    Iteratively split area_mask when islands are inside it and retain only the part closest to start_point.

    Parameters:
    - lake: 2D array (1 for water, 0 for land).
    - area_mask: 2D boolean array (area to process).
    - start_point: (x, y) tuple (column, row).

    Returns:
    - Cleaned area_mask with isolated parts removed.
    """
    # Convert start_point to (row, col)
    sy, sx = start_point[1], start_point[0]

    lake = lake.astype(np.uint8)
    water = (lake == 1)
    land = (lake == 0)

    # Step 1: Label all land blobs (potential islands)
    labeled_islands, num_islands = label(land)
    island_mask = labeled_islands == 169

    i = 0
    for island_label in range(1, num_islands + 1):
        island_mask = labeled_islands == island_label
        # ⛔ Skip if it's not an actual island
        if not is_island(island_mask, lake):
            continue
        i +=1
        if i > 2:
            continue
        # Check if island intersects the area_mask




        if np.any(island_mask & (area_mask+1)):
            print('dont skip')

            # Get island centroid
            props = regionprops(island_mask.astype(int))[0]
            cy, cx = props.centroid
            # print(cx, cy, sx, sy)
            # Get direction vector from start to island
            dy, dx = cy - sy, cx - sx

            # Decide cut direction
            if abs(dy) > abs(dx):  # Vertical direction => horizontal cut
                cut_row = int(np.max(np.argwhere(island_mask)[:, 0]))  # Bottom of island
                area_mask[cut_row:, :] = False
            else:  # Horizontal direction => vertical cut
                cut_col = int(np.max(np.argwhere(island_mask)[:, 1]))  # Right side of island
                area_mask[:, cut_col:] = False

            # Retain only connected part that includes start_point
            # if not area_mask[sy, sx]:
            #     print("Start point removed from area. Skipping this island.")
            #     continue  # Avoid removing everything

            retained_mask = flood(area_mask, (sy, sx))
            area_mask = retained_mask & area_mask  # update

    return area_mask
#
# def is_island(blob_mask, lake):
#     """
#     Determine if a labeled land blob is an island (i.e. fully enclosed by water).
#
#     Parameters:
#     - blob_mask: Boolean mask of the land blob.
#     - lake: Binary array (1 for water, 0 for land).
#
#     Returns:
#     - True if island, False otherwise.
#     """
#     h, w = lake.shape
#
#     # Get bounding box of the blob
#     ys, xs = np.where(blob_mask)
#     min_y, max_y = ys.min(), ys.max()
#     min_x, max_x = xs.min(), xs.max()
#
#     # If touching the edge, it's not an island
#     if min_y == 0 or max_y == h - 1 or min_x == 0 or max_x == w - 1:
#         return False
#
#     # Check surrounding ring: should all be water (1)
#     outer_ring = lake[min_y-1:max_y+2, min_x-1:max_x+2]
#     inner_area = blob_mask[min_y-1:max_y+2, min_x-1:max_x+2]
#
#     border = np.logical_and(~inner_area, outer_ring == 0)  # land outside
#     if np.any(border):
#         return False  # touches other land
#
#     return True


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
def find_boundary_points_along_direction(lake, start_point, target_point, area_mask, boundary_points, distance_threshold=5):
    """
    Finds boundary points that are along the direction from start_point to target_point
    and adjacent to the unfilled area.

    Parameters:
    - lake: 2D binary map (1=water, 0=land)
    - start_point: (x, y) tuple
    - target_point: (x, y) tuple defining the direction
    - area_mask: 2D mask (1=unfilled area)
    - boundary_points: list of (x, y) tuples of filled boundary
    - distance_threshold: max perpendicular distance to the line to be accepted

    Returns:
    - closest_point, farthest_point along the direction
    """
    boundary_points = np.array(boundary_points)

    # Convert start/target to numpy arrays in (y, x) = (row, col)
    p0 = np.array([start_point[1], start_point[0]])
    p1 = np.array([target_point[1], target_point[0]])

    # Direction vector
    v = p1 - p0
    v_norm = v / np.linalg.norm(v)

    # List for adjacent boundary points (that are close to the direction line)
    filtered_boundary_points = []

    for pt in boundary_points:
        y, x = pt
        neighbors = [
            (y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1),
            (y - 1, x - 1), (y - 1, x + 1), (y + 1, x - 1), (y + 1, x + 1)
        ]

        if any(0 <= ny < lake.shape[0] and 0 <= nx < lake.shape[1]
               and lake[ny, nx] == 1 and area_mask[ny, nx] == 1
               for ny, nx in neighbors):
            # Convert to (row, col) = (y, x)
            pt_vector = np.array([y, x]) - p0

            # Project pt_vector onto the direction vector v_norm
            proj_length = np.dot(pt_vector, v_norm)
            proj_point = p0 + proj_length * v_norm

            # Perpendicular distance from pt to the line
            perp_dist = np.linalg.norm(np.array([y, x]) - proj_point)

            if perp_dist <= distance_threshold:
                filtered_boundary_points.append((x, y))

    if not filtered_boundary_points:
        return None, None

    filtered_boundary_points = np.array(filtered_boundary_points)

    # Compute distances from start_point to the filtered points
    distances = cdist([[start_point[1], start_point[0]]], filtered_boundary_points[:, [1, 0]])[0]

    closest_idx = np.argmin(distances)
    farthest_idx = np.argmax(distances)

    return filtered_boundary_points[closest_idx], filtered_boundary_points[farthest_idx]

def find_boundary_points_for_area(lake, start_point, area_mask, boundary_points):
    """
    Finds the closest and farthest boundary points for the given area that are adjacent to the unfilled area.

    Parameters:
    - lake: The original binary map of the lake (1 for water, 0 for land).
    - start_point: Coordinates of the explosion start point.
    - area_mask: Mask representing the unfilled area for which we need to find boundary points.
    - boundary_points: List of points representing the boundary of the 100% filled region.

    Returns:
    - closest_point: The boundary point closest to the start point.
    - farthest_point: The boundary point farthest from the start point.
    """
    # Ensure that boundary_points is a numpy array for easier manipulation
    boundary_points = np.array(boundary_points)
    start_point = [start_point[1], start_point[0]]
    # List to hold adjacent boundary points
    adjacent_boundary_points = []

    # Iterate over the boundary points
    for boundary_point in boundary_points:
        # Reverse x and y for correct handling in the raster data (row, column)
        y, x = boundary_point  # boundary_point = (x, y) -> we swap to (y, x) for raster handling

        # Check the neighboring coordinates (4-connectivity + 8-connectivity)
        neighbors = [
            (y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1),  # Direct neighbors (4-connectivity)
            (y - 1, x - 1), (y - 1, x + 1), (y + 1, x - 1), (y + 1, x + 1)  # Diagonal neighbors (8-connectivity)
        ]

        # Check if any of the neighbors is part of the unfilled area (area_mask == 1) and adjacent to the filled area (lake == 1)
        if any(
                0 <= ny < lake.shape[0] and 0 <= nx < lake.shape[1] and lake[ny, nx] == 1 and area_mask[ny, nx] == 1
                for ny, nx in neighbors
        ):
            adjacent_boundary_points.append((x, y))  # Append in (x, y) order for consistency

    # If no adjacent boundary points were found, return None
    if not adjacent_boundary_points:
        return None, None

    # Convert adjacent_boundary_points to a numpy array for distance calculations
    adjacent_boundary_points = np.array(adjacent_boundary_points)

    # Compute distances from the start point to all adjacent boundary points
    dist_to_adjacent = cdist([start_point], adjacent_boundary_points)

    # Find the index of the closest and farthest boundary points
    closest_idx = np.argmin(dist_to_adjacent)
    farthest_idx = np.argmax(dist_to_adjacent)

    # Return the closest and farthest boundary points (in (x, y) order)
    return adjacent_boundary_points[closest_idx], adjacent_boundary_points[farthest_idx]


def find_boundary_points_for_area_from_line(lake, start_point, end_point, area_mask, boundary_points):
    """
    Finds the closest and farthest boundary points for the given area based on the line defined by start_point and end_point.

    Parameters:
    - lake: The original binary map of the lake (1 for water, 0 for land).
    - start_point: Starting coordinates of the line.
    - end_point: Ending coordinates of the line.
    - area_mask: Mask representing the unfilled area for which we need to find boundary points.
    - boundary_points: List of points representing the boundary of the 100% filled region.

    Returns:
    - closest_point: The boundary point closest to the line.
    - farthest_point: The boundary point farthest from the line.
    """
    # Ensure that boundary_points is a numpy array for easier manipulation

    boundary_points = np.array(boundary_points)

    # Generate the points along the line between start_point and end_point (e.g., Bresenham's algorithm)
    # line_points = bresenham(start_point[0], start_point[1], end_point[0], end_point[1])
    line_points = start_point
    # List to hold adjacent boundary points
    adjacent_boundary_points = []
    # print('start_point:', line_points)

    # Iterate over each point on the line
    for line_point in line_points:
        # y_start, x_start = line_point  # Reverse for raster handling if needed

        # Iterate over the boundary points
        for boundary_point in boundary_points:
            y, x = boundary_point  # boundary_point = (x, y) -> we swap to (y, x) for raster handling

            # Check the neighboring coordinates (4-connectivity + 8-connectivity)
            neighbors = [
                (y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1),  # Direct neighbors (4-connectivity)
                (y - 1, x - 1), (y - 1, x + 1), (y + 1, x - 1), (y + 1, x + 1)  # Diagonal neighbors (8-connectivity)
            ]

            # Check if any of the neighbors is part of the unfilled area (area_mask == 1) and adjacent to the filled area (lake == 1)
            if any(
                    0 <= ny < lake.shape[0] and 0 <= nx < lake.shape[1] and lake[ny, nx] == 1 and area_mask[ny, nx] == 1
                    for ny, nx in neighbors
            ):
                adjacent_boundary_points.append((x, y))  # Append in (x, y) order for consistency

    # If no adjacent boundary points were found, return None
    if not adjacent_boundary_points:
        return None, None

    # Convert adjacent_boundary_points to a numpy array for distance calculations
    adjacent_boundary_points = np.array(adjacent_boundary_points)

    # Compute distances from all points on the line to all adjacent boundary points

    dist_to_adjacent = cdist([[line_points[1],line_points[0]]] , adjacent_boundary_points)

    # Find the index of the closest and farthest boundary points
    closest_row, closest_col = np.unravel_index(np.argmin(dist_to_adjacent), dist_to_adjacent.shape)
    farthest_row, farthest_col = np.unravel_index(np.argmax(dist_to_adjacent), dist_to_adjacent.shape)
    print('All:', line_points, adjacent_boundary_points[closest_col], adjacent_boundary_points[farthest_col])
    # x = boundary_points[:, 0]
    # y = boundary_points[:, 1]
    # plt.figure()
    # plt.scatter(y,x,  s=1, color='k')
    # plt.scatter(start_point[1], start_point[0], color='b')
    # plt.scatter(adjacent_boundary_points[closest_col][0], adjacent_boundary_points[closest_col][1], color='g')
    # plt.scatter(adjacent_boundary_points[farthest_col][0], adjacent_boundary_points[farthest_col][1], color='orange')
    # plt.show()

    # Return the closest and farthest boundary points (in (x, y) order)
    return adjacent_boundary_points[closest_col], adjacent_boundary_points[farthest_col]


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

def merge_contiguous_masks(area_masks, dilation=True):
    n, h, w = area_masks.shape
    used = np.zeros(n, dtype=bool)
    merged_masks = []

    for i in range(n):
        if used[i]:
            continue

        # Start new merged mask
        group_mask = area_masks[i].copy()
        used[i] = True
        changed = True

        while changed:
            changed = False
            for j in range(n):
                if not used[j]:
                    # Check for overlap or touching using binary_dilation
                    overlap = (group_mask & area_masks[j]).any()
                    if not overlap and dilation:
                        overlap = (binary_dilation(group_mask) & area_masks[j]).any()

                    if overlap:
                        group_mask |= area_masks[j]
                        used[j] = True
                        changed = True

        merged_masks.append(group_mask)

    return np.array(merged_masks)

def compute_lake_angles_recursive(lake1, A, B, Angles_matrices, distanceI, island_mask, threshold=20,
                                  recursion_depth=0):
    """
    Recursively computes angles for wave propagation, handling branches and sub-branches.
    """
    # stack = [(lake, A, B, propagated_lake, 0)]  # Initial state
    # while stack:
    #     lake, A, B, propagated_lake, recursion_depth = stack.pop()
    distance = np.round(distanceI, decimals=1)
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
            # plt.figure()
            # plt.imshow(lake)
            # plt.show()


            # print('A, B:  ',A, B)
            structure = np.ones((3, 3))  # 3x3 structure for 8-connected neighborhood
            # plt.figure()
            # plt.imshow(lake)
            # plt.title('lake')
            # plt.scatter(A[1], A[0])
            # plt.scatter(B[1], B[0])
            # plt.show()

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
                # print('Mid2')

                # Label unfilled areas
            if propagated_lake.ndim >= 3:
                print('squeeze')
                propagated_lake = propagated_lake.squeeze()
            # plt.figure()
            # plt.imshow(propagated_lake)
            # plt.title('ProLake')
            # plt.show()

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

                # plt.figure()
                # plt.imshow(area_mask)
                # plt.title('before')
                # plt.figure()
                # plt.imshow(propagated_lake)
                # plt.title('beforeP')
                # plt.show()
                AnTT, BnTT, area_masks, propagated_lakes = split_and_label_area_mask_by_islands(lake, area_mask, propagated_lake, A)
                print(AnTT)


                # if (recursion_depth == 0) and (area_label == 1):
                #     for i in range(len(AnTT)):
                #         plt.figure()
                #         plt.imshow(propagated_lakes[i, :, :])
                #         plt.scatter(AnTT[i][0], AnTT[i][1])
                #         plt.show()
                # if recursion_depth == 1:  # To avoid infinite recursion, limit the recursion depth
                #     break
                # print(AnTT)
                for i in range(len(AnTT)):

                    AnT = AnTT[i]
                    BnT = BnTT[i]
                    # print(AnT)
                    area_maskT = area_masks[i]
                    propagated_lakeT = propagated_lakes[i]
                    area_points = np.argwhere(area_maskT)
                    area_points = area_points[:, [1, 0]]
                    if np.sum(area_maskT) == 0:
                        print('mask0:', AnT)
                        continue
                    # if len(area_points) < threshold:
                    #     if recursion_depth == 0:
                    #         propagated_lakeT[area_maskT] = 1  # Mark the entire area as 100%
                    #     else:
                    #
                    #         propagated_lakeT[area_maskT] = previous_angles[area_maskT]  # Mark the entire area as 100%

                    # else:
                        # lakeSmoothed = gaussian_filter(lake, sigma=1)
                        # area_mask = split_area_mask_around_islands(lake, area_mask, A)
                        # AnT, BnT = split_and_label_area_mask_by_islands(lake, area_mask, propagated_lake, A)

                        # print(area_mask)
                    # plt.figure()
                    # plt.imshow(area_mask)
                    # plt.title('areasMask')
                    # plt.show()

                        # AnT, BnT = find_boundary_points_for_area_from_line(lake, A, B, area_mask, boundary_points)
                        # AnT, BnT = find_boundary_points_for_area(lake, A, area_mask, boundary_points)

                        # plt.figure()
                        # plt.imshow(propagated_lake)
                        # plt.scatter(A[1], A[0])

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
                    # print('An, Bn', An, Bn)
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

                    # if np.any(area_border & island_mask):
                    #     masked_distances = np.where(area_maskT, distance, np.nan)
                    #     unique_distances = np.unique(masked_distances[~np.isnan(masked_distances)])
                    #     for d in unique_distances:
                    #         mask = (distance == d) & area_maskT
                    #         angle_values = angles[mask.astype('bool')]
                    #
                    #         if len(np.unique(angle_values)) > 1:
                    #             min_angle = np.min(angle_values)
                    #             angles[mask] = min_angle
                    # if np.any(area_border & island_mask):
                    #     # Flatten everything
                    #     area_maskT_flat = area_maskT.ravel()
                    #     # distance_flat = np.round(distance.ravel(), 0)
                    #     distance_flat = distance.ravel() // 10 * 10
                    #     angles_flat = angles.ravel()
                    #
                    #     # Indices where area_maskT is True
                    #     valid_indices = np.flatnonzero(area_maskT_flat)
                    #
                    #     # Distances and angles inside the mask
                    #     valid_distances = distance_flat[valid_indices]
                    #     valid_angles = angles_flat[valid_indices]
                    #
                    #     # Group distances
                    #     unique_distances, inverse_indices = np.unique(valid_distances, return_inverse=True)
                    #
                    #     # Work on a copy of the full flattened angle array
                    #     updated_angles_flat = angles_flat.copy()
                    #
                    #     for i, d in enumerate(unique_distances):
                    #         # Mask for all indices where distance == d
                    #         group_mask = (inverse_indices == i)
                    #         group_indices = valid_indices[group_mask]  # Correct size now
                    #
                    #         group_angles = angles_flat[group_indices]
                    #
                    #         if len(np.unique(group_angles)) > 1:
                    #             min_angle = np.min(group_angles)
                    #             updated_angles_flat[group_indices] = min_angle
                    #
                    #     # Reshape back to 2D
                    #     angles[:, :] = updated_angles_flat.reshape(angles.shape)
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


                        # print(max_angle)
                        # plt.figure()
                        # plt.imshow(previous_angles)
                        # plt.plot(rr, cc, color='red')
                        # plt.title('before')
                        # plt.figure()
                        # plt.imshow(area_maskT)
                        # plt.title('mask')
                        area_maskT = area_maskT.astype(bool)
                        previous_angles[area_maskT] = max_angle

                        max_angles.append(max_angle)
                        # plt.figure()
                        # plt.imshow(area_maskT, alpha=0.5)
                        # plt.imshow(previous_angles, alpha=0.5, cmap='jet')
                        # plt.figure()
                        # plt.imshow(previous_angles)
                        # plt.title('after : '+ str(max_angle))

                        # plt.show()





        if recursion_depth == 0:
            Angles_matrices[0, :, :] = angles

        else:
            # plt.figure()
            # plt.imshow(previous_angles)
            # plt.show()
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





def compute_lake_angles(lake1, start_point, distance, threshold=20):
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

    Angles_matrices = compute_lake_angles_recursive(lakeN, [start_point], [start_point], Angles_matrices,distance,island_mask, threshold,
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

