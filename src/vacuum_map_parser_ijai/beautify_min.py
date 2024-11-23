# STRAIGHT REWRITE FROM MI HOME PLUGIN`S CODE
import base64

class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self) -> str:
        return f'Point({self.x},{self.y})'


class BeautifyMap:
    def __init__(self):
        self.map = []
        self.tRect = {
            "x": 0,
            "y": 0,
            "width": 0,
            "height": 0
        }
        self.x_min = 0
        self.x_max = 0
        self.y_min = 0
        self.y_max = 0
        self.resolution = 0
        self.size_x = 0
        self.size_y = 0

    def setSize(self, x_min1, x_max1, y_min1, y_max1, resolution1):
        self.map = []
        self.size_x = 0
        self.size_y = 0
        self.x_min = int(resolution1 * round(x_min1 / resolution1))
        self.y_min = int(resolution1 * round(y_min1 / resolution1))
        self.x_max = int(resolution1 * round(x_max1 / resolution1))
        self.y_max = int(resolution1 * round(y_max1 / resolution1))
        self.size_x = int(round((x_max1 - x_min1) / resolution1))
        self.size_y = int(round((y_max1 - y_min1) / resolution1))
        self.resolution = resolution1

    def setMap(self, tMapStruct):
        self.setSize(tMapStruct.mapHead.minX, tMapStruct.mapHead.maxX,
                     tMapStruct.mapHead.minY, tMapStruct.mapHead.maxY, tMapStruct.mapHead.resolution)

        mapdata_arr = tMapStruct.mapData.mapData

        tempArray = [0] * len(mapdata_arr)

        for i in range(len(mapdata_arr)):
            if (mapdata_arr[i] > 127):
                tempArray[i] = -128
            else:
                tempArray[i] = mapdata_arr[i]

        self.map = tempArray

    def getMap(self):
        return self.map

    def transform(self):
        black_boundary = []
        non_boundary_noise = []
        all_region = []
        self.tRect = self.findRoiMap(self.tRect)
        self.expandBlackRect(4, 4, self.map[0], self.tRect)
        self.expandWhiteRect(4, 4, self.map[0], self.tRect)
        self.refineBoundary(0, 10, self.tRect)
        non_boundary_noise = self.eliminateNonBoundaryNoise(
            non_boundary_noise, self.tRect, 127, -128, 0)
        self.expandSingleConvexBoundary(50, -128, 4, 4, self.tRect)
        black_boundary = self.removeIndependentRegion(
            all_region, black_boundary, 3, self.tRect)
        non_boundary_noise = self.fillNonBoundaryNoise2(
            non_boundary_noise, self.tRect)
        self.refineBoundary(0, 10, self.tRect)
        self.map = self.fillBlackComponent(self.map, black_boundary, -128)

    def findRoiMap(self, rect):
        top_bound = self.size_x
        bottom_bound = 0
        left_bound = self.size_y
        right_bound = 0
        for idx in range(self.size_x):
            for idy in range(self.size_y):
                if (self.map[idy * self.size_x + idx] != 0):
                    if (left_bound > idy - 10):
                        if (idx - 10 >= 0):
                            left_bound = idy - 10
                        else:
                            left_bound = 0
                    else:
                        break
        for _idx in range(self.size_x):
            for _idy in range(self.size_y - 1, 0, -1):
                if (self.map[_idy * self.size_x + _idx] != 0):
                    if (right_bound < _idy + 10):
                        if (_idy + 10 < self.size_y):
                            right_bound = _idy + 10
                        else:
                            right_bound = self.size_y - 1
                    else:
                        break
        for _idy2 in range(self.size_y):
            for _idx2 in range(self.size_x):
                if (self.map[_idy2 * self.size_x + _idx2] != 0):
                    if (top_bound > _idx2 - 10):
                        if (_idx2 - 10 >= 0):
                            top_bound = _idx2 - 10
                        else:
                            top_bound = 0
                    else:
                        break
        for _idy3 in range(self.size_y):
            for _idx3 in range(self.size_x - 1, 0, -1):
                if (self.map[_idy3 * self.size_x + _idx3] != 0):
                    if (bottom_bound < _idx3 + 10):
                        if (_idx3 + 10 < self.size_x):
                            bottom_bound = _idx3 + 10
                        else:
                            bottom_bound = self.size_x - 1
                    else:
                        break
        width = right_bound - left_bound + 1
        height = bottom_bound - top_bound + 1
        if (width > 0 and height > 0 and width < self.size_y and height < self.size_x):
            rect["x"] = top_bound
            rect["y"] = left_bound
            rect["width"] = width
            rect["height"] = height
            return rect
        else:
            rect["x"] = 0
            rect["y"] = 0
            rect["width"] = self.size_y
            rect["height"] = self.size_x
        return rect

    def expandBlackRect(self, kernel_size_x, kernel_size_y, threshold, rect):
        il, ir, jl, jr = (None, None, None, None)

        if (kernel_size_x % 2 == 1):
            ir = kernel_size_x - 1 >> 1
            il = -ir
        else:
            ir = kernel_size_x >> 1
            il = 1 - ir

        if (kernel_size_y % 2 == 1):
            jr = kernel_size_y - 1 >> 1
            jl = -jr
        else:
            jr = kernel_size_y >> 1
            jl = 1 - jr

        dst = []

        for i in range(self.size_y):
            for j in range(self.size_x):
                dst.append(127)

        for _i in range(rect["y"], rect["y"] + rect["width"]):
            for _j in range(rect["x"], rect["x"] + rect["height"]):
                if (self.map[_i * self.size_x + _j] < threshold):
                    for di in range(il, ir + 1):
                        for dj in range(jl, jr + 1):
                            if (_i + di < 0 or _i + di >= rect["y"] + rect["width"] or _j + dj < 0 or _j + dj >= rect["x"] + rect["height"]):
                                continue

                            if (dst[(_i + di) * self.size_x + _j + dj] > self.map[_i * self.size_x + _j]):
                                dst[(_i + di) * self.size_x + _j +
                                    dj] = self.map[_i * self.size_x + _j]

        for _i2 in range(self.size_y):
            for _j2 in range(self.size_x):
                if (dst[_i2 * self.size_x + _j2] == 127):
                    dst[_i2 * self.size_x + _j2] = self.map[_i2 * self.size_x + _j2]

        self.map = dst

    def expandWhiteRect(self, kernel_size_x, kernel_size_y, threshold, rect):
        il, ir, jl, jr = (None, None, None, None)

        if (kernel_size_x % 2 == 1):
            ir = kernel_size_x - 1 >> 1
            il = -ir
        else:
            ir = kernel_size_x >> 1
            il = 1 - ir

        if (kernel_size_y % 2 == 1):
            jr = kernel_size_y - 1 >> 1
            jl = -jr
        else:
            jr = kernel_size_y >> 1
            jl = 1 - jr

        dst = []

        for i in range(self.size_y):
            for j in range(self.size_x):
                dst.append(-128)

        for _i3 in range(rect["y"], rect["y"] + rect["width"]):
            for _j3 in range(rect["x"], rect["x"] + rect["height"]):
                if (self.map[_i3 * self.size_x + _j3] > threshold):
                    for di in range(il, ir + 1):
                        for dj in range(jl, jr + 1):
                            if (_i3 + di < 0 or _i3 + di >= rect["y"] + rect["width"] or _j3 + dj < 0 or _j3 + dj >= rect["x"] + rect["height"]):
                                continue

                            if (dst[(_i3 + di) * self.size_x + _j3 + dj] < self.map[_i3 * self.size_x + _j3] and self.map[(_i3 + di) * self.size_x + _j3 + dj] < threshold):
                                dst[(_i3 + di) * self.size_x + _j3 +
                                    dj] = self.map[_i3 * self.size_x + _j3]

        for _i4 in range(self.size_y):
            for _j4 in range(self.size_x):
                if (dst[_i4 * self.size_x + _j4] == -128):
                    dst[_i4 * self.size_x + _j4] = self.map[_i4 * self.size_x + _j4]

        self.map = dst

    def refineBoundary(self, threshold_black, threshold_white, rect):
        Qx = []
        Qy = []
        hasWhiteNeighbor = None

        for i in range(rect["y"], rect["y"] + rect["width"]):
            for j in range(rect["x"], rect["x"] + rect["height"]):
                if (self.map[i * self.size_x + j] < threshold_black):
                    hasWhiteNeighbor = False

                    for di in range(-1, 2):
                        for dj in range(-1, 2):
                            if (i + di < 0 or i + di >= rect["y"] + rect["width"] or j + dj < 0 or j + dj >= rect["x"] + rect["height"]):
                                continue

                            if (self.map[(i + di) * self.size_x + j + dj] > threshold_white):
                                hasWhiteNeighbor = True

                    if (not hasWhiteNeighbor):
                        Qx.append(i)
                        Qy.append(j)

        for _i5 in range(len(Qx)):
            self.map[Qx[_i5] * self.size_x + Qy[_i5]] = 0

    def eliminateNonBoundaryNoise(self, nonBoundaryNoise, rect, noise_color, border_color, outer_border_color):
        tempnonBoundaryNoise = nonBoundaryNoise

        for i in range(rect["y"], rect["y"] + rect["width"]):
            for j in range(rect["x"], rect["x"] + rect["height"]):
                if (self.map[i * self.size_x + j] == border_color):
                    if (i - 1 < 0 or i + 1 >= rect["y"] + rect["width"] or j - 1 < 0 or j + 1 >= rect["x"] + rect["height"]):
                        continue

                    if (self.map[(i - 1) * self.size_x + j] != outer_border_color and self.map[(i + 1) * self.size_x + j] != outer_border_color and self.map[i * self.size_x + j - 1] != outer_border_color and self.map[i * self.size_x + j + 1] != outer_border_color and self.map[(i - 1) * self.size_x + j - 1] != outer_border_color and self.map[(i - 1) * self.size_x + j + 1] != outer_border_color and self.map[(i + 1) * self.size_x + j - 1] != outer_border_color and self.map[(i + 1) * self.size_x + j + 1] != outer_border_color):
                        self.map[i * self.size_x + j] = noise_color
                        tempnonBoundaryNoise.append(Point(i, j))

        return tempnonBoundaryNoise

    def expandSingleConvexBoundary(self, external_corner_value, fill_value, valid_length, times, rect):
        contour = []
        contour_map = self.map
        result = self.extractExternalContoursNewStrategy(
            contour_map, contour, rect)

        for i in range(times):
            inner_corners = []
            extract_corners = []
            fill_edges = []
            inner_corner_value = external_corner_value + 5
            four_neighbourhood = [[-1, 0], [1, 0], [0, -1], [0, 1]]
            delete_point = []
            contour_map = result["temp_map"]
            contour = result["contour"]
            corner_map = self.map
            result1 = self.extractCorners(
                corner_map, extract_corners, inner_corners, contour, external_corner_value, inner_corner_value, rect)
            result2 = {
                "delete_point": [],
                "fill_edges": []
            }
            for it in result1["extract_corners"]:
                is_valid_length = False
                p_oint = it

                for k in range(4):
                    temp_idy = p_oint.x + four_neighbourhood[k][0]
                    temp_idx = p_oint.y + four_neighbourhood[k][1]

                    if (temp_idy < rect["y"] or temp_idy >= rect["y"] + rect["width"] or temp_idx < rect["x"] or temp_idx >= rect["x"] + rect["height"]):
                        continue

                    if (result1["corner_map"][temp_idy * self.size_x + temp_idx] == inner_corner_value):
                        near_inner_p_oint = Point(temp_idy, temp_idx)
                        is_valid_length = self.statisticalLineLength(
                            result1["corner_map"], near_inner_p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect)
                        break

                    if (k == 3):
                        is_valid_length = True

                result2 = self.fourNeighbourhoodSearchForExtractCorners(
                    result1["corner_map"], p_oint, fill_edges, delete_point, external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect)
                fill_edges = result2["fill_edges"]
            # array = self.updateContour(contour, result2["delete_point"])
            array = contour
            tempContour = self.fillEdges(
                self.map, result1["extract_corners"], array, fill_edges, fill_value)
            contour = tempContour
            delete_point = []
            extract_corners = []
            fill_edges = []

        contour = []

    def extractExternalContoursNewStrategy(self, temp_map, contour, rect):
        gray_region = []
        result = self.findGrayConnectComponent(temp_map, gray_region, rect)
        gray_region = result["gray_region"]
        temp_map = result["temp_map"]
        result1 = self.findExternalContoursNewStrategy(
            temp_map, gray_region, contour, rect)
        return {
            "temp_map": result1["temp_map"],
            "contour": result1["contour"]
        }

    def findGrayConnectComponent(self, temp_map, gray_region, rect):
        four_neighbourhood = [[-1, 0], [0, 1], [1, 0], [0, -1]]
        findOnePoint = False

        for idy in range(rect["y"], rect["y"] + rect["width"]):
            for idx in range(rect["x"], rect["x"] + rect["height"]):
                if (temp_map[idy * self.size_x + idx] == 0):
                    findOnePoint = True
                    p_oint_for_search = []
                    gray_region.append(Point(idy, idx))
                    p_oint_for_search.append(Point(idy, idx))
                    temp_map[idy * self.size_x + idx] = 30

                    while (len(p_oint_for_search) > 0):
                        seed = p_oint_for_search[0]
                        p_oint_for_search.pop(0)

                        for k in range(4):
                            temp_idy = seed.x + four_neighbourhood[k][0]
                            temp_idx = seed.y + four_neighbourhood[k][1]

                            if (temp_idy < rect["y"] or temp_idy >= rect["y"] + rect["width"] or temp_idx < rect["x"] or temp_idx >= rect["x"] + rect["height"]):
                                continue

                            if (temp_map[temp_idy * self.size_x + temp_idx] == 0):
                                temp_map[temp_idy *
                                         self.size_x + temp_idx] = 30
                                p_oint_for_search.append(
                                    Point(temp_idy, temp_idx))
                                gray_region.append(Point(temp_idy, temp_idx))

                if (findOnePoint):
                    break

            if (findOnePoint):
                findOnePoint = False
                break

        return {
            "gray_region": gray_region,
            "temp_map": temp_map
        }

    def findExternalContoursNewStrategy(self, temp_map, gray_region, contour, rect):
        eight_neighbourhood = [[-1, 0], [1, 0], [0, -1],
                               [0, 1], [-1, 1], [1, 1], [1, -1], [-1, -1]]

        for i in range(len(gray_region)):
            for k in range(8):
                temp_idy = gray_region[i].x + eight_neighbourhood[k][0]
                temp_idx = gray_region[i].y + eight_neighbourhood[k][1]

                if (temp_idy < rect["y"] or temp_idy >= rect["y"] + rect["width"] or temp_idx < rect["x"] or temp_idx >= rect["x"] + rect["height"]):
                    continue

                if (temp_map[temp_idy * self.size_x + temp_idx] == -128):
                    temp_map[temp_idy * self.size_x + temp_idx] = 40
                    contour.append(Point(temp_idy, temp_idx))

        return {
            "temp_map": temp_map,
            "contour": contour
        }

    def extractCorners(self, corner_map, extract_corner, inner_corner, contour, external_corner_value, inner_corner_value, rect):
        four_neighbourhood = [[-1, 0], [0, 1], [1, 0], [0, -1]]

        for i in range(len(contour)):
            black_count = 0
            white_count = 0
            gray_count = 0

            for k in range(4):
                temp_idy = contour[i].x + four_neighbourhood[k][0]
                temp_idx = contour[i].y + four_neighbourhood[k][1]

                if (temp_idy < rect["y"] or temp_idy >= rect["y"] + rect["width"] or temp_idx < rect["x"] or temp_idx >= rect["x"] + rect["height"]):
                    continue

                if (self.map[temp_idy * self.size_x + temp_idx] == -128):
                    black_count += 1
                elif (self.map[temp_idy * self.size_x + temp_idx] == 0):
                    gray_count += 1
                elif (self.map[temp_idy * self.size_x + temp_idx] == 127):
                    white_count += 1

                if (gray_count == 2 and black_count == 2):
                    extract_corner.append(Point(contour[i].x, contour[i].y))
                    corner_map[contour[i].x * self.size_x +
                               contour[i].y] = external_corner_value
                elif (white_count == 2 and black_count == 2):
                    corner_map[contour[i].x * self.size_x +
                               contour[i].y] = inner_corner_value

        return {
            "corner_map": corner_map,
            "extract_corners": extract_corner
        }

    def statisticalLineLength(self, temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect):
        if (self.upSearchStatisticalLineLength(temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect)):
            return True
        elif (self.downSearchStatisticalLineLength(temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect)):
            return True
        elif (self.leftSearchStatisticalLineLength(temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect)):
            return True
        elif (self.rightSearchStatisticalLineLength(temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect)):
            return True

        return False

    def upSearchStatisticalLineLength(self, temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect):
        if (p_oint.x + 1 < rect["y"] + rect["width"] and self.map[(p_oint.x + 1) * self.size_x + p_oint.y] == 127):
            idy = p_oint.x + 1
            idx = p_oint.y
            line = []
            line.append(Point(idy, idx))

            for j in range(idy, rect["y"] + rect["width"]):
                if (temp_map[j * self.size_x + idx] == 127):
                    black_count = 0
                    left_and_right_neighbourhood = [[0, -1], [0, 1]]

                    for k in range(2):
                        tmp_idy = j + left_and_right_neighbourhood[k][0]
                        tmp_idx = idx + left_and_right_neighbourhood[k][1]

                        if (tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"] or tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(j, idx))
                    else:
                        break
                else:
                    break

            if (len(line) > valid_length):
                return True
            else:
                return False

        return False

    def downSearchStatisticalLineLength(self, temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect):
        if (p_oint.x - 1 > rect["y"] and self.map[(p_oint.x - 1) * self.size_x + p_oint.y] == 127):
            idy = p_oint.x - 1
            idx = p_oint.y
            line = []
            line.append(Point(idy, idx))

            for j in range(idy, rect["y"], -1):
                if (temp_map[j * self.size_x + idx] == 127):
                    black_count = 0
                    left_and_right_neighbourhood = [[0, -1], [0, 1]]

                    for k in range(2):
                        tmp_idy = j + left_and_right_neighbourhood[k][0]
                        tmp_idx = idx + left_and_right_neighbourhood[k][1]

                        if (tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"] or tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(j, idx))
                    else:
                        break
                else:
                    break

            if (len(line) > valid_length):
                return True
            else:
                return False
        return False

    def leftSearchStatisticalLineLength(self, temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect):
        if (p_oint.y + 1 < rect["x"] + rect["height"] and self.map[p_oint.x * self.size_x + p_oint.y + 1] == 127):
            idy = p_oint.x
            idx = p_oint.y + 1
            line = []
            line.append(Point(idy, idx))

            for j in range(idx, rect["x"] + rect["height"]):
                if (temp_map[idy * self.size_x + j] == 127):
                    black_count = 0
                    up_and_down_neighbourhood = [[-1, 0], [1, 0]]

                    for k in range(2):
                        tmp_idy = idy + up_and_down_neighbourhood[k][0]
                        tmp_idx = j + up_and_down_neighbourhood[k][1]

                        if (tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"] or tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(idy, j))
                    else:
                        break
                else:
                    break

            if (len(line) > valid_length):
                return True
            else:
                return False
        return False

    def rightSearchStatisticalLineLength(self, temp_map, p_oint, external_corner_value, inner_corner_value, fill_value, valid_length, rect):
        if (p_oint.y - 1 > rect["x"] and self.map[p_oint.x * self.size_x + p_oint.y - 1] == 127):
            idy = p_oint.x
            idx = p_oint.y - 1
            line = []
            line.append(Point(idy, idx))

            for j in range(idx, rect["x"], -1):
                if (temp_map[idy * self.size_x + j] == 127):
                    black_count = 0
                    up_and_down_neighbourhood = [[-1, 0], [1, 0]]

                    for k in range(2):
                        tmp_idy = idy + up_and_down_neighbourhood[k][0]
                        tmp_idx = j + up_and_down_neighbourhood[k][1]

                        if (tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"] or tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(idy, j))
                    else:
                        break
                else:
                    break

            if (len(line) > valid_length):
                return True
            else:
                return False

        return False

    def fourNeighbourhoodSearchForExtractCorners(self, temp_map, p_oint, fill_edges, delete_point, external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect):
        result1 = self.upSearchForExtractCorners(temp_map, p_oint, fill_edges, delete_point,
                                                 external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect)
        result2 = self.downSearchForExtractCorners(
            temp_map, p_oint, result1["fill_edges"], result1["delete_point"], external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect)
        result3 = self.leftSearchForExtractCorners(
            temp_map, p_oint, result2["fill_edges"], result2["delete_point"], external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect)
        result4 = self.rightSearchForExtractCorners(
            temp_map, p_oint, result3["fill_edges"], result3["delete_point"], external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect)
        return result4

    def upSearchForExtractCorners(self, temp_map, p_oint, fill_edges, delete_point, external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect):
        if (p_oint.x + 1 < rect["y"] + rect["width"] and self.map[(p_oint.x + 1) * self.size_x + p_oint.y] == 0):
            idy = p_oint.x + 1
            idx = p_oint.y
            line = []
            line.append(Point(idy, idx))

            for j in range(idy, rect["y"] + rect["width"]):
                if (temp_map[j * self.size_x + idx] == 0):
                    black_count = 0
                    left_and_right_neighbourhood = [[0, -1], [0, 1]]

                    for k in range(2):
                        tmp_idy = j + left_and_right_neighbourhood[k][0]
                        tmp_idx = idx + left_and_right_neighbourhood[k][1]

                        if (tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"] or tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(j, idx))
                    else:
                        break
                else:
                    break

            if (is_valid_length and len(line) > 1):
                line.append(p_oint)
                fill_edges.append(line)

                for i in range(len(line)):
                    _left_and_right_neighbourhood = [[0, -1], [0, 1]]

                    for _k in range(2):
                        _tmp_idy = line[i].x + \
                            _left_and_right_neighbourhood[_k][0]

                        _tmp_idx = line[i].y + \
                            _left_and_right_neighbourhood[_k][1]

                        if (_tmp_idx < rect["x"] or _tmp_idx >= rect["x"] + rect["height"] or _tmp_idy < rect["y"] or _tmp_idy >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[_tmp_idy * self.size_x + _tmp_idx] == -128 or temp_map[_tmp_idy * self.size_x + _tmp_idx] == inner_corner_value):
                            self.map[_tmp_idy * self.size_x + _tmp_idx] = 127
                            delete_point.append(Point(_tmp_idy, _tmp_idx))
            elif (len(line) > valid_length):
                line.append(p_oint)
                fill_edges.append(line)

                for _i6 in range(len(line)):
                    _left_and_right_neighbourhood2 = [[0, -1], [0, 1]]

                    for _k2 in range(2):
                        _tmp_idy2 = line[_i6].x + \
                            _left_and_right_neighbourhood2[_k2][0]

                        _tmp_idx2 = line[_i6].y + \
                            _left_and_right_neighbourhood2[_k2][1]

                        if (_tmp_idx2 < rect["x"] or _tmp_idx2 >= rect["x"] + rect["height"] or _tmp_idy2 < rect["y"] or _tmp_idy2 >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[_tmp_idy2 * self.size_x + _tmp_idx2] == -128 or temp_map[_tmp_idy2 * self.size_x + _tmp_idx2] == inner_corner_value):
                            self.map[_tmp_idy2 * self.size_x + _tmp_idx2] = 127
                            delete_point.append(Point(_tmp_idy2, _tmp_idx2))
            else:
                line = []

        return {
            "delete_point": delete_point,
            "fill_edges": fill_edges
        }

    def downSearchForExtractCorners(self, temp_map, p_oint, fill_edges, delete_point, external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect):
        if (p_oint.x - 1 > rect["y"] and self.map[(p_oint.x - 1) * self.size_x + p_oint.y] == 0):
            idy = p_oint.x - 1
            idx = p_oint.y
            line = []
            line.append(Point(idy, idx))

            for j in range(idy, rect["y"], -1):
                if (temp_map[j * self.size_x + idx] == 0):
                    black_count = 0
                    left_and_right_neighbourhood = [[0, -1], [0, 1]]

                    for k in range(2):
                        tmp_idy = j + left_and_right_neighbourhood[k][0]
                        tmp_idx = idx + left_and_right_neighbourhood[k][1]

                        if (tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"] or tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(j, idx))
                    else:
                        break
                else:
                    break

            if (is_valid_length and len(line) > 1):
                line.append(p_oint)
                fill_edges.append(line)

                for i in range(len(line)):
                    _left_and_right_neighbourhood3 = [[0, -1], [0, 1]]

                    for _k3 in range(2):
                        _tmp_idy3 = line[i].x + \
                            _left_and_right_neighbourhood3[_k3][0]

                        _tmp_idx3 = line[i].y + \
                            _left_and_right_neighbourhood3[_k3][1]

                        if (_tmp_idy3 < rect["y"] or _tmp_idy3 >= rect["y"] + rect["width"] or _tmp_idx3 < rect["x"] or _tmp_idx3 >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[_tmp_idy3 * self.size_x + _tmp_idx3] == -128 or temp_map[_tmp_idy3 * self.size_x + _tmp_idx3] == inner_corner_value):
                            self.map[_tmp_idy3 * self.size_x + _tmp_idx3] = 127
                            delete_point.append(Point(_tmp_idy3, _tmp_idx3))
            elif (len(line) > valid_length):
                line.append(p_oint)
                fill_edges.append(line)

                for _i7 in range(len(line)):
                    _left_and_right_neighbourhood4 = [[0, -1], [0, 1]]

                    for _k4 in range(2):
                        _tmp_idy4 = line[_i7].x + \
                            _left_and_right_neighbourhood4[_k4][0]

                        _tmp_idx4 = line[_i7].y + \
                            _left_and_right_neighbourhood4[_k4][1]

                        if (_tmp_idy4 < rect["y"] or _tmp_idy4 >= rect["y"] + rect["width"] or _tmp_idx4 < rect["x"] or _tmp_idx4 >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[_tmp_idy4 * self.size_x + _tmp_idx4] == -128 or temp_map[_tmp_idy4 * self.size_x + _tmp_idx4] == inner_corner_value):
                            self.map[_tmp_idy4 * self.size_x + _tmp_idx4] = 127
                            delete_point.append(Point(_tmp_idy4, _tmp_idx4))
            else:
                line = []

        return {
            "delete_point": delete_point,
            "fill_edges": fill_edges
        }

    def leftSearchForExtractCorners(self, temp_map, p_oint, fill_edges, delete_point, external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect):
        if (p_oint.y + 1 < rect["x"] + rect["height"] and self.map[p_oint.x * self.size_x + p_oint.y + 1] == 0):
            idy = p_oint.x
            idx = p_oint.y + 1
            line = []
            line.append(Point(idy, idx))

            for j in range(idx, rect["x"] + rect["height"]):
                if (temp_map[idy * self.size_x + j] == 0):
                    black_count = 0
                    up_and_down_neighbourhood = [[-1, 0], [1, 0]]

                    for k in range(2):
                        tmp_idy = idy + up_and_down_neighbourhood[k][0]
                        tmp_idx = j + up_and_down_neighbourhood[k][1]

                        if (tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"] or tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(idy, j))
                    else:
                        break
                else:
                    break

            if (is_valid_length and len(line) > 1):
                line.append(p_oint)
                fill_edges.append(line)

                for i in range(len(line)):
                    _up_and_down_neighbourhood = [[-1, 0], [1, 0]]

                    for _k5 in range(2):
                        _tmp_idy5 = line[i].x + \
                            _up_and_down_neighbourhood[_k5][0]

                        _tmp_idx5 = line[i].y + \
                            _up_and_down_neighbourhood[_k5][1]

                        if (_tmp_idy5 < rect["y"] or _tmp_idy5 >= rect["y"] + rect["width"] or _tmp_idx5 < rect["x"] or _tmp_idx5 >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[_tmp_idy5 * self.size_x + _tmp_idx5] == -128 or temp_map[_tmp_idy5 * self.size_x + _tmp_idx5] == inner_corner_value):
                            self.map[_tmp_idy5 * self.size_x + _tmp_idx5] = 127
                            delete_point.append(Point(_tmp_idy5, _tmp_idx5))
            elif (len(line) > valid_length):
                line.append(p_oint)
                fill_edges.append(line)

                for _i8 in range(len(line)):
                    _up_and_down_neighbourhood2 = [[-1, 0], [1, 0]]

                    for _k6 in range(2):
                        _tmp_idy6 = line[_i8].x + \
                            _up_and_down_neighbourhood2[_k6][0]

                        _tmp_idx6 = line[_i8].y + \
                            _up_and_down_neighbourhood2[_k6][1]

                        if (_tmp_idy6 < rect["y"] or _tmp_idy6 >= rect["y"] + rect["width"] or _tmp_idx6 < rect["x"] or _tmp_idx6 >= rect["x"] + rect["height"]):
                            continue

                        if (temp_map[_tmp_idy6 * self.size_x + _tmp_idx6] == -128 or temp_map[_tmp_idy6 * self.size_x + _tmp_idx6] == inner_corner_value):
                            self.map[_tmp_idy6 * self.size_x + _tmp_idx6] = 127
                            delete_point.append(Point(_tmp_idy6, _tmp_idx6))
            else:
                line = []

        return {
            "delete_point": delete_point,
            "fill_edges": fill_edges
        }

    def rightSearchForExtractCorners(self, temp_map, p_oint, fill_edges, delete_point, external_corner_value, inner_corner_value, fill_value, valid_length, is_valid_length, rect):
        if (p_oint.y - 1 > rect["x"] and self.map[p_oint.x * self.size_x + p_oint.y - 1] == 0):
            idy = p_oint.x
            idx = p_oint.y - 1
            line = []
            line.append(Point(idy, idx))

            for j in range(idx, rect["x"], -1):
                if (temp_map[idy * self.size_x + j] == 0):
                    black_count = 0
                    up_and_down_neighbourhood = [[-1, 0], [1, 0]]

                    for k in range(2):
                        tmp_idy = idy + up_and_down_neighbourhood[k][0]
                        tmp_idx = j + up_and_down_neighbourhood[k][1]

                        if (tmp_idx < rect["x"] or tmp_idx >= rect["x"] + rect["height"] or tmp_idy < rect["y"] or tmp_idy >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[tmp_idy * self.size_x + tmp_idx] == -128 or temp_map[tmp_idy * self.size_x + tmp_idx] == external_corner_value or temp_map[tmp_idy * self.size_x + tmp_idx] == inner_corner_value):
                            black_count += 1

                    if (black_count == 1):
                        line.append(Point(idy, j))
                    else:
                        break
                else:
                    break

            if (is_valid_length and len(line) > 1):
                line.append(p_oint)
                fill_edges.append(line)

                for i in range(len(line)):
                    _up_and_down_neighbourhood3 = [[-1, 0], [1, 0]]

                    for _k7 in range(2):
                        _tmp_idy7 = line[i].x + \
                            _up_and_down_neighbourhood3[_k7][0]

                        _tmp_idx7 = line[i].y + \
                            _up_and_down_neighbourhood3[_k7][1]

                        if (_tmp_idx7 < rect["x"] or _tmp_idx7 >= rect["x"] + rect["height"] or _tmp_idy7 < rect["y"] or _tmp_idy7 >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[_tmp_idy7 * self.size_x + _tmp_idx7] == -128 or temp_map[_tmp_idy7 * self.size_x + _tmp_idx7] == inner_corner_value):
                            self.map[_tmp_idy7 * self.size_x + _tmp_idx7] = 127
                            delete_point.append(Point(_tmp_idy7, _tmp_idx7))
            elif (len(line) > valid_length):
                line.append(p_oint)
                fill_edges.append(line)

                for _i9 in range(len(line)):
                    _up_and_down_neighbourhood4 = [[-1, 0], [1, 0]]

                    for _k8 in range(2):
                        _tmp_idy8 = line[_i9].x + \
                            _up_and_down_neighbourhood4[_k8][0]

                        _tmp_idx8 = line[_i9].y + \
                            _up_and_down_neighbourhood4[_k8][1]

                        if (_tmp_idx8 < rect["x"] or _tmp_idx8 >= rect["x"] + rect["height"] or _tmp_idy8 < rect["y"] or _tmp_idy8 >= rect["y"] + rect["width"]):
                            continue

                        if (temp_map[_tmp_idy8 * self.size_x + _tmp_idx8] == -128 or temp_map[_tmp_idy8 * self.size_x + _tmp_idx8] == inner_corner_value):
                            self.map[_tmp_idy8 * self.size_x + _tmp_idx8] = 127
                            delete_point.append(Point(_tmp_idy8, _tmp_idx8))
            else:
                line = []

        return {
            "delete_point": delete_point,
            "fill_edges": fill_edges
        }

    def updateContour(self, contour, delete_points):  # does nothing, really
        return contour

    def fillEdges(self, temp_map, corners, contour, fill_edges, value):
        for i in range(len(fill_edges)):
            edge = fill_edges[i]

            for j in range(len(edge)):
                self.map[edge[j].x * self.size_x + edge[j].y] = value
                contour.append(edge[j])
        return contour

    def removeIndependentRegion(self, all_region, black_boundary, valid_length, rect):
        temp_black_boundary = black_boundary
        return temp_black_boundary

    def fillBlackComponent(self, temp_map, black_region, value):
        for i in range(len(black_region)):
            temp_map[black_region[i].x *
                     self.size_x + black_region[i].y] = value

        return temp_map

    def fillNonBoundaryNoise2(self, nonBoundaryNoise, rect):
        temp_map = (self.map)

        for i in range(len(nonBoundaryNoise)):
            temp_map[nonBoundaryNoise[i].x *
                     self.size_x + nonBoundaryNoise[i].y] = 28

        four_neighbourhood = [[5, 0, 4, 0, 3, 0, 2, 0, 1, 0], [0, 5, 0, 4, 0, 3, 0, 2, 0, 1], [
            -5, 0, -4, 0, -3, 0, -2, 0, -1, 0], [0, -5, 0, -4, 0, -3, 0, -2, 0, -1]]

        for _i12 in range(len(nonBoundaryNoise)):
            for k in range(4):
                tmp_idx5 = nonBoundaryNoise[_i12].y + four_neighbourhood[k][0]
                tmp_idy5 = nonBoundaryNoise[_i12].x + four_neighbourhood[k][1]
                tmp_idx4 = nonBoundaryNoise[_i12].y + four_neighbourhood[k][2]
                tmp_idy4 = nonBoundaryNoise[_i12].x + four_neighbourhood[k][3]
                tmp_idx3 = nonBoundaryNoise[_i12].y + four_neighbourhood[k][4]
                tmp_idy3 = nonBoundaryNoise[_i12].x + four_neighbourhood[k][5]
                tmp_idx2 = nonBoundaryNoise[_i12].y + four_neighbourhood[k][6]
                tmp_idy2 = nonBoundaryNoise[_i12].x + four_neighbourhood[k][7]
                tmp_idx1 = nonBoundaryNoise[_i12].y + four_neighbourhood[k][8]
                tmp_idy1 = nonBoundaryNoise[_i12].x + four_neighbourhood[k][9]

                if (tmp_idy5 < rect["y"] or tmp_idy5 >= rect["y"] + rect["width"] or tmp_idx5 < rect["x"] or tmp_idx5 >= rect["x"] + rect["height"] or tmp_idy4 < rect["y"] or tmp_idy4 >= rect["y"] + rect["width"] or tmp_idx4 < rect["x"] or tmp_idx4 >= rect["x"] + rect["height"] or tmp_idy3 < rect["y"] or tmp_idy3 >= rect["y"] + rect["width"] or tmp_idx3 < rect["x"] or tmp_idx3 >= rect["x"] + rect["height"] or tmp_idy2 < rect["y"] or tmp_idy2 >= rect["y"] + rect["width"] or tmp_idx2 < rect["x"] or tmp_idx2 >= rect["x"] + rect["height"] or tmp_idy1 < rect["y"] or tmp_idy1 >= rect["y"] + rect["width"] or tmp_idx1 < rect["x"] or tmp_idx1 >= rect["x"] + rect["height"]):
                    continue

                if (temp_map[tmp_idy5 * self.size_x + tmp_idx5] == -128 and temp_map[tmp_idy4 * self.size_x + tmp_idx4] == 127 and temp_map[tmp_idy3 * self.size_x + tmp_idx3] == 127 and temp_map[tmp_idy2 * self.size_x + tmp_idx2] == 127 and temp_map[tmp_idy1 * self.size_x + tmp_idx1] == 127):
                    nonBoundaryNoise.append(Point(tmp_idy4, tmp_idx4))
                    nonBoundaryNoise.append(Point(tmp_idy3, tmp_idx3))
                    nonBoundaryNoise.append(Point(tmp_idy2, tmp_idx2))
                    nonBoundaryNoise.append(Point(tmp_idy1, tmp_idx1))
                    break
                elif (temp_map[tmp_idy4 * self.size_x + tmp_idx4] == -128 and temp_map[tmp_idy3 * self.size_x + tmp_idx3] == 127 and temp_map[tmp_idy2 * self.size_x + tmp_idx2] == 127 and temp_map[tmp_idy1 * self.size_x + tmp_idx1] == 127):
                    nonBoundaryNoise.append(Point(tmp_idy3, tmp_idx3))
                    nonBoundaryNoise.append(Point(tmp_idy2, tmp_idx2))
                    nonBoundaryNoise.append(Point(tmp_idy1, tmp_idx1))
                    break
                elif (temp_map[tmp_idy3 * self.size_x + tmp_idx3] == -128 and temp_map[tmp_idy2 * self.size_x + tmp_idx2] == 127 and temp_map[tmp_idy1 * self.size_x + tmp_idx1] == 127):
                    nonBoundaryNoise.append(Point(tmp_idy2, tmp_idx2))
                    nonBoundaryNoise.append(Point(tmp_idy1, tmp_idx1))
                    break
                elif (temp_map[tmp_idy2 * self.size_x + tmp_idx2] == -128 and temp_map[tmp_idy1 * self.size_x + tmp_idx1] == 127):
                    nonBoundaryNoise.append(Point(tmp_idy1, tmp_idx1))
                    break

        for _i13 in range(len(nonBoundaryNoise)):
            self.map[nonBoundaryNoise[_i13].x *
                     self.size_x + nonBoundaryNoise[_i13].y] = -128

        return nonBoundaryNoise

    def roomColorByChain(self, roomChain):
        for row in range(self.size_y):
            for col in range(self.size_x):
                current_map_value = self.map[row * self.size_x + col]

                if (current_map_value == -128):
                    self.map[row * self.size_x + col] = -1
                elif (current_map_value == 127):
                    self.map[row * self.size_x + col] = 1

        for i in range(len(roomChain)):
            self.floodFillSingleChain(
                roomChain[i].points, roomChain[i].roomId)

    def floodFillSingleChain(self, chain_point, value):
        contour_chain_point = []
        dst = []
        row, col = None, None
        init_seed = Point(1, 1)
        contour_chain_point = self.getContourInforChainPoint(
            contour_chain_point, chain_point)

        for i in range(self.size_y):
            for j in range(self.size_x):
                dst.append(value)

        for _i14 in range(len(contour_chain_point)):
            row = contour_chain_point[_i14].y
            col = contour_chain_point[_i14].x
            dst[row * self.size_x + col] = 0

        dst = self.scanLineFloodFill(dst, init_seed, value, 0)

        for _i15 in range(len(contour_chain_point)):
            row = contour_chain_point[_i15].y
            col = contour_chain_point[_i15].x
            dst[row * self.size_x + col] = value

        for _row in range(self.size_y):
            for _col in range(self.size_x):
                current_map_value = self.map[_row * self.size_x + _col]

                if (dst[_row * self.size_x + _col] == value and current_map_value != -1 and current_map_value != 0 and current_map_value != -9):
                    self.map[_row * self.size_x +
                             _col] = dst[_row * self.size_x + _col]

        if (len(contour_chain_point) > 3):
            for _i16 in range(1, len(contour_chain_point) - 1):
                row = contour_chain_point[_i16].y
                col = contour_chain_point[_i16].x

                for di in range(-2, 3):
                    for dj in range(-2, 3):
                        if (row + di < 0 or row + di >= self.size_y or col + dj < 0 or col + dj >= self.size_x):
                            continue
                        else:
                            if (self.map[(row + di) * self.size_x + col + dj] == 1):
                                self.map[(row + di) * self.size_x +
                                         col + dj] = value

    def getContourInforChainPoint(self, contour_chain_point, chain_point):
        for i in range(len(chain_point)):
            point = Point(0, 0)
            point.x = chain_point[i].x
            point.y = chain_point[i].y
            contour_chain_point.append(point)

        return contour_chain_point

    def scanLineFloodFill(self, dst, initial_seed, raw_value, new_value):
        scan_line_seed = []
        scan_line_seed.append(initial_seed)
        tempDst = None

        while (len(scan_line_seed) > 0):
            seed = scan_line_seed[0]
            scan_line_seed.pop(0)
            result1 = self.floodFillLine(dst, seed, -1, raw_value, new_value)
            x_left = result1["boundary"]
            result2 = self.floodFillLine(
                result1["dst"], seed, 1, raw_value, new_value)
            x_right = result2["boundary"]
            tempDst = result2["dst"]
            scan_line_seed = self.searchLineForNewSeed(
                result2["dst"], x_left, x_right, seed.y - 1, raw_value, scan_line_seed)
            scan_line_seed = self.searchLineForNewSeed(
                result2["dst"], x_left, x_right, seed.y + 1, raw_value, scan_line_seed)

        return tempDst

    def floodFillLine(self, dst, initial_seed, direction, raw_value, new_value):
        row = initial_seed.y
        col = initial_seed.x
        boundary = col

        if (direction > 0):
            col += direction

        while (col >= 0 and col < self.size_x):
            if (dst[row * self.size_x + col] == raw_value):
                boundary = col
                dst[row * self.size_x + col] = new_value
                col += direction
            else:
                break

        return {
            "dst": dst,
            "boundary": boundary
        }

    def searchLineForNewSeed(self, dst, x_left, x_right, line_row, raw_value, scan_line_seed):
        if (line_row < 0 or line_row > self.size_y - 1):
            return scan_line_seed

        x_right_copy = x_right
        is_find_seed = False

        while (x_right_copy >= x_left):
            if (dst[line_row * self.size_x + x_right_copy] == raw_value):
                if (not is_find_seed):
                    seed = Point(x_right_copy, line_row)
                    scan_line_seed.append(seed)
                    is_find_seed = True
            else:
                is_find_seed = False

            x_right_copy -= 1

        return scan_line_seed

    def fillInternalObstacles(self):
        if (self.tRect["width"] == 0 and self.tRect["height"] == 0):
            self.tRect = self.findRoiMap(self.tRect)

        contour = []
        internal_obstacles = []
        contour_map = self.map
        result1 = self.extractExternalContoursNewStrategy(
            contour_map, contour, self.tRect)
        contour_map = result1["temp_map"]
        contour = result1["contour"]
        result2 = self.findContourConnectComponent(
            contour_map, contour, self.tRect)
        contour_map = result2["temp_map"]
        contour = result2["contour"]
        contour_map = self.fillBlackComponent(contour_map, contour, 30)
        internal_obstacles = self.findInternalObstacles(
            contour_map, internal_obstacles, self.tRect)
        self.map = self.fillBlackComponent(self.map, internal_obstacles, -9)

    def findContourConnectComponent(self, temp_map, contour, rect):
        eight_neighbourhood = [[-1, 0], [1, 0], [0, -1],
                               [0, 1], [-1, 1], [1, 1], [1, -1], [-1, -1]]
        temp_contour = contour

        while (len(temp_contour) != 0):
            seed = temp_contour[0]
            temp_contour.pop(0)

            for k in range(8):
                temp_idy = seed.x + eight_neighbourhood[k][0]
                temp_idx = seed.y + eight_neighbourhood[k][1]
                if (temp_idy < rect["y"] or temp_idy >= rect["y"] + rect["width"] or temp_idx < rect["x"] or temp_idx >= rect["x"] + rect["height"]):
                    continue

                if (temp_map[temp_idy * self.size_x + temp_idx] == -128):
                    temp_map[temp_idy * self.size_x + temp_idx] = 30
                    temp_contour.append(Point(temp_idy, temp_idx))
                    contour.append(Point(temp_idy, temp_idx))

        return {
            "temp_map": temp_map,
            "contour": contour
        }

    def findInternalObstacles(self, temp_map, point_deque, rect):
        for idy in range(rect["y"], rect["y"] + rect["width"]):
            for idx in range(rect["x"], rect["x"] + rect["height"]):
                if (temp_map[idy * self.size_x + idx] == -128):
                    point_deque.append(Point(idy, idx))
        return point_deque
