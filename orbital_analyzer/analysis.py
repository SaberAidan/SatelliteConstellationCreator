from satellite_constellation.utils import *
import time


class orbital_analysis:

    def __init__(self, constellation):

        self.constellation = constellation

    def calculate_revisit(self, target_long, target_lat):

        d_prop = 1
        tracking = True
        start_time = -1
        finish_time = -1
        gap_times = []
        view_times = []
        overflow_ctr = 0
        overflow = 10000

        num_sats = []

        while len(gap_times) < 10:

            overflow_ctr += 1

            if overflow_ctr > overflow:
                if len(gap_times) > 0:
                    return np.mean(gap_times), np.std(gap_times)
                else:
                    print("Constant coverage")
                    return 0, 0

            d_theta = d_prop * 0.5 * math.pi / 180

            seconds = (d_theta / (2 * math.pi)) * self.constellation.orbital_period

            longitudinal_drift = seconds * constants["wE"] * 180 / math.pi

            satellites = self.constellation.propagate(d_theta, radians=True)

            in_view = False
            num_sats_round = 0
            start = time.process_time()

            for satellite in satellites:

                if self.check_sat_FOV(satellite, target_lat, target_long, longitudinal_drift):
                    in_view = True
                    num_sats_round += 1

            num_sats.append(num_sats_round)

            if in_view and tracking == False:  # Has just entered FOV
                tracking = True
                start_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                gap_times.append(start_time - finish_time)

            if not in_view and tracking == True:  # Has just left FOV
                tracking = False
                finish_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                view_times.append(finish_time - start_time)

            d_prop += 1
        # return gap_times
        return np.mean(gap_times), np.std(gap_times), np.mean(view_times), np.std(view_times)

    def check_sat_FOV(self, satellite, target_lat, target_lon, drift, threshold=0):  # Convert to numpy

        coords = sat_to_xyz(satellite)
        spherical_coords = cart2polar(coords[0], coords[1], coords[2])
        geographic_coords = spherical2geographic(spherical_coords[1], spherical_coords[2], radians=True)

        lat_i = geographic_coords[0]
        long_i = geographic_coords[1]
        long_i = long_i - drift

        if long_i > 180:
            long_i = long_i - 180
        elif long_i < -180:
            long_i = long_i + 360

        # print('normal', lat_i, long_i)

        distance = geographic_distance(target_lat, target_lon, lat_i, long_i,
                                       heavenly_body_radius[self.constellation.focus],
                                       radians=False)

        threshold, theta = self.calculate_satellite_coverage(satellite)

        # if threshold == 0:
        #     if self.constellation.earth_coverage_radius == 0:
        #         threshold, theta = self.calculate_satellite_coverage(satellite)
        #     else:
        #         threshold = self.constellation.earth_coverage_radius

        if distance < threshold:
            in_view = True
        else:
            in_view = False

        return in_view

    def calculate_satellite_coverage(self, satellite):  # Convert to numpy

        if satellite.eccentricity > 0:

            a = satellite.semi_major
            b = a * math.sqrt(1 - math.pow(self.constellation.eccentricity, 2))
            f = (self.constellation.altitude + heavenly_body_radius[self.constellation.focus]) * 10 ** 3
            disp = a - f

            ang = deg_2_rad(satellite.ta + 180)
            x_i, y_i, z_i = disp + a * np.cos(ang), b * np.sin(ang), 0
            coords = np.array([x_i, y_i, z_i]) * 10 ** -3
            coords = rotate(coords, deg_2_rad(satellite.right_ascension), 'z')
            coords = rotate(coords, deg_2_rad(satellite.inclination), 'x')

            r = math.sqrt(np.sum(coords ** 2))

        else:

            r = satellite.true_alt

        half_width = deg_2_rad(satellite.beam / 2)

        # print('normal hw',half_width)

        max_width = math.atan(heavenly_body_radius[self.constellation.focus] / (
                self.constellation.altitude + heavenly_body_radius[self.constellation.focus]))

        if half_width > max_width:
            half_width = max_width

        x = r * math.tan(half_width)
        if x > heavenly_body_radius[self.constellation.focus]:
            x = heavenly_body_radius[self.constellation.focus]

        theta = math.asin(x / heavenly_body_radius[self.constellation.focus])

        coverage_radius = heavenly_body_radius[self.constellation.focus] * theta

        return coverage_radius, theta

    def find_links(self, custom_satellites=None):  # Convert to numpy

        sat_coords = self.constellation.as_cartesian(self.constellation.satellites)
        average_links = np.array([])

        for idy in range(0, sat_coords.shape[0]):  # Draw line of sight between satellites
            ctr = 0
            for idz in range(0, sat_coords.shape[0]):
                if idz != idy:
                    temp_coords = np.append([sat_coords[idy, :]], [sat_coords[idz, :]], axis=0)
                    if not sphere_intercept(temp_coords[0], temp_coords[1],
                                            heavenly_body_radius[self.constellation.focus]):
                        ctr += 1
            average_links = np.append(average_links, ctr)

        return math.floor(np.mean(average_links))

    def calculate_constellation_coverage(self, resolution=0.5, custom_satellites=None):

        d_lat = -90
        area = 0

        coords = self.constellation.as_geographic()
        sats = self.constellation.satellites

        while d_lat < 90:
            d_long = -180
            while d_long < 180:
                centre_lat = (2 * d_lat + resolution) / 2
                centre_long = (2 * d_long + resolution) / 2
                d_sat = 0

                for coord in coords:

                    distance = geographic_distance(centre_lat, centre_long, coord[0], coord[1],
                                                   heavenly_body_radius[self.constellation.focus], radians=False)

                    if self.constellation.earth_coverage_radius == 0:
                        earth_coverage_radius, theta = self.constellation.calculate_satellite_coverage(sats[d_sat])
                    else:
                        earth_coverage_radius = self.constellation.earth_coverage_radius

                    if distance < earth_coverage_radius:
                        area += geographic_area(d_lat, d_long, d_lat + resolution, d_long + resolution,
                                                heavenly_body_radius[self.constellation.focus],
                                                radians=False)

                        break
                    d_sat += 1

                d_long += resolution

            d_lat += resolution

        return area  # Convert to numpy

    def check_ground_FOV(self, satellite, lat, lon, FOV, drift, radians=False):

        if not radians:
            FOV = deg_2_rad(FOV)

        lon = lon + drift

        if lon > 180:
            lon = lon - 180
        elif lon < -180:
            lon = lon + 360

        lon = (lon + 360) % 360

        if lon > 180:
            lon -= 360

        spherical = geographic2spherical(lat, lon, heavenly_body_radius[self.constellation.focus])

        coords = sat_to_xyz(satellite)
        alt = satellite.true_alt
        max_radius = alt * FOV / 2

        cartesian = polar2cart(alt, spherical[0], spherical[1])

        distance = math.sqrt(
            math.pow(coords[0] - cartesian[0], 2) + math.pow(coords[1] - cartesian[1], 2) + math.pow(
                coords[2] - cartesian[2], 2))

        if distance > max_radius:
            return False
        else:
            return True  # Convert to numpy

    def ground_station_visibility(self, station_lat, station_long, station_fov, radians=False):

        d_prop = 1
        tracking = True
        start_time = -1
        finish_time = -1
        gap_times = []
        view_times = []
        overflow_ctr = 0
        overflow = 10000

        num_sats = []

        while len(gap_times) < 10:

            overflow_ctr += 1

            if overflow_ctr > overflow:
                if len(gap_times) > 0:
                    return np.mean(gap_times), np.std(gap_times)
                else:
                    print("Constant coverage")
                    return 0, 0

            d_theta = d_prop * 0.5 * math.pi / 180

            seconds = (d_theta / (2 * math.pi)) * self.constellation.orbital_period
            longitudinal_drift = seconds * constants["wE"] * 180 / math.pi

            satellites = self.constellation.propagate(d_theta, radians=True)

            in_view = False
            num_sats_round = 0

            for satellite in satellites:

                if self.check_ground_FOV(satellite, station_lat, station_long, station_fov, longitudinal_drift,
                                         radians):
                    in_view = True
                    num_sats_round += 1

            num_sats.append(num_sats_round)

            if in_view and tracking == False:  # Has just entered FOV
                tracking = True
                start_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                gap_times.append(start_time - finish_time)

            if not in_view and tracking == True:  # Has just left FOV
                tracking = False
                finish_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                view_times.append(finish_time - start_time)

            d_prop += 1

        return np.mean(num_sats), np.mean(gap_times), np.std(gap_times), np.mean(view_times), np.std(
            view_times)  # Convert to numpy

    def revisit_numpy(self, target_long, target_lat, step_res=0.5, num_iterations=10000, num_samples=10):

        """

        Function to numerically determine the revisit time of the satellite constellation, i.e, how long until a
        position on the orbited bodies surface is in view of a satellite again. \n
        Works by propagating the orbit forward by incrementing the true anomaly of each satellite and checking if the
        point is within the view of the satellites beam.

        :param target_long: Longitude of the target
        :param target_lat: Latitude of the target
        :param step_res: Angle the satellites true anomaly is incremented by [radians]
        :param num_iterations: Maximum number of times the orbit is propagated forward
        :param num_samples: Number of complete exits from view to be recorded


        """

        d_prop = 1
        tracking = True
        start_time = -1
        finish_time = -1
        gap_times = []
        view_times = []
        overflow_ctr = 0
        overflow = num_iterations

        num_sats = []

        while len(gap_times) < num_samples:

            overflow_ctr += 1

            if overflow_ctr > overflow:
                if len(gap_times) > 0:
                    return np.mean(gap_times), np.std(gap_times)
                else:
                    print("Constant coverage")
                    return 0, 0

            d_theta = d_prop * step_res * math.pi / 180

            seconds = (d_theta / (2 * math.pi)) * self.constellation.orbital_period

            longitudinal_drift = seconds * constants["wE"] * 180 / math.pi

            satellites = np.array(self.constellation.propagate_np(d_theta, radians=True))

            in_view = self.check_sat_FOV_np(satellites, target_lat, target_long, longitudinal_drift)

            if in_view.any():
                num_sats.append(np.sum(in_view))

            if in_view.any() and tracking == False:  # Has just entered FOV
                tracking = True
                start_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                gap_times.append(start_time - finish_time)
            elif not in_view.any() and tracking == True:  # Has just left FOV
                tracking = False
                finish_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                view_times.append(finish_time - start_time)

            d_prop += 1

        res = {
            "gap_durations": gap_times,
            "mean_gap": np.mean(gap_times),
            "std_gap": np.std(gap_times),
            "cover_durations": view_times,
            "mean_cover": np.mean(view_times),
            "std_cover": np.std(view_times),
            "gap_cover_ratio": np.mean(gap_times) / np.mean(view_times),
            "sats_in_view": num_sats,
            "mean_sats_in_view": np.mean(num_sats),
            "std_sats_in_view": np.std(num_sats)
        }

        return res

    def check_sat_FOV_np(self, satellites, target_lat, target_lon, drift, threshold=0):

        """

        Function to determine if a position on the orbited bodies surface is within view of any of the passed satellites.\n
        Works by determining the centre of the satellites beam on the earths surface and calculating the geodesic
        distance to the point of interest and determining if it is within the satellites beam.

        :param satellites: List of satellites :param target_lon: Longitude of the target
        :param target_lat: Latitude of the target
        :param target_lon: Longitude of the target
        :param drift: How far the earth has rotated relative to its initial rotation in the satellite frame
        :param threshold: Allowable distance from beam centre to point of interest. If set to 0, this will be determined by the beam width, altitude of the satellite and radius of the heavenly body.

        """

        coords = self.constellation.as_geographic_np(satellites)

        coords[:, 1] = coords[:, 1] - drift
        coords[coords > 180] = coords[coords > 180] - 180
        coords[coords < -180] = coords[coords < -180] + 360

        coords = coords * math.pi / 180

        dist = geographic_distance_np(target_lat, target_lon, coords[:, 0], coords[:, 1],
                                      heavenly_body_radius[self.constellation.focus], radians=True)

        if threshold == 0:
            threshold = self.calculate_satellite_coverage_np(self.constellation.as_numpy(satellites))[0]

        in_view = dist < threshold

        return in_view

    def calculate_satellite_coverage_np(self, satellites):

        """

        Function to determine the individual coverage of a list of satellites on the surface of the orbited body.\n

        :param satellites: List of satellites to get coverage for
        """

        coords = sat_to_xyz_np(satellites)
        r = np.sqrt(np.sum(np.square(coords), axis=1))

        true_alt = satellites[:, 1]
        alt = satellites[:, 0]

        half_width = (satellites[:, 7] / 2) * math.pi / 180

        max_width = np.arctan((true_alt - alt) / true_alt)

        half_width[half_width > max_width] = max_width[half_width > max_width]

        x = r * np.tan(half_width)

        x[x > (true_alt - alt)] = (true_alt - alt)[x > (true_alt - alt)]

        theta = np.arcsin(x / (true_alt - alt))

        coverage_radius = (true_alt - alt) * theta

        return coverage_radius, theta

    def find_links_np(self, custom_satellites=None):

        """

          Function to determine how many other satellites in the constellation each satellite can see.\n
          Needs to be reworked to propagate the orbit forward and calculate the number of links.

          :param custom_satellites: List of satellites to check

          """

        if custom_satellites is None:
            sats = np.array(self.constellation.satellites)
        else:
            sats = np.array(custom_satellites)

        sat_coords = self.constellation.as_cartesian_np(sats)

        link_nums = np.zeros(len(sat_coords))
        d_sat = 0

        for coords in sat_coords:
            link_nums[d_sat] = np.sum(sphere_intercept_np(coords, sat_coords, 6371) == False)
            d_sat += 1

        return np.mean(link_nums)

    def calculate_constellation_coverage_np(self, custom_satellites=None, resolution=0.5):

        """

          Function to determine how much surface area of the orbited body is covered by the constellation.\n
          Works by splitting the surface of the earth into latitude/longitude blocks and determining if any of the
          satellites beams cover the centre of the block. The area of every covered block is then combined. \n
          Needs to be reworked to propagate the orbit and get the average coverage.

          :param custom_satellites: List of satellites to use in coverage calculation
          :param resolution:

          """

        area = 0
        coords = self.constellation.as_geographic()

        lat = np.arange(-90, 90, resolution)
        lon = np.arange(-180, 180, resolution)

        centres = np.array(np.meshgrid(lat, lon)).T.reshape(180, 360, 2) + 0.5
        print(np.shape(centres))

        threshold = self.calculate_satellite_coverage_np(self.constellation.as_numpy())
        # print(threshold[0])

        for i in centres:
            for j in i:
                distance = geographic_distance_np(j[0], j[1], coords[:, 0], coords[:, 1],
                                                  heavenly_body_radius[self.constellation.focus], radians=False)
                # hits = np.where(distance < threshold)

                if np.any(distance < threshold[0]):
                    area_i = geographic_area(j[0] - 0.5, j[1] - 0.5, j[0] + 0.5, j[1] + 0.5,
                                             heavenly_body_radius[self.constellation.focus], radians=False)
                    area += area_i
        # print(threshold[0])

        return area  # Convert to numpy

    def check_ground_FOV_np(self, satellites, lat, lon, FOV, drift, radians=False):

        """

        Function to determine if any of the passed satellites are in view of a given ground station.\n
        Works by projecting the beam of the ground station onto a celestial sphere matching the altitude of each
        satellite. The geodesic distance from the centre of this projected circle to the position of the satellite is
        then calculated and checked against the radius of the circle to determine if it is in view. \n
        Needs to be reworked to take into account both the beam of ground station and satellite?

        :param satellites: List of satellites :param target_lon: Longitude of the target
        :param lat: Latitude of the ground station
        :param lon: Longitude of the ground station
        :param drift: How far the earth has rotated relative to its initial rotation in the satellite frame
        :param FOV: Field of view of the ground station
        :param radians: Indicates whether the FOV is passed in radians or not

        """

        if not radians:
            FOV = deg_2_rad(FOV)

        lon = lon + drift

        if lon > 180:
            lon = lon - 180
        elif lon < -180:
            lon = lon + 360

        lon = (lon + 360) % 360

        if lon > 180:
            lon -= 360

        spherical = geographic2spherical(lat, lon, heavenly_body_radius[self.constellation.focus])

        coords = sat_to_xyz_np(satellites)
        alt = np.sqrt(np.sum(np.square(coords), axis=1))
        max_radius = alt * FOV / 2

        cartesian = polar2cart(alt, spherical[0], spherical[1])

        distance = np.sqrt(
            np.power(coords[:, 0] - cartesian[0], 2) + np.power(coords[:, 1] - cartesian[1], 2) + np.power(
                coords[:, 2] - cartesian[2], 2))

        return distance < max_radius

    def ground_station_visibility_np(self, station_lat, station_long, station_fov, radians=False, step_res=0.5,
                                     num_iterations=10000, num_samples=10):

        """

        Function to numerically determine the revisit time of the satellite constellation, i.e, how long until a
        position on the orbited bodies surface is in view of a satellite again. \n
        Works by propagating the orbit forward by incrementing the true anomaly of each satellite and checking if the
        point is within the view of the satellites beam.

        :param station_lat: Longitude of the target
        :param station_long: Latitude of the target
        :param station_fov: Station field of view
        :param radians: Indicates whether the FOV is passed in radians or not
        :param step_res: Angle the satellites true anomaly is incremented by [radians]
        :param num_iterations: Maximum number of times the orbit is propagated forward
        :param num_samples: Number of complete exits from view to be recorded


        """

        d_prop = 1
        tracking = True
        start_time = -1
        finish_time = -1
        gap_times = []
        view_times = []
        overflow_ctr = 0
        overflow = num_iterations

        num_sats = []

        while len(gap_times) < num_samples:

            overflow_ctr += 1

            if overflow_ctr > overflow:
                if len(gap_times) > 0:
                    break
                    # return np.mean(gap_times), np.std(gap_times)
                else:
                    print("Constant coverage")
                    return 0, 0

            d_theta = d_prop * step_res * math.pi / 180

            seconds = (d_theta / (2 * math.pi)) * self.constellation.orbital_period
            longitudinal_drift = seconds * constants["wE"] * 180 / math.pi

            satellites = self.constellation.as_numpy(self.constellation.propagate_np(d_theta, radians=True))

            in_view = self.check_ground_FOV_np(satellites, station_lat, station_long, station_fov, longitudinal_drift,
                                               radians)

            if in_view.any():
                num_sats.append(np.sum(in_view))

            if in_view.any() and tracking == False:  # Has just entered FOV
                tracking = True
                start_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                gap_times.append(start_time - finish_time)
            elif not in_view.any() and tracking == True:  # Has just left FOV
                tracking = False
                finish_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                view_times.append(finish_time - start_time)

            d_prop += 1

        res = {
            "gap_durations": gap_times,
            "mean_gap": np.mean(gap_times),
            "std_gap": np.std(gap_times),
            "cover_durations": view_times,
            "mean_cover": np.mean(view_times),
            "std_cover": np.std(view_times),
            "gap_cover_ratio": np.mean(gap_times) / np.mean(view_times),
            "sats_in_view": num_sats,
            "mean_sats_in_view": np.mean(num_sats),
            "std_sats_in_view": np.std(num_sats)
        }
        return res
