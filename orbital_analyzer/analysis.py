from satellite_constellation.Constellation import *


class orbital_analysis:

    def __init__(self, constellation):

        self.constellation = constellation
        self.num_satellites = constellation.num_satellites
        self.orbital_period = constellation.orbital_period
        self.altitude = constellation.altitude
        self.beam_width = constellation.beam_width
        self.eccentricity = constellation.eccentricity
        self.inclination = constellation.inclination
        self.focus = constellation.focus
        self.name = constellation.name
        self.satellites = []
        self.earth_coverage_radius = 0

    def calculate_revisit(self, target_long, target_lat):

        d_prop = 1
        tracking = True
        start_time = -1
        finish_time = -1
        gap_times = []
        view_times = []
        overflow_ctr = 0
        overflow = 10000

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
            coords = self.constellation.as_geographic(satellites)

            in_view = False

            d_sat = 0

            for coord in coords:

                lat_i = coord[0]
                long_i = coord[1]
                long_i = long_i - longitudinal_drift

                if long_i > 180:
                    long_i = long_i - 180
                elif long_i < -180:
                    long_i = long_i + 360

                distance = geographic_distance(target_lat, target_long, lat_i, long_i, radians=False)

                if self.constellation.earth_coverage_radius == 0:
                    earth_coverage_radius, theta = self.calculate_satellite_coverage(satellites[d_sat])
                else:
                    earth_coverage_radius = self.constellation.earth_coverage_radius

                if distance < earth_coverage_radius:
                    in_view = True

                d_sat += 1

            if in_view and tracking == False:  # Has just entered FOV
                tracking = True
                start_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                gap_times.append(start_time - finish_time)

            if not in_view and tracking == True:  # Has just left FOV
                tracking = False
                finish_time = self.constellation.orbital_period * d_theta / (2 * math.pi)
                view_times.append(finish_time - start_time)

            d_prop += 1

        return np.mean(gap_times), np.std(gap_times)

    def calculate_satellite_coverage(self, satellite):

        if satellite.eccentricity > 0:

            a = satellite.semi_major
            b = a * math.sqrt(1 - math.pow(self.constellation.eccentricity, 2))
            f = (self.constellation.altitude + heavenly_body_radius[self.constellation.focus]) * 10 ** 3
            disp = a - f

            ang = (satellite.ta + 180) * math.pi / 180
            x_i, y_i, z_i = disp + a * np.cos(ang), b * np.sin(ang), 0
            coords = np.array([x_i, y_i, z_i]) * 10 ** -3
            coords = rotate(coords, satellite.right_ascension * math.pi / 180, 'z')
            coords = rotate(coords, satellite.inclination * math.pi / 180, 'x')

            r = math.sqrt(np.sum(coords ** 2))

        else:

            r = satellite.true_alt

        half_width = (satellite.beam / 2) * math.pi / 180
        max_width = math.atan(heavenly_body_radius[self.constellation.focus] / (self.constellation.altitude + heavenly_body_radius[self.constellation.focus]))
        if half_width > max_width:
            half_width = max_width
        x = r * math.tan(half_width)
        if x > heavenly_body_radius[self.constellation.focus]:
            x = heavenly_body_radius[self.constellation.focus]

        theta = math.asin(x / heavenly_body_radius[self.constellation.focus])
        coverage_radius = heavenly_body_radius[self.constellation.focus] * theta
        # print(coverage_radius)
        return coverage_radius, theta

    def find_links(self, custom_satellites=None):

        # sat_coords = np.array([[0, 0, 0]])
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

                    distance = geographic_distance(centre_lat, centre_long, coord[0], coord[1], radians=False)

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

        return area
