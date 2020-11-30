"""
Class for holding Constellations of Satellites within it.
"""
from .Satellite import Satellite
from .utils import *
from scipy import optimize
import warnings
import math
import numpy as np


class Constellation:
    """
    Class to implement a layer of abstraction for satellite constellations
    Contains parameters general to all types of satellites

    :param num_satellites: Number of satellites used in the constellation
    :param orbital_period: Orbital period of the satellites
    :param altitude: altitude of the satellites
    :param beam_width: Sensor beam width of the satellites
    :param eccentricity: Eccentricity of satellite orbits
    :param focus: Heavenly body at the focus of the orbit, defaults to 'Earth'
    """

    def __init__(self, num_satellites, orbital_period, altitude, beam_width, eccentricity, focus='Earth'):
        self.num_satellites = num_satellites
        self.orbital_period = orbital_period
        self.altitude = altitude
        self.beam_width = beam_width
        self.eccentricity = eccentricity
        self.focus = focus
        self.satellites = []

    def __str__(self):
        sat_string = ""
        for sat in self.satellites:
            sat_string += sat.__str__() + '\n'

        return sat_string.rstrip()

    def as_dict(self):
        constellation = {}
        for sat in self.satellites:
            if sat.name not in constellation:
                constellation[sat.name] = sat.as_dict()
        constellation['Type'] = 'constellation'
        return constellation

    def as_xml(self):
        warnings.warn("XML support is depreciated and not supported from PIGI 0.8.5 onward", DeprecationWarning)
        return self.as_pigi_output()

    def as_pigi_output(self):
        short_scene = ""
        for sat in self.satellites:
            short_scene += sat.as_xml()

        return short_scene


class WalkerConstellation(Constellation):
    """
    Class for describing and holding a walker constellation of satellites

    :param num_sats: Number of satellites used in the constellation "t"
    :param num_planes: The number of different orbit planes in the constellation "p"
    :param phasing: Dictates the spacing between equivalent satellites in neighbouring orbital planes "f"
    :param inclination: Inclination of orbit relative to equatorial plane "i"
    :param altitude: Altitude of satellites in orbit
    :param eccentricity: Eccentricity of satellite orbits
    :param beam_width: Sensor beam width of the satellites
    :param name: Name of constellation, defaults to "Sat"
    :param focus: Heavenly body at the focus of the orbit, defaults to 'Earth'
    """

    def __init__(self, num_sats, num_planes, phasing, inclination, altitude,
                 eccentricity, beam_width, name="Sat", focus="earth", starting_number=0):
        self.num_sats = num_sats
        self.num_planes = num_planes
        self.phasing = phasing
        self.inclination = inclination
        self.altitude = altitude
        self.eccentricity = eccentricity
        self.beam_width = beam_width
        self.start_num = starting_number
        self.constellation_name = name
        self.focus = focus
        self.sats_per_plane, self.correct_phasing = self.__corrected_planes()
        self.perigee_positions = self.__perigee_positions()
        self.raan = self.__calculate_raan()
        self.ta = self.__calculate_ta()
        self.satellites = self.__build_satellites()

    def __corrected_planes(self):
        sats_per_plane = int(self.num_sats / self.num_planes)
        corrected_phasing = 360 * self.phasing / self.num_sats
        return sats_per_plane, corrected_phasing

    def __perigee_positions(self):
        perigees = list(range(0, 360, int(360 / self.sats_per_plane)))
        all_perigees = []
        for i in range(self.num_sats):
            all_perigees.extend(perigees)
        return all_perigees

    def __calculate_raan(self):
        raan = [0] * self.num_sats
        for i in range(self.sats_per_plane, self.num_sats):
            raan[i] = raan[i - self.sats_per_plane] + 360 / self.num_planes
        return raan

    def __calculate_ta(self):
        ta = [0] * self.num_sats
        for i in range(self.sats_per_plane, self.num_sats):
            ta[i] = ta[i - self.sats_per_plane] + self.correct_phasing
        return ta

    def __build_satellites(self):
        satellites = []
        for i in range(self.num_sats):
            sat_num = i + self.start_num + 1
            sat_name = self.constellation_name + " " + str(sat_num)
            satellites.append(Satellite(sat_name, self.altitude, self.eccentricity, self.inclination, self.raan[i],
                                        self.perigee_positions[i], self.ta[i], self.beam_width, focus=self.focus,
                                        rads=False))
        return satellites

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}, name={7}, starting_number={8}".format(self.num_sats, self.num_planes,
                                                                                         self.phasing, self.inclination,
                                                                                         self.altitude,
                                                                                         self.eccentricity,
                                                                                         self.beam_width,
                                                                                         self.constellation_name,
                                                                                         self.start_num)

    def representation(self):
        return "{0}:{1}/{2}/{3} Walker".format(self.inclination, self.num_sats, self.num_planes, self.phasing)


class SOCConstellation(Constellation):  # Needs to be cleaned up

    """
    Class for describing and holding a walker constellation of satellites

    :param num_streets: Number of satellites used in the constellation
    :param street_width: Angular width of street
    :param altitude: Altitude of satellites in orbit
    :param beam_width: Sensor beam width of the satellites
    :param raan: List of right ascension for each street
    :param eccentricity: Eccentricity of satellite orbits
    :param revisit_time: Time until latitude re-enters satellite coverage
    :param name: Name of constellation, defaults to "Sat"
    :param focus: Heavenly body at the focus of the orbit, defaults to 'Earth'
    """

    def __init__(self, num_streets, street_width, altitude, beam_width, raan, eccentricity, revisit_time, name="Sat",
                 focus="earth", starting_number=0):  # Start off with just a single polar orbit
        self.inclination = 90  # Polar Orbit
        self.num_streets = num_streets
        self.altitude = altitude
        self.beam = beam_width
        self.start_num = starting_number
        self.raan = raan
        self.eccentricity = eccentricity
        self.constellation_name = name
        self.focus = focus
        self.street_width = street_width
        self.revisit_time = revisit_time
        self.earth_coverage_radius, self.earth_coverage_angle = self.__calculate_earth_coverage()
        self.linear_spacing, self.angular_spacing = self.__calculate_spacing()
        self.perigee, self.semi_major, self.orbital_period = self.__calculate_orbit_params()
        self.num_sats, self.sats_per_street, self.true_spacing = self.__calculate_required_satellites()
        self.perigee_positions = self.__perigee_positions()
        self.ta = self.__calculate_ta()
        self.satellites = self.__build_satellites()
        self.longitudinal_drift = self.__calculate_longitudinal_drift()

    def __calculate_earth_coverage(self):
        x = self.altitude * math.tan((math.pi / 180) * self.beam / 2)
        theta = math.asin(x / heavenly_body_radius[self.focus])
        r = heavenly_body_radius[self.focus] * theta
        return r, theta

    def __calculate_spacing(self):
        street_width = heavenly_body_radius[self.focus] * self.street_width * math.pi / 180
        y = math.sqrt(math.pow(self.earth_coverage_radius, 2) - math.pow(street_width / 2, 2))
        ang_spacing = 2 * y / heavenly_body_radius[self.focus]
        return y, ang_spacing

    def __calculate_orbit_params(self):
        perigee = heavenly_body_radius[self.focus] + self.altitude  # [km]
        semi_major = perigee / (1 - self.eccentricity)  # [km]
        orbital_period = 2 * math.pi * math.sqrt(
            (semi_major * 10 ** 3) ** 3 / (heavenly_body_mass[self.focus] * constants['G']))  # [s]
        return perigee, semi_major, orbital_period

    def __calculate_required_satellites(self):
        num_satellites_a = self.orbital_period / self.revisit_time  # Calculated from revisit time
        num_satellites_b = 2 * math.pi / self.angular_spacing  # Total coverage

        if num_satellites_a > num_satellites_b:
            num_satellites = num_satellites_b
        else:
            num_satellites = num_satellites_a

        upper = math.ceil(num_satellites)
        lower = math.floor(num_satellites)

        if num_satellites - lower <= 0.01:
            true_spacing = 360 / lower
            total_satellites = lower * self.num_streets
            return total_satellites, lower, true_spacing
        else:
            true_spacing = 360 / upper
            total_satellites = upper * self.num_streets
            return total_satellites, upper, true_spacing

    def __calculate_ta(self):
        phase = 360 / self.num_sats
        ta = [0] * self.num_sats
        for i in range(self.num_streets):
            for j in range(self.sats_per_street):
                ta[i * self.sats_per_street + j] = i * phase
        return ta

    def __perigee_positions(self):
        perigees = list(range(0, 360, int(360 / self.sats_per_street)))
        all_perigees = []
        for i in range(self.num_sats):
            all_perigees.extend(perigees)
        return all_perigees

    def __build_satellites(self):
        satellites = []
        for i in range(self.num_streets):
            for j in range(self.sats_per_street):
                sat_num = i * self.sats_per_street + j + self.start_num + 1
                sat_name = self.constellation_name + " " + str(sat_num)
                satellites.append(Satellite(sat_name, self.altitude, self.eccentricity, self.inclination, self.raan[i],
                                            self.perigee_positions[j], self.ta[i * self.sats_per_street + j], self.beam,
                                            focus=self.focus, rads=False))
        return satellites

    def __calculate_longitudinal_drift(self):
        drift = self.orbital_period * constants["wE"] * 180 / math.pi
        return drift

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}".format(self.num_sats,
                                                self.altitude, self.eccentricity,
                                                self.beam,
                                                self.constellation_name)


class FlowerConstellation():
    def __init__(self, num_petals, num_days, num_satellites, phasing_n, phasing_d, perigee_argument,
                 inclination, perigee_altitude, focus='earth'):
        self.num_petals = num_petals
        self.num_days = num_days
        self.num_satellites = num_satellites
        self.phasing_n = phasing_n
        self.phasing_d = phasing_d
        self.perigee_argument = perigee_argument
        self.inclination = inclination
        self.perigee_altitude = perigee_altitude
        self.focus = focus
        self.max_sats_per_orbit, self.max_sats = self.__calculate_max_satellites()
        self.raan_spacing, self.mean_anomaly_spacing = self.__calculate_spacing()
        self.num_orbits = self.__calculate_num_orbits()
        self.raan, self.mean_anomaly = self.__calculate_orbits()
        self.revisit_time = self.__calculate_revisit_time()

        # print(self.equations(self.solve()))
        # T, e, a = self.__calculate_orbit_params()
        # temp = optimize.least_squares(self.equations, (e,a), bounds=(( 0, 0), ( 1, np.inf)))
        # print(temp.x)
        # print(temp.cost)
        # print(self.__calculate_orbit_params())

    def __calculate_max_satellites(self):
        return self.num_days, self.phasing_d * self.num_days

    def __calculate_num_orbits(self):
        return math.floor(self.num_satellites / self.num_days)

    def __calculate_spacing(self):
        raan_spacing = -360 * self.phasing_n / self.phasing_d
        mean_anomaly_spacing = -1 * raan_spacing * self.num_petals / self.num_days
        if abs(mean_anomaly_spacing) > 360:
            mean_anomaly_spacing = mean_anomaly_spacing % 360
        if abs(raan_spacing) > 360:
            raan_spacing %= 360

        return raan_spacing, mean_anomaly_spacing

    def __calculate_orbits(self):
        raan = [0]
        M = [0]

        for idx in range(1, min(self.max_sats, self.num_satellites)):
            raan_i = (raan[idx - 1] + self.raan_spacing)
            # if (abs(raan_i) > 360):
            # raan_i = raan_i%360
            raan.append(raan_i)
            M_i = M[idx - 1] + self.mean_anomaly_spacing
            if (abs(M_i) > 360):
                M_i = M_i % 360
            M.append(M_i)
        return raan, M

    def __calculate_orbit_params(self):
        M = heavenly_body_mass[self.focus]
        G = constants["G"]
        R = heavenly_body_radius[self.focus]
        T = self.num_days / self.num_petals * 24 * 60 * 60  # [s]
        a = math.pow(G * M * math.pow(T / (2 * math.pi), 2), 1 / 3)
        e = 1 - (R + self.perigee_altitude) / a
        return T, e, a

    def equations(self, args):
        e, a = args
        R = heavenly_body_radius[self.focus]
        M = heavenly_body_mass[self.focus]
        G = constants["G"]
        J2 = constants["J2E"]
        wE = constants["wE"]
        # print(1)
        T = 28800
        p = + 2 * (R + self.perigee_altitude) - math.pow(+self.perigee_altitude, 2) / a
        epsilon = 3 * math.pow(R, 2) * J2 / (4 * math.pow(p, 2))
        a1 = math.pow(1 + 2 * epsilon * (2 * math.pi / T) * math.cos(self.inclination * math.pi / 180) / wE, -1)
        b = epsilon * (4 + 2 * math.sqrt(1 - math.pow(e, 2)))
        c = (5 + 3 * math.sqrt(1 - math.pow(e, 2))) * math.pow(math.sin(self.inclination * math.pi / 180), 2)
        x = e + (R + self.perigee_altitude) / a1 - 1
        z = -T + 2 * math.pi * math.sqrt(math.pow(a, 3) / (M * G))
        v = (2 * math.pi * self.num_days / (wE * self.num_petals)) * a1 + b - c

        return x, z, v

    def solve(self):
        # T, e, a, p, epsilon = optimize.fsolve(self.equations, (1000, 0, 1000, 1000, 1000))
        T, e, a = optimize.least_squares(self.equations, (1000, 0, 1000, 1000, 1000),
                                         bounds=((0, 0, 0, 0, 0), (1000000, 1, 1000000, 1000000, 1000000)))
        return T, e, a

    def __calculate_revisit_time(self):
        return (self.num_days/self.num_satellites)*24*60*60

    def representation(self, metric = "flower"):
        if metric == "flower":
            i = self.inclination
            return "{0}-{1}-{2}-{3}-{4} Flower".format(self.num_petals, self.num_days, self.num_satellites, self.phasing_n,
                                                   self.phasing_d)
        elif metric == "walker":
            i = self.inclination
            No = self.phasing_d
            G = self.phasing_d*self.num_days/self.max_sats
            Nso = self.max_sats_per_orbit
            t = No * Nso
            p = No
            Nc = (self.num_petals*self.phasing_n+self.phasing_d)/G
            f = int(Nc%No)

            return "{0}:{1}/{2}/{3} Walker".format(i, t, p, f)

