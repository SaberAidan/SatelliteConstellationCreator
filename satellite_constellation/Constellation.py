"""
Class for holding Constellations of Satellites within it.
"""
from .Satellite import Satellite
from .utils import *
import warnings
import math


class Constellation:
    """
    Class to implement a layer of abstraction for satellite constellations
    Contains parameters general to all types of satellites

    :param num_satellites: Number of satellites used in the constellation
    :param orbital_period: Orbital period of the satellites [s]
    :param altitude: altitude of the satellites [km]
    :param beam_width: Sensor beam width of the satellites [degrees]
    :param eccentricity: Eccentricity of satellite orbits
    :param inclination: Inclination of orbit relative to equatorial plane "i" [degrees]
    :param focus: Heavenly body at the focus of the orbit, defaults to 'Earth'
    """

    def __init__(self, num_satellites, orbital_period, altitude, beam_width, eccentricity, inclination, focus='Earth',
                 name="constellation"):
        self.num_satellites = num_satellites
        self.orbital_period = orbital_period
        self.altitude = altitude
        self.beam_width = beam_width
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.focus = focus
        self.name = name
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


class WalkerConstellation(Constellation):  # Walker delta pattern
    """
    Class for describing and holding a walker constellation of satellites

    :param num_sats: Number of satellites used in the constellation "t"
    :param num_planes: The number of different orbit planes in the constellation "p"
    :param phasing: Dictates the spacing between equivalent satellites in neighbouring orbital planes "f"
    :param inclination: Inclination of orbit relative to equatorial plane "i" [degrees]
    :param altitude: Altitude of satellites in orbit [km]
    :param eccentricity: Eccentricity of satellite orbits
    :param beam_width: Sensor beam width of the satellites [degrees]
    :param name: Name of constellation, defaults to "Sat"
    :param focus: Heavenly body at the focus of the orbit, defaults to 'Earth'
    """

    def __init__(self, num_sats, num_planes, phasing, inclination, altitude,
                 eccentricity, beam_width, name="Sat", focus="earth", starting_number=0):
        super(WalkerConstellation, self).__init__(num_sats, 0, altitude, beam_width, eccentricity, inclination, focus,
                                                  name)
        self.num_planes = num_planes
        self.phasing = phasing
        self.start_num = starting_number
        self.constellation_name = name
        self.sats_per_plane, self.correct_phasing = self.__corrected_planes()
        self.perigee_positions = self.__perigee_positions()
        self.raan = self.__calculate_raan()
        self.ta = self.__calculate_ta()
        self.satellites = self.__build_satellites()

    def __corrected_planes(self):
        sats_per_plane = int(self.num_satellites / self.num_planes)
        corrected_phasing = 360 * self.phasing / self.num_satellites
        return sats_per_plane, corrected_phasing

    def __perigee_positions(self):
        perigees = list(range(0, 360, int(360 / self.sats_per_plane)))
        all_perigees = []
        for i in range(self.num_planes):
            all_perigees.extend(perigees)
        return all_perigees

    def __calculate_raan(self):
        raan = [0] * self.num_satellites
        for i in range(self.sats_per_plane, self.num_satellites):
            raan[i] = raan[i - self.sats_per_plane] + 360 / self.num_planes
        return raan

    def __calculate_ta(self):
        ta = [0] * self.num_satellites
        for i in range(self.sats_per_plane, self.num_satellites):
            ta[i] = ta[i - self.sats_per_plane] + self.correct_phasing
        return ta

    def __build_satellites(self):
        satellites = []
        for i in range(self.num_satellites):
            sat_num = i + self.start_num + 1
            sat_name = self.constellation_name + " " + str(sat_num)
            satellites.append(Satellite(sat_name, self.altitude, self.eccentricity, self.inclination, self.raan[i],
                                        self.perigee_positions[i], self.ta[i], self.beam_width, focus=self.focus,
                                        rads=False))
        return satellites

    # def __calculate_revisit_time(self):  # Method from https://arxiv.org/pdf/1807.02021.pdf
    # longitudinal_
    # return revisit_time

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}, name={7}, starting_number={8}".format(self.num_satellites,
                                                                                         self.num_planes,
                                                                                         self.phasing, self.inclination,
                                                                                         self.altitude,
                                                                                         self.eccentricity,
                                                                                         self.beam_width,
                                                                                         self.constellation_name,
                                                                                         self.start_num)

    def representation(self):
        return "{0}:{1}/{2}/{3} Walker".format(self.inclination, self.num_satellites, self.num_planes, self.phasing)


class SOCConstellation(Constellation):  # This is really just a part of a walker star pattern

    """
    Class for describing and holding a walker constellation of satellites

    :param num_streets: Number of satellites used in the constellation
    :param street_width: Angular width of street [degrees]
    :param altitude: Altitude of satellites in orbit [km]
    :param beam_width: Sensor beam width of the satellites [degrees]
    :param raan: List of right ascension for each street [degrees]
    :param eccentricity: Eccentricity of satellite orbits
    :param revisit_time: Time until latitude re-enters satellite coverage [s]
    :param name: Name of constellation, defaults to "Sat"
    :param focus: Heavenly body at the focus of the orbit, defaults to 'Earth'
    """

    def __init__(self, num_streets, street_width, altitude, beam_width, raan, eccentricity, revisit_time, name="Sat",
                 focus="earth", starting_number=0):

        super(SOCConstellation, self).__init__(0, 0, altitude, beam_width, eccentricity, 90, focus, name)

        self.num_streets = num_streets
        self.start_num = starting_number
        self.raan = raan
        self.street_width = street_width
        self.revisit_time = revisit_time
        self.earth_coverage_radius, self.earth_coverage_angle = self.__calculate_earth_coverage()
        self.linear_spacing, self.angular_spacing = self.__calculate_spacing()
        self.perigee, self.semi_major, self.orbital_period = self.__calculate_orbit_params()
        self.num_satellites, self.sats_per_street, self.true_spacing = self.__calculate_required_satellites()
        self.perigee_positions = self.__perigee_positions()
        self.ta = self.__calculate_ta()
        self.satellites = self.__build_satellites()
        self.longitudinal_drift = self.__calculate_longitudinal_drift()

    def __calculate_earth_coverage(self):
        x = self.altitude * math.tan((math.pi / 180) * self.beam_width / 2)
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
        phase = 360 / self.num_satellites
        ta = [0] * self.num_satellites
        for i in range(self.num_streets):
            for j in range(self.sats_per_street):
                ta[i * self.sats_per_street + j] = i * phase
        return ta

    def __perigee_positions(self):
        perigees = list(range(0, 360, int(360 / self.sats_per_street)))
        all_perigees = []
        for i in range(self.num_satellites):
            all_perigees.extend(perigees)
        return all_perigees

    def __build_satellites(self):
        satellites = []
        for i in range(self.num_streets):
            for j in range(self.sats_per_street):
                sat_num = i * self.sats_per_street + j + self.start_num + 1
                sat_name = self.name + " " + str(sat_num)
                satellites.append(Satellite(sat_name, self.altitude, self.eccentricity, self.inclination, self.raan[i],
                                            self.perigee_positions[j], self.ta[i * self.sats_per_street + j],
                                            self.beam_width, focus=self.focus, rads=False))
        return satellites

    def __calculate_longitudinal_drift(self):
        drift = self.orbital_period * constants["wE"] * 180 / math.pi
        return drift

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}".format(self.num_satellites,
                                                self.altitude, self.eccentricity,
                                                self.beam_width,
                                                self.name)

    def representation(self):
        return "{0}:{1}/{2}/{3} Walker".format(self.inclination, self.num_satellites, self.num_streets, 0)


class FlowerConstellation(Constellation):
    """
    Class for describing and holding a flower constellation of satellites

    :param num_petals: Number of petals formed when viewing the relative orbits
    :param num_days: The number of days for the constellation to completely repeat its orbit
    :param num_satellites: The desired number of satellites involved in the constellation
    :param phasing_n: Phasing parameter n, effects allowable satellite positions
    :param phasing_d: Phasing parameter d, effects the maximum number of satellites
    :param perigee_argument: Argument of perigee for satellite orbits [degrees]
    :param inclination: Inclination of orbit relative to equatorial plane [degrees]
    :param perigee_altitude: Altitude of perigee for satellite orbits [km]
    :param beam_width: Angular width of satellite sensor beam [degrees]
    :param focus: Object at focus of orbit, defaults to 'earth'
    """

    def __init__(self, num_petals, num_days, num_satellites, phasing_n, phasing_d, perigee_argument,
                 inclination, perigee_altitude, beam_width, focus='earth', name="constellation"):
        super(FlowerConstellation, self).__init__(num_satellites, 0, perigee_altitude, beam_width, 0, inclination,
                                                  focus, name)
        self.num_petals = num_petals
        self.num_days = num_days
        self.phasing_n = phasing_n
        self.phasing_d = phasing_d
        self.perigee_argument = perigee_argument
        self.max_sats_per_orbit, self.max_sats = self.__calculate_max_satellites()
        self.raan_spacing, self.mean_anomaly_spacing = self.__calculate_spacing()
        self.num_orbits = self.__calculate_num_orbits()
        self.orbital_period, self.eccentricity, self.semi_major = self.__calculate_orbit_params()
        self.raan, self.mean_anomaly, self.true_anomaly = self.__calculate_orbits()
        self.revisit_time = self.__calculate_revisit_time()
        self.minimum_revisit_time = self.__calculate_minimum_revisit_time()
        self.satellites = self.__build_satellites()

    def __calculate_max_satellites(self):
        max_sats = self.phasing_d*self.num_days
        if max_sats < self.num_satellites:
            self.num_satellites = max_sats
        return self.num_days, max_sats

    def __calculate_num_orbits(self):
        # return math.floor(self.num_satellites / self.num_days)
        return self.phasing_d

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
        v = [0]
        testList = [[0,0]]
        for idx in range(1, min(self.max_sats, self.num_satellites)):
            raan_i = (raan[idx - 1] + self.raan_spacing)
            if raan_i < 0:
                raan_i = 360 + raan_i
            raan.append(raan_i)
            M_i = M[idx - 1] + self.mean_anomaly_spacing
            if (abs(M_i) > 360):
                M_i = M_i % 360
            if M_i < 0:
                M_i = 360 + M_i
            M.append(M_i)
            v_i = M_i + (2 * self.eccentricity - 0.25 * math.pow(self.eccentricity, 3)) * math.sin(M_i * math.pi / 180)
            v.append(v_i)
        return raan, M, v

    def __calculate_orbit_params(self):
        M = heavenly_body_mass[self.focus]
        G = constants["G"]
        wE = constants["wE"]
        T = (2 * math.pi / wE) * (self.num_days / self.num_petals)
        a = math.pow(G * M * math.pow(T / (2 * math.pi), 2), 1 / 3)
        e = 1 - (heavenly_body_radius[self.focus] * 10 ** 3 + self.altitude) / a

        return T, e, a

    def __calculate_revisit_time(self):
        revisit_time = self.num_days / self.num_satellites
        return revisit_time

    def __calculate_minimum_revisit_time(self):
        min_revisit_time = self.num_days / self.max_sats
        return min_revisit_time

    def __build_satellites(self):
        satellites = []
        for i in range(self.num_satellites):
            sat_name = i
            satellites.append(Satellite(sat_name,self.altitude,self.eccentricity,self.inclination,self.raan[i],
                                        self.true_anomaly[i],self.true_anomaly[i],self.beam_width,self.focus,rads = True))

        return satellites

    def representation(self, representation_type="flower"):
        if representation_type == "flower":
            i = self.inclination
            return "{0}-{1}-{2}-{3}-{4} Flower".format(self.num_petals, self.num_days, self.num_satellites,
                                                       self.phasing_n,
                                                       self.phasing_d)
        elif representation_type == "walker":
            i = self.inclination
            No = self.phasing_d
            G = self.phasing_d * self.num_days / self.max_sats
            Nso = self.max_sats_per_orbit
            t = No * Nso
            p = No
            Nc = (self.num_petals * self.phasing_n + self.phasing_d) / G
            f = int(Nc % No)

            return "{0}:{1}/{2}/{3} Walker".format(i, t, p, f)
