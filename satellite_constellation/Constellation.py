"""
Class for holding Constellations of Satellites within it.
"""
from .Satellite import Satellite
from .utils import heavenly_body_radius, heavenly_body_mass, proper_round
import warnings
import math


class Constellation(object):
    """
    Class for describing and holding a constellation of satellites
    """
    def __init__(self, num_sats, num_planes, phasing, inclination, altitude,
                 eccentricity, beam_width, name="Sat", focus="earth", starting_number=0):
        self.num_sats = num_sats
        self.num_planes = num_planes
        self.phasing = phasing
        self.inclination = inclination
        self.altitude = altitude
        self.e = eccentricity
        self.beam = beam_width
        self.start_num = starting_number
        self.constellation_name = name
        self.focus = focus
        self.sats_per_plane, self.correct_phasing = self.__corrected_planes()
        self.perigee_positions = self.__perigee_positions()
        self.raan = self.__calculate_raan()
        self.ta = self.__calculate_ta()
        self.satellites = self.__build_satellites()
        print(self.__str__())

    def __corrected_planes(self):
        sats_per_plane = int(self.num_sats / self.num_planes)
        corrected_phasing = 360 * self.phasing / self.num_sats
        return sats_per_plane, corrected_phasing

    def __perigee_positions(self):
        perigees = list(range(0, 360, int(360/self.sats_per_plane)))
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
            # print('raan', self.raan[i] , ", perigee", self.perigee_positions[i] , ", ta" , self.ta[i])
            satellites.append(Satellite(sat_name, self.altitude, self.e, self.inclination, self.raan[i],
                                        self.perigee_positions[i], self.ta[i], self.beam, focus=self.focus,rads = False))
        return satellites

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}, name={7}, starting_number={8}".format(self.num_sats, self.num_planes,
                                                                                         self.phasing, self.inclination,
                                                                                         self.altitude, self.e,
                                                                                         self.beam,
                                                                                         self.constellation_name,
                                                                                         self.start_num)

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

class SOCConstellation:
    def __init__(self, num_streets, street_width ,altitude, beam_width, raan, eccentricity, revisit_time, name="Sat", focus="earth", starting_number=0): #Start off with just a single polar orbit
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
        self.num_sats, self.sats_per_street = self.__calculate_required_satellites()
        self.perigee = self.__perigee_positions()
        self.ta = self.__calculate_ta()
        print(self.num_sats, "satellites w/", self.angular_spacing , "degrees spacing")
        self.satellites = self.__build_satellites()
        print(self.__str__())

        # self.sats_per_plane, self.correct_phasing = self.__corrected_planes()
        # self.perigee_positions = self.__perigee_positions()
        # self.satellites = self.__build_satellites()

    def __calculate_earth_coverage(self):
        x = self.altitude * math.tan((math.pi/180)*self.beam/2)
        theta = math.asin(x/heavenly_body_radius[self.focus])
        r = heavenly_body_radius[self.focus]*theta
        return r, theta

    def __calculate_spacing(self):
        y = math.sqrt(math.pow(self.earth_coverage_radius,2)-math.pow(self.street_width/2,2))
        ang_spacing = 2*y/heavenly_body_radius[self.focus]
        return y, ang_spacing

    def __calculate_orbit_params(self):
        perigee = heavenly_body_radius[self.focus] + self.altitude  # [km]
        semi_major = perigee/(1-self.eccentricity)  # [km]
        orbital_period = 2*math.pi*math.sqrt((semi_major*10**3)**3/(heavenly_body_mass[self.focus]*6.67408*10**(-11)))  # [s]
        return perigee, semi_major, orbital_period

    def __calculate_required_satellites(self):
        num_satellites_a = self.orbital_period / self.revisit_time  # Calculated from revisit time
        num_satellites_b = math.pi * 2 / self.angular_spacing  # Total coverage

        if num_satellites_a > num_satellites_b:
            num_satellites = num_satellites_b
        else:
            num_satellites = num_satellites_a

        upper = math.ceil(num_satellites)
        lower = math.floor(num_satellites)

        if num_satellites - lower <= 0.01:
            self.angular_spacing = 2*math.pi/lower
            # print(self.angular_spacing)
            self.angular_spacing = 360/lower
            total_satellites = lower*self.num_streets
            return total_satellites, lower
        else:
            # self.angular_spacing = 2*math.pi/upper
            # print(self.angular_spacing)
            self.angular_spacing = 360/upper
            total_satellites = upper*self.num_streets
            return total_satellites, upper

    # def __calculate_ta(self):
    #     ta = [0] * self.num_sats
    #     for idx, anom in enumerate(ta):
    #         ta[idx] = idx*self.angular_spacing
    #     return ta

    def __calculate_ta(self):
        phase = 360/self.num_sats
        print(phase)
        ta = [0] * self.num_sats
        for i in range(self.num_streets):
            for j in range(self.sats_per_street):
                ta[i*self.sats_per_street + j] = i*phase
        return ta

    def __perigee_positions(self):
        perigees = list(range(0, 360, int(360/self.sats_per_street)))
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
                                            self.perigee[j], self.ta[i*self.sats_per_street+j], self.beam, focus=self.focus,rads = False))
        return satellites

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}".format(self.num_sats,
                                                self.altitude, self.eccentricity,
                                                self.beam,
                                                self.constellation_name)

    def __str__(self):
        sat_string = ""
        for sat in self.satellites:
            sat_string += sat.__str__() + '\n'

        return sat_string.rstrip()