"""
Class for holding Constellations of Satellites within it.
"""
from .Satellite import Satellite
import warnings


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
            satellites.append(Satellite(sat_name, self.altitude, self.e, self.inclination, self.raan[i],
                                        self.perigee_positions[i], self.ta[i], self.beam, focus=self.focus))
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

