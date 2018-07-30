"""
A class for creating a satellite object, describing the characteristics of it.
"""

from math import pi
from .utils import heavenly_body_radius
import warnings


class Satellite(object):

    def __init__(self, name, altitude, eccentricity, inclination, right_ascension, perigee, ta, beam,
                 focus="earth", rads=True):
        self._name = name
        self._altitude = altitude
        self._focus = focus
        self._true_alt = self.altitude + self.__get_radius()
        self._eccentricity = eccentricity
        self._beam = beam

        if not rads:
            self.inclination = inclination
            self.right_ascension = right_ascension
            self.perigee = perigee
            self.ta = ta
            self.inclination_r, self.right_ascension_r, self.perigee_r, self.ta_r = self.__convert_to_rads()
        else:
            self.inclination_r = inclination
            self.right_ascension_r = right_ascension
            self.perigee_r = perigee
            self.ta_r = ta
            self.inclination, self.right_ascension, self.perigee, self.ta = self.__convert_to_degs()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def altitude(self):
        return self._altitude

    @altitude.setter
    def altitude(self, new_alt):
        if new_alt < 100:
            return ValueError("Satellite's orbital altitude must be over the Karman line.")
        else:
            self._altitude = new_alt
            self._true_alt = new_alt + self.__get_radius()

    @property
    def true_alt(self):
        return self._true_alt

    @property
    def eccentricity(self):
        return self._eccentricity

    @eccentricity.setter
    def eccentricity(self, new_e):
        if new_e < 0:
            return ValueError("Eccentricity can't be set below a perfect circle.")
        else:
            self._eccentricity = new_e

    @property
    def beam(self):
        return self._beam

    @beam.setter
    def beam(self, new_beam):
        if new_beam < 0:
            return ValueError("Beam width must be between 0 and 180 degrees")
        elif new_beam > 180:
            return ValueError("Beam width must be between 0 and 180 degrees")
        else:
            self._beam = new_beam

    def __convert_to_rads(self, value=None):
        to_rad = pi / 180
        if value:
            return value * to_rad
        else:
            return self.inclination * to_rad, self.right_ascension * to_rad, self.perigee * to_rad, self.ta * to_rad

    def __convert_to_degs(self, value=None):
        to_deg = 180 / pi
        if value:
            return value * to_deg
        else:
            return self.inclination_r * to_deg, self.right_ascension_r * to_deg, self.perigee_r * to_deg, \
                   self.ta_r * to_deg

    def __get_radius(self):
        return heavenly_body_radius[self._focus.lower()]

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}".format(self.name, self.altitude, self.eccentricity,
                                                          self.inclination, self.right_ascension, self.perigee, self.ta)

    def __str__(self):
        return "Satellite Name: {0}, Alt: {1}, e: {2}, " \
               "Inclination: {3}, RA: {4}, Periapsis: {5}, Anomaly: {6}".format(self.name, self.altitude,
                                                                                self.eccentricity, self.inclination,
                                                                                self.right_ascension, self.perigee,
                                                                                self.ta)

    def as_dict(self, rads=True):
        if rads:
            sat = {"Name": self.name,
                   "Orbital Elements": {
                       "Eccentricity": self.eccentricity,
                       "Right Ascension": self.right_ascension_r,
                       "Semi-major Axis": self.true_alt,
                       "Arg. Periapsis": self.perigee_r,
                       "Mean Anomaly": self.ta_r,
                       "Inclination": self.inclination_r
                   },
                   "Beam Width": self.beam}
        else:
            sat = {"Name": self.name,
                   "Orbital Elements": {
                       "Eccentricity": self.eccentricity,
                       "Right Ascension": self.right_ascension,
                       "Semi-major Axis": self.true_alt,
                       "Arg. Periapsis": self.perigee,
                       "Mean Anomaly": self.ta,
                       "Inclination": self.inclination
                   },
                   "Beam Width": self.beam}
        sat['Focus'] = self._focus
        sat['Type'] = 'satellite'

        return sat

    def as_xml(self, epoch_date='2017-Jan-18 00:00:00', fov=1):
        warnings.warn("XML support is depreciated and not supported from PIGI 0.8.5 onward", DeprecationWarning)
        return '\t\t< Entity Type = "Satellite" Name = "{0}" >\n' \
               '\t\t\t<PropertySection Name="UserProperties">\n' \
               '\t\t\t\t<StringPropertyValue name="PlanetName" value="Earth"/>\n' \
               '\t\t\t\t<StringPropertyValue name="CatalogName" value="{0}"/>\n' \
               '\t\t\t\t<StringPropertyValue name="MeshName" value="SaberBox.mesh"/>\n' \
               '\t\t\t\t<StringPropertyValue name="BindingsFile" value=""/>\n' \
               '\t\t\t\t<IntPropertyValue name="ManualOrbitalElements" value="0"/>\n' \
               '\t\t\t\t<StringPropertyValue name="AssemblyFile" value=""/>\n' \
               '\t\t\t\t<StringPropertyValue name="SystemMapSourceId" value=""/>\n' \
               '\t\t\t\t<StringPropertyValue name="ResourceGroup" value="Autodetect"/>\n' \
               '\t\t\t\t<StringPropertyValue name="MetricSourceIds" value=""/>\n' \
               '\t\t\t\t<FloatPropertyValue name="BeamWidth" value="{1}"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="SGP4 Parameters">\n' \
               '\t\t\t\t<FloatPropertyValue name="B Star" value="-1.1606e-005"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Eccentricity" value="{2}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="RAAN" value="{3}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Semi-major axis" value="{4}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Arg. Perigee" value="{5}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Mean Anomaly" value="{6}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Inclination" value="{7}"/>\n' \
               '\t\t\t\t<TimestampPropertyValue name="Epoch" value="{8}"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Fov">\n' \
               '\t\t\t\t<EnumPropertyValue name="Enabled" value="{9}"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="TargetSceneNode">\n' \
               '\t\t\t\t<ArrayPropertyValue name="Position" value="[0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Orientation" value="[1, 0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Scale" value="[200, 200, 200]"/>\n' \
               '\t\t\t\t<EnumPropertyValue name="Debug" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Billboard">\n' \
               '\t\t\t\t<ArrayPropertyValue name="Position" value="[0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Orientation" value="[1, 0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Scale" value="[200, 200, 200]"/>\n' \
               '\t\t\t\t<EnumPropertyValue name="Debug" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Favourite">\n' \
               '\t\t\t\t<EnumPropertyValue name="favourite" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t</Entity>\n'.format(self.name, self.beam, self.eccentricity, self.right_ascension_r,
                                        self.true_alt, self.perigee_r, self.ta_r, self.inclination_r,
                                        epoch_date, fov)
