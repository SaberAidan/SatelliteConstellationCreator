"""
Recreation of the Octave SceneCreator script for use within Flask or other python projects
@creator: Aidan O'Brien
"""

from .Constellation import Constellation, SOCConstellation
from .GroundStation import GroundStation
from itertools import zip_longest
import warnings
from .utils import mod, heavenly_body_radius
from .ConstellationExceptions import *

def scene_xml_generator(scene):
    warnings.warn("XML support is depreciated and not supported from PIGI 0.8.5 onward", DeprecationWarning)
    scene_start = '<Pigi>\n' \
                  '\t<Entities>\n' \
                  '\t\t<Entity Type="Planet" Name="Mercury">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Mercury"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="0.1"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="0.244"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="0"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Venus">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Venus"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="177"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="0.605"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="1"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Earth">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Earth"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="23.4"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="0.6371"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="1"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Mars">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Mars"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="25"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="0.3371"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="1"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Jupiter">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Jupiter"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="3"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="6.99"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="0"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Saturn">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Saturn"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="27"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="6.033"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="0"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Neptune">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Neptune"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="30"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="2.46"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="0"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n' \
                  '\t\t<Entity Type="Planet" Name="Pluto">\n' \
                  '\t\t\t<PropertySection Name="UserProperties">\n' \
                  '\t\t\t\t<StringPropertyValue name="PlanetName" value="Pluto"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="obliquity" value="120"/>\n' \
                  '\t\t\t\t<FloatPropertyValue name="scale" value="0.186"/>\n' \
                  '\t\t\t\t<IntPropertyValue name="Atmosphere" value="0"/>\n' \
                  '\t\t\t</PropertySection>\n' \
                  '\t\t</Entity>\n'

    scene_end = '\t</Entities>\n' \
                '</Pigi>'

    scene_xml = ""
    for item in scene:
        scene_xml += item.as_xml()

    return "{0} {1} {2}".format(scene_start, scene_xml, scene_end)

def constellation_creator_general(num_constellations,constellation_types, constellation_vars):

    """
    :param num_constellations: Integer of the number of constellations that are for the scene
    :param constellation_types: List of types of constellation, i.e. "walker", "streets" or "flower"
    :param constellation_vars: List of dicts of constellation parameters

    """
    scene = []

    for idx in range(num_constellations):

        vars = constellation_vars[idx]

        if constellation_types[idx] == "walker":

            if not single_walker_errors(vars["satellite_nums"], vars["satellite_planes"], vars["plane_phasing"],

                             vars["inclination"], vars["altitude"], vars["eccentricity"], vars['constellation_beam_width'],
                             vars["focus"]):

                test_constellation = Constellation(vars["satellite_nums"],vars["satellite_planes"], vars["plane_phasing"],vars["inclination"],
                                                  vars["altitude"], vars["eccentricity"], vars['constellation_beam_width'], vars["focus"])

        elif constellation_types[idx] == "streets":

            temp_constellation = SOCConstellation(vars["num_streets"], vars["street_width"], vars["altitude"],
                             vars["beam_width"], vars["raan"], vars["eccentricity"], vars['revisit_time'], vars["focus"])

def single_walker_errors( satellite_nums, satellite_planes, plane_phasing, inclination, altitude, eccentricity,
                  constellation_beam_width, sat_name="Sat", focus="earth"):

    error_found = False

    if satellite_nums%satellite_planes:
        error_found = True
        raise ConstellationPlaneMismatchError("Number of satellites not compatible with planes in constellation")

    if plane_phasing < 0:
        error_found = True
        raise PhaseError("Negative number of phases not allowed")
    # elif phase >= satellite_planes[idx]:
    #     raise PhaseError("Number of phases in constellation larger than number of planes")
    # elif (satellite_nums[idx] % phase) != 0:
    #     raise PhaseError("Number of satellites must be divisible by the number of phases")
    elif plane_phasing % 1:
        error_found = True
        raise PhaseError("Number of phases not an integer")
    elif plane_phasing % satellite_planes == 0:
        error_found = True
        raise PhaseError("Phase cannot equal number of planes times a constant")

    if abs(inclination) > 90:
        error_found = True
        raise InclinationError("Inclination greater than 90 degrees from equitorial orbit")
        #Could implement wrap around of angles here

    if altitude < 0:
        error_found = True
        raise AltitudeError("Negative altitude not allowed")
    elif (altitude < 100):
        error_found = True
        raise AltitudeError("Altitude below Karman line")

    if (eccentricity < 0):
        error_found = True
        raise EccentricityError("Invalid eccentricity. Eccentricity cannot be less than 0")
    elif (eccentricity >= 1):
        error_found = True
        raise EccentricityError("Invalid eccentricity. Hyperbolic trajectories not allowed")

    if constellation_beam_width <= 0 or constellation_beam_width > 180:
        error_found = True
        raise BeamError("Beam width error. Must be within 0 to 180 degrees")

    if focus.lower() not in heavenly_body_radius:
        error_found = True
        raise FocusError("'" + focus.capitalize() + "' not supported as a celestial body origin.")

    if error_found:
        return True
    else:
        return False

def constellation_creator(num_constellations, satellite_nums, satellite_planes, plane_phasing, inclination, altitude,
                          eccentricity, constellation_beam_width, sat_name="Sat", focus="earth"):
    """

    Need to add the streets of coverage method to the scene creator
    Change this to use 3 param
    Number of constellations
    List of constellation types
    List of dict with constellation vars, dict change with different type

    :param num_constellations: Integer of the number of constellations that are for the scene
    :param satellite_nums: List of numbers of satellites for each constellation
    :param satellite_planes: List of number of planes of satellites for each constellation.
    :param plane_phasing: The planar phases for their respective constellations
    :param inclination: List of inclinations for constellations
    :param altitude: List of altitudes for constellations
    :param eccentricity: List of eccentricities for each constellation
    :param constellation_beam_width: List of Beam widths for satellites in each constellation
    :param sat_name: Root name for satellites, defaults to Sat
    :param focus: The focus of the satellite, i.e. Earth, Luna, Mars, in lower case.
    :return: Returns a list of constellations that have been formatted
    """

    walker_errors(num_constellations, satellite_nums, satellite_planes, plane_phasing, inclination, altitude,
                          eccentricity, constellation_beam_width,focus = "earth")

    scene = []
    for idx in range(num_constellations):
        start_num = sum(satellite_nums[0:idx])
        scene.append(Constellation(satellite_nums[idx], satellite_planes[idx], plane_phasing[idx], inclination[idx],
                                   altitude[idx], eccentricity[idx], constellation_beam_width[idx], name=sat_name,
                                   starting_number=start_num, focus=focus))

    return scene

def ground_array_creator(num_ground_stations, latitudes, longitudes, elevations, beam_widths, name="GS"):
    """

    :param num_ground_stations:
    :param latitudes:
    :param longitudes:
    :param elevations:
    :param beam_widths:
    :param name:
    :return:
    """

    if num_ground_stations < 1 or (num_ground_stations % 1):
        return -1

    if num_ground_stations != len(latitudes):
        return -1

    try:
        zip_longest(latitudes, longitudes, elevations)
    except IndexError:
        return -1

    if (len(latitudes) != len(beam_widths)) or (len(beam_widths) == 1):
        return -1

    scene = []
    for idx in range(num_ground_stations):
        gs_name = name + str(idx)
        if len(beam_widths) > 1:
            beam = beam_widths[idx]
        else:
            beam = beam_widths
        scene.append(GroundStation(gs_name, latitudes[idx], longitudes[idx], elevations[idx], beam))

    return scene

def walker_errors(num_constellations, satellite_nums, satellite_planes, plane_phasing, inclination, altitude, eccentricity,
                  constellation_beam_width, sat_name="Sat", focus="earth"):

    error_found = False

    if num_constellations < 1 or (num_constellations % 1):
        error_found = True
        raise ConstellationNumberError("Invalid integer number of constellations")

    if num_constellations != len(satellite_nums):
        error_found = True
        raise ConstellationConfigurationError("Number of constellations does not match number of satellites provided")

    if any(mod(satellite_nums, satellite_planes)):
        error_found = True
        raise ConstellationPlaneMismatchError("Number of satellites not compatible with planes in constellation")

    for idx, phase in enumerate(plane_phasing):
        if phase < 0:
            error_found = True
            raise PhaseError("Negative number of phases not allowed")
        # elif phase >= satellite_planes[idx]:
        #     raise PhaseError("Number of phases in constellation larger than number of planes")
        # elif (satellite_nums[idx] % phase) != 0:
        #     raise PhaseError("Number of satellites must be divisible by the number of phases")
        elif phase % 1:
            error_found = True
            raise PhaseError("Number of phases not an integer")
        elif phase % satellite_planes[idx] == 0:
            error_found = True
            raise PhaseError("Phase cannot equal number of planes times a constant")

    if any(abs(x) > 90 for x in inclination):
        error_found = True
        raise InclinationError("Inclination greater than 90 degrees from equitorial orbit")
        #Could implement wrap around of angles here

    if any(x < 0 for x in altitude):
        error_found = True
        raise AltitudeError("Negative altitude not allowed")
    elif any(x < 100 for x in altitude):
        error_found = True
        raise AltitudeError("Altitude below Karman line")

    if any(x < 0 for x in eccentricity):
        error_found = True
        raise EccentricityError("Invalid eccentricity. Eccentricity cannot be less than 0")
    elif any(x >= 1 for x in eccentricity):
        error_found = True
        raise EccentricityError("Invalid eccentricity. Hyperbolic trajectories not allowed")

    if any(x <= 0 or x > 180 for x in constellation_beam_width):
        error_found = True
        raise BeamError("Beam width error. Must be within 0 to 180 degrees")

    if focus.lower() not in heavenly_body_radius:
        error_found = True
        raise FocusError("'" + focus.capitalize() + "' not supported as a celestial body origin.")

    if error_found:
        return True
    else:
        return False

def create_scene(num_constellations, num_sats, sat_planes, plane_phasing, sat_inclination, sat_alt,
                 sat_eccentricity, const_beam_width,
                 num_ground_stations, latitudes, longitudes, elevations, beam_widths,
                 sat_name="Sat", gsname="GS"):
    """
    Creates the scene list from constellations and ground stations
    :param num_constellations:
    :param num_sats:
    :param sat_planes:
    :param plane_phasing:
    :param sat_inclination:
    :param sat_altitude:
    :param sat_eccentricity:
    :param const_beam_width:
    :param num_ground_stations:
    :param latitudes:
    :param longitudes:
    :param elevations:
    :param beam_widths: List of beam widths for the ground stations
    :param sat_name: Root of satellite names
    :param gsname: Root of ground station names
    :return: List of constellation and ground station objects
    """

    scene = constellation_creator(num_constellations, num_sats, sat_planes, plane_phasing, sat_inclination, sat_alt,
                                  sat_eccentricity, const_beam_width, sat_name=sat_name)
    scene.extend(ground_array_creator(num_ground_stations, latitudes, longitudes, elevations, beam_widths, name=gsname))

    return scene