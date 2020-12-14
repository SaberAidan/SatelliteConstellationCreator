from satellite_constellation.ConstellationExceptions import *
from .utils import *


# General constellation errors
def constellation_errors(num_satellites, orbital_period, altitude, beam_width, eccentricity, inclination, focus='Earth',
                         name="constellation"):
    if altitude < 0:
        raise AltitudeError("Negative altitude not allowed")
    elif altitude < 100:
        raise AltitudeError("Altitude below Karman line")

    if eccentricity < 0:
        raise EccentricityError("Invalid eccentricity. Eccentricity cannot be less than 0")
    elif eccentricity >= 1:
        raise EccentricityError("Invalid eccentricity. Hyperbolic trajectories not allowed")

    if abs(inclination) > 90:
        raise InclinationError("Inclination greater than 90 degrees from equatorial orbit")

    if beam_width <= 0 or beam_width > 180:
        raise BeamError("Beam width error. Must be within 0 to 180 degrees")

    if focus not in heavenly_body_radius:
        raise FocusError("'" + focus.capitalize() + "' not supported as a celestial body origin.")


# Walker specific constellation errors
def walker_errors(satellite_nums, satellite_planes, plane_phasing, inclination, altitude,
                  eccentricity,
                  constellation_beam_width, sat_name="Sat", focus="earth"):
    if satellite_nums % satellite_planes:
        raise ConstellationPlaneMismatchError("Number of satellites not compatible with planes in constellation")

    if plane_phasing < 0:
        raise PhaseError("Negative number of phases not allowed")
    if plane_phasing % 1:
        raise PhaseError("Number of phases not an integer")
    if plane_phasing % satellite_planes == 0:
        raise PhaseError("Phase cannot equal number of planes times an integer constant")


# Flower specific constellation errors
def flower_errors(num_petals, num_days, num_satellites, phasing_n, phasing_d, perigee_argument,
                  inclination, perigee_altitude, beam_width, focus='earth', name="constellation"):
    if num_satellites > phasing_d * num_days :
        raise MaxSatellitesExceededError("Number of satellites specified greater than maximum possible: {0}".format(phasing_d*num_days))

    if num_petals < 1 or (num_petals % 1):
        raise ConstellationNumberError("Invalid integer number of petals")

    if num_days < 1 or (num_days % 1):
        raise ConstellationNumberError("Invalid integer number of days to repeat")

    if num_satellites < 1 or (num_satellites % 1):
        raise ConstellationNumberError("Invalid integer number of satellites")

    if phasing_n < 1 or (phasing_n % 1):
        raise ConstellationNumberError("Invalid integer number of phasing_n")

    if phasing_d < 1 or (phasing_d % 1):
        raise ConstellationNumberError("Invalid integer number of phasing_d")


# Street specific constellation errors
def street_errors(num_streets, street_width, altitude, beam_width, raan, eccentricity, revisit_time,
                  name="Streets", focus="earth", starting_number=0):

    if num_streets != len(raan):
        raise ConstellationConfigurationError("Number of streets not compatible with number of right ascensions "
                                              "provided")

    if num_streets < 1 or (num_streets % 1):
        raise ConstellationNumberError("Invalid integer number of streets")

    if street_width > 180 or street_width <= 0:
        raise StreetWidthError("Street width error. Must be within 0 to 180 degrees")

    if revisit_time < 0:
        raise RevisitError("Revisit time cannot be less than 0 s")
