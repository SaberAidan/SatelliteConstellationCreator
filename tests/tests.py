import unittest
from satellite_constellation.Satellite import Satellite
from satellite_constellation.Constellation import SOCConstellation
from satellite_constellation.ConstellationExceptions import *
from satellite_constellation.SceneCreator import constellation_creator
import warnings


class testSatellite(unittest.TestCase):  # Test for errors in the constellation creator

    def setUp(self):
        testCons = SOCConstellation(1, 10, 1500, 60, [20], 0.8, 100)
        # self.constellation = constellation_creator(1, [18], [3], [1], [30], [1000], [0], [20])

    def test_0_constellations(self):  # Check for number error with 0 constellations
        constellation_num = 0
        T, P, F = 1, 1, 1
        with self.assertRaises(ConstellationNumberError) as error:
            constellation = constellation_creator(constellation_num, [T], [P], [F], [30], [1000], [0.5], [20])

    def test_constellation_mismatch(self):  # Check for mismatch error catching with an empty satellite nums list
        constellation_num = 1
        T, P, F = [18, 18], [3, 3], [1, 1]
        satellite_nums = [T]
        with self.assertRaises(ConstellationConfigurationError) as error:
            constellation = constellation_creator(constellation_num, T, P, F, [30], [1000], [0.5], [20])

    def test_plane_mismatch(
            self):  # Check for mismatch error if the satellites cannot be evenly spread across the planes
        T, P, F = 1, 2, 1
        plane_nums = [P]
        satellite_nums = [T]
        with self.assertRaises(ConstellationPlaneMismatchError):
            constellation = constellation_creator(1, satellite_nums, plane_nums, [F], [30], [1000], [0.5], [20])

    def test_phasing(self):
        T, P, F = 18, 3, 3  # This would put the satellites right on top of eachother
        with self.assertRaises(PhaseError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [20])
        T, P, F = 18, 3, 6  # This would also put the satellites right on top of eachother
        with self.assertRaises(PhaseError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [20])

    def test_inclination(self):
        T, P, F = 2, 2, 1
        inclination = 100
        with self.assertRaises(InclinationError):
            constellation = constellation_creator(1, [T], [P], [F], [inclination], [1000], [0.5], [20])

    def test_altitude(self):
        T, P, F = 2, 2, 1
        altitude = -1
        with self.assertRaises(AltitudeError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [altitude], [0.5], [20])
        altitude = 50
        with self.assertRaises(AltitudeError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [altitude], [0.5], [20])

    def test_eccentricity(self):
        T, P, F = 2, 2, 1
        eccentricity = -1
        with self.assertRaises(EccentricityError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [eccentricity], [20])
        eccentricity = 2
        with self.assertRaises(EccentricityError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [eccentricity], [20])

    def test_beam_width(self):
        T, P, F = 2, 2, 1
        beam_width = 190
        with self.assertRaises(BeamError):  # Beam width > 180
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [beam_width])
        T, P, F = 2, 2, 1
        beam_width = 0
        with self.assertRaises(BeamError):  # Beam width = 0
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [beam_width])

    def test_focus(self):
        T, P, F = 2, 2, 1
        with self.assertRaises(FocusError):  # Beam width > 180
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [50], focus="notafocus")

class testStreets(unittest.TestCase):

    def setUp(self):
        self.streets_constellation = testCons = SOCConstellation(1, 10, 1500, 60, [20], 0.8, 100)

    def test_perigee(self):
        self.assertAlmostEqual(self.streets_constellation.perigee,7871,3)

    def test_semi_major(self):
        self.assertAlmostEqual(self.streets_constellation.semi_major,39355,3)

    def test_orbital_period(self):
        self.assertAlmostEqual(self.streets_constellation.orbital_period, 77700.55,1)

    def test_earth_radial_coverage(self):
        self.assertAlmostEqual(self.streets_constellation.earth_coverage_radius, 868.7 , 1)

    def test_earth_angular_coverage(self):
        self.assertAlmostEqual(self.streets_constellation.earth_coverage_angle, 0.136, 3)

    def test_required_satellites_by_coverage(self):
        self.assertAlmostEqual(self.streets_constellation.num_sats, 30)

    def test_required_satellites_by_period(self):
        streets_constellation = SOCConstellation(1, 10, 1500, 60, [20], 0.8, 7770)
        self.assertAlmostEqual(streets_constellation.num_sats,10,1)




if __name__ == '__main__':
    unittest.main()
