from __init__ import add_parent_directory  # Needed for running tests in pipelines

if __name__ == "__main__":
    add_parent_directory()  # This must stay prior to the satellite_constellation imports
from satellite_constellation.Constellation import *
from satellite_constellation.SceneCreator import *
import unittest
import math


class TestConstellationCreator(unittest.TestCase):  # Test for errors in the constellation creator

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
        T, P, F = 18, 3, 3  # This would put the satellites right on top of each other
        with self.assertRaises(PhaseError):
            constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [20])
        T, P, F = 18, 3, 6  # This would also put the satellites right on top of each other
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


class TestSatellite(unittest.TestCase):

    def setUp(self):
        self.satellite = Satellite("testSat", 1000, 0, 10, 20, 30, 40, 30, rads=False)

    def test_true_altitude(self):
        self.assertEqual(self.satellite.true_alt, 7371)

    def test_deg_to_rad(self):
        self.assertAlmostEqual(0.1745, self.satellite.inclination_r, 4)
        self.assertAlmostEqual(0.3491, self.satellite.right_ascension_r, 4)
        self.assertAlmostEqual(0.5236, self.satellite.perigee_r, 4)
        self.assertAlmostEqual(0.6981, self.satellite.ta_r, 4)

    def test_rad_to_degree(self):
        rad_sat = Satellite("testSat", 1000, 0, math.pi / 5, math.pi / 4, math.pi / 3, math.pi / 2, math.pi / 5,
                            rads=True)
        self.assertAlmostEqual(36, rad_sat.inclination, 4)
        self.assertAlmostEqual(45, rad_sat.right_ascension, 4)
        self.assertAlmostEqual(60, rad_sat.perigee, 4)
        self.assertAlmostEqual(90, rad_sat.ta, 4)


class TestWalker(unittest.TestCase):

    def setUp(self):
        self.walker_constellation = WalkerConstellation(18, 3, 1, 30, 1000, 0, 20)

    def test_sats_per_plane(self):
        self.assertEqual(6, self.walker_constellation.sats_per_plane)

    def test_phasing(self):
        self.assertEqual(20, self.walker_constellation.correct_phasing)


class TestStreets(unittest.TestCase):

    def setUp(self):
        self.streets_constellation = SOCConstellation(1, 10, 1500, 60, [20], 0.8, 100)

    def test_perigee(self):
        self.assertAlmostEqual(self.streets_constellation.perigee, 7871, 3)

    def test_semi_major(self):
        self.assertAlmostEqual(self.streets_constellation.semi_major, 39355, 3)

    def test_orbital_period(self):
        self.assertAlmostEqual(self.streets_constellation.orbital_period, 77700.55, 1)

    def test_earth_radial_coverage(self):
        self.assertAlmostEqual(self.streets_constellation.earth_coverage_radius, 868.7, 1)

    def test_earth_angular_coverage(self):
        self.assertAlmostEqual(self.streets_constellation.earth_coverage_angle, 0.136, 3)

    def test_required_satellites_by_coverage(self):
        self.assertAlmostEqual(self.streets_constellation.num_satellites, 30)

    def test_required_satellites_by_period(self):
        streets_constellation = SOCConstellation(1, 10, 1500, 60, [20], 0.8, 7770)
        self.assertAlmostEqual(streets_constellation.num_satellites, 10, 1)


class TestFlower(unittest.TestCase):

    # Test cases from "The Flower Constellations - Theory, Design Process and Applications" - Wilkins, P.M
    def setUp(self):
        self.flower_suite = [FlowerConstellation(8, 1, 9, 1, 9, 0, 0, 2500, 0),
                             FlowerConstellation(769, 257, 4, 1, 4, 0, 0, 600, 0),
                             FlowerConstellation(4, 1, 4, 1, 4, 0, 0, 600, 0),
                             FlowerConstellation(3, 1, 4, 1, 4, 0, 0, 600, 0),
                             FlowerConstellation(3, 1, 4, 1, 7, 0, 0, 600, 0),
                             FlowerConstellation(3, 1, 4, 1, 2, 0, 0, 600, 0),
                             FlowerConstellation(3, 2, 4, 1, 2, 0, 0, 600, 0),
                             FlowerConstellation(31, 11, 30, 7, 10, 0, 0, 9000, 0),
                             FlowerConstellation(37, 18, 57, 6, 19, 0, 0, 19702, 0),
                             FlowerConstellation(15, 7, 49, 23, 49, 0, 0, 19702, 0)]

    def test_num_satellites(self):
        num_sats = []
        num_sats_test = [9, 4, 4, 4, 4, 2, 4, 30, 57, 49]
        for idx in range(len(self.flower_suite)):
            num_sats.append(self.flower_suite[idx].num_satellites)
        self.assertEqual(num_sats, num_sats_test)

    def test_max_satellites(self):
        max_sats = []
        max_sats_test = [9, 1028, 4, 4, 7, 2, 4, 110, 342, 343]
        for idx in range(len(self.flower_suite)):
            max_sats.append(self.flower_suite[idx].max_sats)
        self.assertEqual(max_sats, max_sats_test)

    def test_max_sats_per_orbit(self):
        max_sats_per_orbit = []
        max_sats_per_orbit_test = [1, 257, 1, 1, 1, 1, 2, 11, 18, 7]
        for idx in range(len(self.flower_suite)):
            max_sats_per_orbit.append(self.flower_suite[idx].max_sats_per_orbit)
        self.assertEqual(max_sats_per_orbit, max_sats_per_orbit_test)

    def test_raan_spacing(self):
        raan_spacing_test = [-40, -90, -90, -90, -51.42, -180, -180, -252, -113.68, -168.97]
        for idx in range(len(self.flower_suite)):
            self.assertAlmostEqual(raan_spacing_test[idx], self.flower_suite[idx].raan_spacing, delta=0.1)

    def test_mean_anomaly_spacing(self):
        mean_anomaly_spacing_test = [320, 269.2, 360, 270, 154.28, 180, 270, 350.18, 233.68, 2.099]
        for idx in range(len(self.flower_suite)):
            self.assertAlmostEqual(mean_anomaly_spacing_test[idx], self.flower_suite[idx].mean_anomaly_spacing,
                                   delta=0.1)

    def test_raan_anomaly_simple(self):
        param_list = []
        param_list_test = [[0, 0], [40, 40], [80, 80], [120, 120], [160, 160], [200, 200], [240, 240], [280, 280],
                           [320, 320]]
        for idx in range(len(self.flower_suite[0].raan)):
            param_list.append([self.flower_suite[0].raan[idx], self.flower_suite[0].mean_anomaly[idx]])

        self.assertTrue(len(param_list) == len(param_list_test))

        for idx in range(len(param_list)):
            self.assertTrue(param_list_test[idx] in param_list)

    def test_raan_anomaly_complex(self):
        param_list = []
        param_list_test = [[0, 0], [0, 180], [180, 90], [180, 270]]
        for idx in range(len(self.flower_suite[6].raan)):
            param_list.append([self.flower_suite[6].raan[idx], self.flower_suite[6].mean_anomaly[idx]])

        for idx in range(len(param_list)):
            self.assertTrue(param_list_test[idx] in param_list)


if __name__ == '__main__':
    unittest.main()
