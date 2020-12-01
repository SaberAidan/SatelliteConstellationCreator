from __init__ import add_parent_directory  # Needed for running tests in pipelines

add_parent_directory()
import unittest
from satellite_constellation.Satellite import Satellite
from satellite_constellation.Constellation import *
from satellite_constellation.ConstellationExceptions import *
from satellite_constellation.SceneCreator import *
import math


class TestDummy(unittest.TestCase):
    flower_dummy = FlowerConstellation(3, 1, 4, 1, 4, 10, 10, 600)

    print(flower_dummy.representation(), flower_dummy.representation(representation_type="walker"))
    print("Np : {0}, Nd : {1}, Ns : {2} , Fn : {3}, Fd : {4}".format(flower_dummy.num_petals,
                                                                     flower_dummy.num_days,
                                                                     flower_dummy.num_satellites,
                                                                     flower_dummy.phasing_n,
                                                                     flower_dummy.phasing_d))
    print("Max sats per obit : {0}, Max sats : {1}, Num orbits : {2}, Num sats : {3}".format(
        flower_dummy.max_sats_per_orbit,
        flower_dummy.max_sats,
        flower_dummy.num_orbits,
        flower_dummy.num_satellites))

    print("RAAN spacing : {0}, Mean Anomaly Spacing : {1}".format(flower_dummy.raan_spacing,
                                                                  flower_dummy.mean_anomaly_spacing))

    print('Eccentricity : {0}, Orbital Period : {1}, Semi Major {2}'.format(flower_dummy.eccentricity,
                                                                            flower_dummy.orbital_period,
                                                                            flower_dummy.semi_major))

    print("RAAN : ", flower_dummy.raan)
    print("Mean Anomaly : ", flower_dummy.mean_anomaly)
    print("True Anomaly : ", flower_dummy.true_anomaly)
    print("Revisit Time : {0}, Minimum Revisit Time : {1}".format(flower_dummy.revisit_time,
                                                                  flower_dummy.minimum_revisit_time))

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

    # def test_focus(self):
    #     T, P, F = 2, 2, 1
    #     with self.assertRaises(FocusError):  # Beam width > 180
    #         constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [50], focus="notafocus")


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


if __name__ == '__main__':

    unittest.main()

