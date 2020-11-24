import unittest
from satellite_constellation.Satellite import Satellite
from satellite_constellation.Constellation import Constellation
from satellite_constellation.ConstellationExceptions import *
from satellite_constellation.SceneCreator import constellation_creator
import warnings

class testSatellite(unittest.TestCase):

    # def setUp(self):
    #     self.constellation_creator = constellation_creator(1,[1],[1],[1],[30],[500],[0.5],[20])

    def test_0_constellations(self): #Check for number error with 0 constellations
        constellation_num = 0
        T, P, F = 1, 1, 1
        with self.assertRaises(ConstellationNumberError) as error:
            self.constellation = constellation_creator(constellation_num,[T],[P],[F],[30],[1],[0.5],[20])

    def test_constellation_mismatch(self): #Check for mismatch error catching with an empty satellite nums list
        constellation_num = 1
        T, P , F = [18,18],[3,3],[1,1]
        satellite_nums = [T]
        with self.assertRaises(ConstellationConfigurationError) as error:
            self.constellation = constellation_creator(constellation_num,T,P,F,[30],[1000],[0.5],[20])

    def test_plane_mismatch(self): #Check for mismatch error if the satellites cannot be evenly spread across the planes
        T, P, F = 1, 2, 1
        plane_nums = [P]
        satellite_nums = [T]
        with self.assertRaises(ConstellationPlaneMismatchError):
            self.constellation = constellation_creator(1, satellite_nums, plane_nums, [F], [30], [1000], [0.5], [20])

    def test_phasing(self):
        T, P, F = 18, 3, 3 #This would put the satellites right on top of eachother
        with self.assertRaises(PhaseError):
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [20])
        T, P, F = 18, 3, 6 #This would also put the satellites right on top of eachother
        with self.assertRaises(PhaseError):
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [20])

    def test_inclination(self):
        T, P, F = 2, 2, 1
        inclination = 100
        with self.assertRaises(InclinationError):
            self.constellation = constellation_creator(1, [T], [P], [F], [inclination], [1000], [0.5], [20])

    def test_altitude(self):
        T, P, F = 2, 2, 1
        altitude = -1
        with self.assertRaises(AltitudeError):
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [altitude], [0.5], [20])
        altitude = 50
        with self.assertRaises(AltitudeError):
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [altitude], [0.5], [20])

    def test_eccentricity(self):
        T, P, F = 2, 2, 1
        eccentricity = -1
        with self.assertRaises(EccentricityError):
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [eccentricity], [20])
        eccentricity = 2
        with self.assertRaises(EccentricityError):
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [eccentricity], [20])

    def test_beam_width(self):
        T, P, F = 2, 2, 1
        beam_width = 190
        with self.assertRaises(BeamError): # Beam width > 180
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [beam_width])
        T, P, F = 2, 2, 1
        beam_width = 0
        with self.assertRaises(BeamError): # Beam width = 0
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [beam_width])

    def test_focus(self):
        T, P, F = 2, 2, 1
        with self.assertRaises(FocusError): # Beam width > 180
            self.constellation = constellation_creator(1, [T], [P], [F], [30], [1000], [0.5], [50], focus = "notafocus")

if __name__ == '__main__':
    unittest.main()









