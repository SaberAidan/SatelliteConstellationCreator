from satellite_constellation.Constellation import *
from plotter import visualiser
from orbital_analyzer.analysis import *
from satellite_constellation.Errors import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(3, 3, 1, 80, 35785, 0, 5)

    start = time.process_time()
    myWalker.as_cartesian_np()
    print("NP cartesian calc time", time.process_time() - start, 's')

    start = time.process_time()
    myWalker.as_cartesian()
    print("Normal cartesian calc time", time.process_time() - start, 's')

    start = time.process_time()
    myWalker.as_geographic_np()
    print("NP geo calc time", time.process_time() - start, 's')

    start = time.process_time()
    myWalker.as_geographic()
    print("Normal geo calc time", time.process_time() - start, 's')

    myAnalyzer = orbital_analysis(myWalker)

    start = time.process_time()
    print(myAnalyzer.calculate_revisit_numpyg(0,0))
    print("NP revisit calc time", time.process_time() - start, 's')

    start = time.process_time()
    print(myAnalyzer.calculate_revisit(0,0))
    print("Normal revisit calc time", time.process_time() - start, 's')

    # visualiser.draw_walker_plotly(myWalker, sensor_regions=True, links=False)

    # myFlower = FlowerConstellation(8, 1, 9, 1, 9, 0, 45, 2500, 10)
    # visualiser.draw_flower_plotly(myFlower, links=False)

    # myStreets = SOCConstellation(1, 20, 15000, 20, [0], 0, 10)
    # visualiser.draw_soc_plotly(myStreets, links=False, sensor_regions=True)

    # print(myWalker.as_spherical_np())
    # start = time.process_time()
    # myWalker.as_geographic_np()
    # print("NP geo calc time", time.process_time() - start, 's')

    # sat_to_xyz_np(np.array(myWalker.satellites))
    # print(myWalker.satellites
    # print(len(np.arange(0,math.pi/2,(math.pi/2)/3000)[0:3000]))
    # sat_to_xyz_np(myWalker.as_numpy())
    # myWalker.as_geographic()

    # angs = np.arange(0, math.pi / 2, (math.pi / 2) / myWalker.num_satellites)[0:myWalker.num_satellites]
    # rpy = np.column_stack([angs, angs, angs])
    # rotate_np(myWalker.as_cartesian(), angs, ax='custom', rpy=rpy, basis = rpy)
    # sat_to_xyz_np(myWalker.as_numpy())

    # print(myWalker.as_cartesian_np())
    # print(myWalker.as_cartesian())
    # print(abs(myWalker.as_cartesian_np() - myWalker.as_cartesian()))

