from satellite_constellation.Constellation import *
from plotter import visualiser
from orbital_analyzer.analysis import *
from satellite_constellation.Errors import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(4, 2, 1, 45, 35785, 0, 5)
    visualiser.draw_walker_plotly(myWalker, sensor_regions=True)

    myAnalyzer = orbital_analysis(myWalker)

    # start = time.process_time()
    # myAnalyzer.find_links_np()
    # print(myAnalyzer.calculate_constellation_coverage_np(resolution=1))
    # print(myAnalyzer.ground_station_visibility(0,0,5))
    # print("NP revisit calc time", time.process_time() - start, 's')

    myAnalyzer.check_ground_FOV_np(myAnalyzer.constellation.as_numpy(), 0, 0, 20,0)

    # myStreets = SOCConstellation(6, 5, 38758, 5, [0, 30, 60, 90, 120, 150], 0, 1)
    # visualiser.draw_soc_plotly(myStreets, sensor_regions=True)
    #
    # myAnalyzer = orbital_analysis(myStreets)
    #
    # start = time.process_time()
    # # myAnalyzer.find_links_np()
    # print(myAnalyzer.calculate_constellation_coverage_np(resolution=1))
    # print("NP revisit calc time", time.process_time() - start, 's')

    # start = time.process_time()
    # # myAnalyzer.find_links_np()
    # print(myAnalyzer.calculate_constellation_coverage(resolution=1))
    # print("Normal revisit calc time", time.process_time() - start, 's')

    # start = time.process_time()
    # myAnalyzer.find_links()
    # print("Normal revisit calc time", time.process_time() - start, 's')
