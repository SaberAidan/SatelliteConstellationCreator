from satellite_constellation.Constellation import *
from plotter import visualiser
from orbital_analyzer.analysis import *
from satellite_constellation.Errors import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(40, 2, 1, 45, 35785, 0, 5)
    # visualiser.draw_walker_plotly(myWalker, sensor_regions=True)

    myAnalyzer = orbital_analysis(myWalker)

    start = time.process_time()
    res = myAnalyzer.ground_station_visibility_np(0, 0, 20)
    for key in res:
        print(key, ':', res[key])
    print("NP station revisit calc time", time.process_time() - start, 's \n')

    # start = time.process_time()
    # print(myAnalyzer.ground_station_visibility(0, 0, 20))
    # print("Normal station revisit calc time", time.process_time() - start, 's')

    start = time.process_time()
    res = myAnalyzer.revisit_numpy(0, 0)
    for key in res:
        print(key, ':', res[key])
    print("NP revisit calc time", time.process_time() - start, 's \n')

    # start = time.process_time()
    # print(myAnalyzer.calculate_revisit(0, 0))
    # print("Normal revisit calc time", time.process_time() - start, 's')
