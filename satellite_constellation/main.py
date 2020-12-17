from satellite_constellation.Constellation import *
from plotter import visualiser
from orbital_analyzer.analysis import *
from satellite_constellation.Errors import *
import time
import json

if __name__ == '__main__':
    myWalker = WalkerConstellation(30, 15, 1, 60, 35786, 0, 20)
    visualiser.draw_walker_plotly(myWalker, sensor_regions=True)

    myAnalyzer = orbital_analysis(myWalker)

    start = time.process_time()
    res = myAnalyzer.revisit_numpy(0, 0)
    print("\nNP revisit calc time", time.process_time() - start, 's')
    for key in res:
        print(key, ':', res[key])

    start = time.process_time()
    res = myAnalyzer.ground_station_visibility_np(0, 0, 40)
    print("\nNP ground revisit calc time", time.process_time() - start, 's')
    for key in res:
        print(key, ':', res[key])

    start = time.process_time()
    res = myAnalyzer.calculate_constellation_coverage_np(resolution=1)
    print("\nNP coverage calc time", time.process_time() - start, 's')
    print(res)

    start = time.process_time()
    res = myAnalyzer.find_links_np(60, 6)
    print("\n NP links calc time", time.process_time() - start, 's')
    for key in res:
        print(key, ':', res[key])
