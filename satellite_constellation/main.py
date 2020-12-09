from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import plotly as ply
import plotly.graph_objects as go

if __name__ == '__main__':
    # myWalker = WalkerConstellation(8,  4, 3, 60, 36000, 0, 45)
    # draw_walker_plotly(myWalker, links=True, sensor_regions=True)
    # myFlower = FlowerConstellation(37, 18, 57, 6, 19, 0, 0, 19702, 0)
    # draw_flower_plotly(myFlower, links = False)
    myStreets = SOCConstellation(1, 60, 12000, 15, [0,  60,  120], 0, 10)
    draw_soc_plotly(myStreets, links=False, sensor_regions=True)
