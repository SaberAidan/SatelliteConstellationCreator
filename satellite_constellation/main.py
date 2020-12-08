from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import plotly as ply
import plotly.graph_objects as go

if __name__ == '__main__':
    myWalker = WalkerConstellation(9, 3, 1, 60, 36000, 0, 45)
    draw_walker_plotly(myWalker, links=True)
    print(myWalker.average_links)
    myFlower = FlowerConstellation(37, 18, 57, 6, 19, 0, 0, 19702, 0)
    draw_flower_plotly(myFlower)
    print(myFlower.average_links)
    myStreets = SOCConstellation(6, 15, 15000, 60, [0, 30, 60, 90, 120, 150], 0, 1000)
    draw_soc_plotly(myStreets, links = True)
    print(myStreets.average_links)
