from satellite_constellation.visualiser import *
from satellite_constellation.Constellation import *
import plotly as ply
import plotly.graph_objects as go

if __name__ == '__main__':
    myWalker = WalkerConstellation(8, 8, 2, 60, 36000, 0, 45)
    draw_walker_plotly(myWalker)
    myFlower = FlowerConstellation(4, 1, 4, 1, 4, 0, 30, 5000, 0)
    draw_flower_plotly(myFlower)
