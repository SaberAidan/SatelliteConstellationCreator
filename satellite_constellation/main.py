from satellite_constellation.visualiser import *
import matplotlib.pyplot as plt
import argparse

def argParser():
    """
    Args:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', dest='test',
                        action="store_true",default=False,
                        help="Run test constellations")
    return parser.parse_args()

def main():
    p = argParser()
    if p.test:
        draw_walker(WalkerConstellation(8, 2, 1, 45, 1000, 0, 20))
        draw_flower(FlowerConstellation(37, 18, 57, 6, 19, 0, 0, 19702, 30))
        plt.show()
    else:
        print('Nothing to do') #TODO add more functionality here

if __name__ == '__main__':
    main()
