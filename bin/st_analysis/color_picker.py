from matplotlib import cm
import matplotlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cmap', type=str, default='Spectral')
parser.add_argument('-n', '--ncolors', type=int, default=20)
args = parser.parse_args()

cmap = cm.get_cmap(args.cmap, args.ncolors)    # PiYG

for i in range(cmap.N-1):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    print("'", matplotlib.colors.rgb2hex(rgba), "'", sep='', end=', ')
for i in range(cmap.N-1, cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    print("'", matplotlib.colors.rgb2hex(rgba), "'", sep='')