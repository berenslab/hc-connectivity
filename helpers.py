import matplotlib
import seaborn as sns

def sns_styleset():
    sns.set_context('paper')
    sns.set_style('ticks')
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
    matplotlib.rcParams['mathtext.rm'] = 'Arial'
    matplotlib.rcParams['axes.linewidth']    = .75
    matplotlib.rcParams['lines.linewidth']   = .75
    matplotlib.rcParams['xtick.major.width'] = .75
    matplotlib.rcParams['ytick.major.width'] = .75
    matplotlib.rcParams['xtick.major.size'] = 3
    matplotlib.rcParams['ytick.major.size'] = 3
    matplotlib.rcParams['xtick.minor.size'] = 2
    matplotlib.rcParams['ytick.minor.size'] = 2
    matplotlib.rcParams['font.size']       = 9
    matplotlib.rcParams['axes.titlesize']  = 9
    matplotlib.rcParams['axes.labelsize']  = 9
    matplotlib.rcParams['legend.fontsize'] = 8
    matplotlib.rcParams['xtick.labelsize'] = 9
    matplotlib.rcParams['ytick.labelsize'] = 9
