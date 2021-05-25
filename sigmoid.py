# %%
from pylab import np, plt, mpl
# # 字体
# mpl.rcParams['mathtext.fontset'] = 'stix'
# mpl.rcParams['font.family'] = 'STIXGeneral'
# plt.rcParams['font.weight'] = '900'
# plt.rcParams['font.size'] = 12  # 13 is too big, commenting out is too small

# # 颜色设置壹/叁
# plt.rcParams['axes.facecolor'] = '#272822'   # '#3d444c'
# plt.rcParams['axes.labelcolor'] = '#d5d1c7'  # tint grey
# plt.rcParams['axes.axisbelow'] = 'False'   # 'line'
# plt.rcParams['ytick.color'] = '#d5d1c7'  # 'white'
# plt.rcParams['xtick.color'] = '#d5d1c7'  # 'white'
# plt.rcParams['text.color'] = '#d5d1c7'  # 'white'
# plt.rcParams['grid.color'] = '#d5d1c7'  # grid color
# plt.rcParams['grid.linestyle'] = '--'      # solid
# plt.rcParams['grid.linewidth'] = 0.3       # in points
# plt.rcParams['grid.alpha'] = 0.8       # transparency, between 0.0 and 1.0

def get_fig():
    # return plt.figure(figsize=(12, 6), facecolor='#272822') #'#3d444c')
    return plt.figure(figsize=(12, 6))


def sm_sigmoid(x, a):
    return 2.0 / (1.0 + np.exp(-a*x)) - 1.0


# %%
get_fig()
x = np.arange(-1, 1, 0.001)
for a in [0.1, 1, 5, 20, 200, 500]:
    y = sm_sigmoid(x, a)
    plt.plot(x, y, label=f'a={a}')
plt.legend()
plt.grid()

# %%
