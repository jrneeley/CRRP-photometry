import numpy as np
import lightcurves
from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure


target = 'NGC3201'
optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'
folder = '/Users/jrneeley/CRRP/'+target+'/'
dtype1 = np.dtype([('id', int), ('period', float)])
data = np.loadtxt(folder+target+'-clement.txt', dtype=dtype1, usecols=(0,3))
ind=34
lcv_file = optical_folder+target+'lcvs/'+target+'V'+str(data['id'][ind])+'.lcv'
lcv = 'V'+str(data['id'][ind])
print lcv, data['period'][ind]
U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file, old=0)
best_period = data['period'][ind]
error_threshold = 0.05

xx = np.array(V[2][V[1] < error_threshold], dtype=float)
y = np.array(V[0][V[1] < error_threshold], dtype=float)
er = np.array(V[1][V[1] < error_threshold], dtype=float)

best_freq = 1/best_period
x = np.mod(xx*best_freq, 1)

source = ColumnDataSource(data=dict(x=x, y=y))

# Set up plot
plot = figure(plot_height=400, plot_width=600, title=lcv,
    tools="save",
    x_range=[0, 1], y_range=[y.max()+0.05, y.min()-0.05])

plot.scatter('x', 'y', source=source)

# Set up widgets
freq = Slider(title="frequency", value=best_period, start=best_period-0.0001, end=best_period+0.0001, step=0.00000001)
def update_data(attrname, old, new):

    # Get the current slider values
    f = freq.value
    # Generate the new curve
    x = np.mod(xx/f, 1)
    source.data = dict(x=x, y=y)

for w in [freq]:
    w.on_change('value', update_data)

# Set up layouts and add to document
inputs = widgetbox(freq)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
