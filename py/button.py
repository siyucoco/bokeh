from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import Figure, output_file, show
from bokeh.models.widgets import Button

x = [x*0.05 for x in range(0, 200)]
y = x

source = ColumnDataSource(data=dict(x=x, y=y))

plot = Figure(plot_width = 400, plot_height = 400)
plot.line('x', 'y', source = source, line_width=3, line_alpha = 0.6)

callback = CustomJS(args=dict(source= source), code="""
    var data = source.data;
    var x = data['x']
    var y = data['y']
    for(var i =0; i<x.length; i++){
        y[i]=Math.pow(x[i], 4)
    }
    source.change.emit();
    """)

btn = Button(label="click here", callbacks = callback, name = "1")
layout = column(btn, plot)
show(layout)