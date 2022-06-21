from bokeh.core.property.instance import Instance
from bokeh.io import save
from bokeh.layouts import column
from bokeh.model import Model
from bokeh.models import CustomJS, Slider, Callback
from bokeh.plotting import ColumnDataSource, figure, show

source_k1 = ColumnDataSource(data=dict(t=[], s=[]))
source_k2 = ColumnDataSource(data=dict(t=[], i=[]))
source_k3 = ColumnDataSource(data=dict(t=[], r=[]))
# print(source.data)
plot = figure(plot_width=400, plot_height=400)
plot.line('t', 's', source=source_k1, line_width=3, line_alpha=0.6)
plot.line('t', 'i', source=source_k2, line_width=3, line_alpha=0.6, color='orange')
plot.line('t', 'r', source=source_k3, line_width=3, line_alpha=0.6, color='green')

callback_k1 = CustomJS(args=dict(source=source_k1), code="""\
    const N = 200;
    let s = 489599 / 489609;
    let i = 10 / 489609;
    let r = 0 / 489609;
    const bet = cb_obj.value;
    const gam = 0.189;
    const tlst = source.data.t = [];
    const slst = source.data.s = [];
    for (let t = 0; t < N; ++t) {
        s -= bet * s * i;
        i += bet * s * i - gam * i;
        r += gam * i;
        tlst.push(t);
        slst.push(s);
    }
    source.change.emit();
""")

callback_k2 = CustomJS(args=dict(source=source_k2), code="""\
    const N = 200;
    let s = 489599 / 489609;
    let i = 10 / 489609;
    let r = 0 / 489609;
    const bet = cb_obj.value;
    const gam = 0.189;
    const tlst = source.data.t = [];
    const ilst = source.data.i = [];
    for (let t = 0; t < N; ++t) {
        i += bet * s * i - gam * i;
        tlst.push(t);
        ilst.push(i);
    }
    source.change.emit();
""")

callback_k3 = CustomJS(args=dict(source=source_k3), code="""\
    const N = 200;
    let s = 489599 / 489609;
    let i = 10 / 489609;
    let r = 0 / 489609;
    const bet = cb_obj.value;
    const gam = 0.189;
    const tlst = source.data.t = [];
    const rlst = source.data.r = [];
    for (let t = 0; t < N; ++t) {
        r += gam * i;
        tlst.push(t);
        rlst.push(r);
    }
    source.change.emit();
""")

slider_k1 = Slider(start=0.1, end=2, value=1, step=.1, title="k1 ")
slider_k2 = Slider(start=0.1, end=2, value=1, step=.1, title="k2 ")
slider_k3 = Slider(start=0.1, end=2, value=1, step=.1, title="k3 ")
# # slider_k4 = Slider(start=0.1, end=4, value=1, step=.1, title="k4 ")

slider_k1.js_on_change('value', callback_k1)
slider_k2.js_on_change('value', callback_k2)
slider_k3.js_on_change('value', callback_k3)
# # slider_k4.js_on_change('value', callback)


layout = column(slider_k1, slider_k2, slider_k3, plot)
show(layout)


