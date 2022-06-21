from bokeh.core.property.instance import Instance
from bokeh.io import save
from bokeh.layouts import column
from bokeh.model import Model
from bokeh.models import CustomJS, Slider, Callback
from bokeh.plotting import ColumnDataSource, figure, show

source = ColumnDataSource(data=dict(t=[], s=[], i=[], r=[]))

plot = figure(plot_width=400, plot_height=400)
plot.line('t', 's', source=source, line_width=3, line_alpha=0.6)
# plot.line('t', 'i', source=source, line_width=3, line_alpha=0.6, color='orange')
# plot.line('t', 'r', source=source, line_width=3, line_alpha=0.6, color='green')
slider1 = Slider(start=0.1, end=4, value=1, step=.1, title="Alpha ")
beta_slider = Slider(start=0.1, end=4, value=1, step=.1, title="Beta ")
callback = CustomJS(args=dict(source=source, b = beta_slider, a = slider1), code="""\
    const N = 200;
    let s = 489599 / 489609;
    let i = 10 / 489609;
    let r = 0 / 489609;
    const bet = b.value;
    const gam = a.value;
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

beta_slider.js_on_change('value', callback)
slider1.js_on_change('value', callback)


layout = column(beta_slider, slider1, plot)
show(layout)


# class IdleDocObserver(Model):
#     """Work around https://github.com/bokeh/bokeh/issues/4272."""
#     on_idle = Instance(Callback)

#     # language=TypeScript
#     __implementation__ = """\
#         import {View} from "core/view"
#         import {Model} from "model"
#         import * as p from "core/properties"

#         export class IdleDocObserverView extends View {}

#         export namespace IdleDocObserver {
#             export type Attrs = p.AttrsOf<Props>
#             export type Props = Model.Props & {on_idle: p.Property<any>}
#         }

#         export interface IdleDocObserver extends IdleDocObserver.Attrs {}

#         export class IdleDocObserver extends Model {
#             static init_IdleDocObserver(): void {
#                 this.prototype.default_view = IdleDocObserverView
#                 this.define<IdleDocObserver.Props>({on_idle: [p.Any]})
#             }

#             _doc_attached(): void {
#                 super._doc_attached()
#                 const execute = () => this.on_idle!.execute(this)
#                 if (this.document!.is_idle)
#                     execute();
#                 else
#                     this.document!.idle.connect(execute);
#             }
#         }
#     """


# idle_doc_observer = IdleDocObserver(on_idle=CustomJS(args=dict(callback=callback, slider=slider),
#                                                      code="callback.execute(slider);"))

# save([idle_doc_observer, layout])