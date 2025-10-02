import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
from ZZ_basicUtility import pathFixer
from reader import ParquetChunkReader
import time

def init(title, meta):
    """Initializes:\n
    Qt: `app`, `win`, `plot`\n
    data: `circs_outer`, `circs_inner`, `lines_item`"""
    global app, win, plot, circs_outer, circs_inner, lines_item, text_items, reader
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication([])
    win = pg.GraphicsLayoutWidget(show=True, title=title)
    plot = win.addPlot()
    plot.setAspectLocked(True)
    plot.setXRange(-meta["A1"], meta["A1"], padding=0)
    plot.setYRange(-meta["A2"], meta["A2"], padding=0)
    plot.setMouseEnabled(False, False)

    # Initialisierung
    circs_outer = []
    circs_inner = []
    lines_item = pg.PlotDataItem(pen=pg.mkPen('w', width=2), connect='pairs')
    text_items = []

    text = pg.TextItem("halo", anchor=(0.5, 0.5))
    text.setPos(-0.8*reader.metaData["A1"], 0.8*reader.metaData["A2"])
    plot.addItem(text)
    text_items.append(text)

    for i in range(meta["NN"]):
        co = QtWidgets.QGraphicsEllipseItem(-r_outer, -r_outer, 2*r_outer, 2*r_outer)
        co.setBrush(pg.mkBrush("w"))
        co.setPen(pg.mkPen(None))
        circs_outer.append(co)

        ci = QtWidgets.QGraphicsEllipseItem(-r_inner, -r_inner, 2*r_inner, 2*r_inner)
        ci.setPen(pg.mkPen(None))
        circs_inner.append(ci)

        plot.addItem(co)
        plot.addItem(ci)

    plot.addItem(lines_item)

def draw_frame(reader):
    global r_inner, r_outer, circs_outer, circs_inner, lines_item, scientists
    try:
        frame = reader.__next__()
    except RuntimeError as e:
        timer.disconnect()
        raise e
    except:
        timer.disconnect()
        return


    lines_x, lines_y = [], []
    ws, xs, ys, signs = frame
    for i in range(len(xs)):
        w, x, y, sign = ws[i], xs[i], ys[i], signs[i]
        circs_outer[i].setPos(x, y)
        circs_inner[i].setPos(x, y)

        color = 'b' if sign > 0 else 'r'
        circs_inner[i].setBrush(pg.mkBrush(color))

        x_end = x + r_inner * np.cos(w)
        y_end = y + r_inner * np.sin(w)
        lines_x += [x, x_end]
        lines_y += [y, y_end]

    text_items[0].setText(str(round(reader.globalTime * metaData["dt_log"], 1)) + " s")


    lines_item.setData(lines_x, lines_y)

path = r".\outData\staticMem32"
path = pathFixer(path) + "spinner_T1800_LC21_density1.500000_s1_0.parquet"

reader = ParquetChunkReader(path, skipFrames=2, preFormatFunction=lambda x: [x[i * len(x)//4 : (i+1) * len(x)//4] for i in range(4)])
metaData = reader.metaData
layout = np.zeros(shape=reader.shape)
layout[0, :reader.metaData["NN"]] = 1
layout[1, :reader.metaData["NN"]] = 1
layout[2, :reader.metaData["NN"]] = 1
layout[metaData["gapPosition1"], :reader.metaData["NN"]] = 1
reader.mask = np.array(layout, dtype=bool)

r_outer = metaData["r_out"]
r_inner = metaData["r_inner"]
dt = metaData["dt_log"]

init("Spinner Animation", metaData)
timer = pg.QtCore.QTimer()
timer.timeout.connect(lambda: draw_frame(reader))
timer.start(1 + 4) # ms

app.exec_()