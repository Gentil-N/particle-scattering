import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib.widgets  import SpanSelector
from matplotlib.widgets  import Button
import tkinter as tk
from tkinter import filedialog as fd

WINDOW_TITLE = 'Extractor'
IMG = np.array([])
VMIN = 0
VMAX = 0
SELECTION = []

def choose_file():
    global IMG
    global LABEL
    global SELECTION
    filename = fd.askopenfilename()
    if filename != () and filename != "":
        IMG = mpimg.imread(filename)
        LABEL.config(text = filename)
        SELECTION.clear()
        print("Selected file: ", filename)

def select_areas():
    global IMG
    if IMG.size == 0:
        print("Unable to select areas! Choose a file first.")
        return
    fig, ax = plt.subplots()
    plt.imshow(IMG)

    def on_select(vmin, vmax):
        global VMIN
        global VMAX
        VMIN = round(vmin)
        VMAX = round(vmax)
    def add_selection(event):
        global SELECTION
        SELECTION.append((VMIN, VMAX))
        print("Selection added: (", VMIN, ",", VMAX, ")")

    rs = SpanSelector(ax, on_select, 'vertical', props=dict(facecolor='red', alpha=0.35), interactive=True)
    axeselect = plt.axes([0.8, 0.0, 0.2, 0.075])
    bselect = Button(axeselect, 'Add Selection')
    bselect.on_clicked(add_selection)
    plt.show()

def extract():
    global SELECTION
    global IMG
    if len(SELECTION) == 0:
        print("No selection.")
        return
    if IMG.size == 0:
        print("Unable to extract! Choose a file first.")
        return
    for v in SELECTION:
        print("Extraction of selection", v, "...", end='', flush=True)
        data = []
        offset = v[1] - v[0]
        for i in range(len(IMG[0])):
            average = 0
            for j in range(offset):
                average += (int(IMG[v[0] + j][i][0]) + int(IMG[v[0] + j][i][1]) + int(IMG[v[0] + j][i][2])) / 3
            average /= offset
            data.append(average)
        print(" done")


root = tk.Tk(className=WINDOW_TITLE)
root.wm_title(WINDOW_TITLE)
root.geometry('600x150+50+50')
LABEL = tk.Label(root, text='No file selected')
LABEL.pack()
tk.Button(root, text='Choose File', command=choose_file, justify=tk.CENTER).pack()
tk.Button(root, text='Select Areas', command=select_areas, justify=tk.CENTER).pack()
tk.Button(root, text='Extract', command=extract, justify=tk.CENTER).pack()
root.mainloop()

print("haha")