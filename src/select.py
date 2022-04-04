import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets  import SpanSelector
from matplotlib.widgets  import Button
import tkinter as tk
from tkinter import filedialog as fd

class AscData:
    def __init__(self):
        self.data = np.zeros((1))
        self.min_w = 0.0
        self.max_w = 0.0
        self.max_val = 0.0
    def to_image(self):
        img = np.zeros((256, 1024, 3))
        for i in range(256):
            for j in range(1024):
                img[i][j][0] = self.data[i][j] / self.max_val
                img[i][j][1] = self.data[i][j] / self.max_val
                img[i][j][2] = self.data[i][j] / self.max_val
        return img
    def extract_curve(self, vmin, vmax, bck, ref):
        x = np.linspace(self.min_w, self.max_w, 1024)
        y = []
        offset = vmax - vmin
        for i in range(1024):
            average = 0
            for j in range(offset):
                average += (self.data[vmin + j][i] - bck.data[vmin + j][i]) / (ref.data[vmin + j][i] - bck.data[vmin + j][i])
            average /= offset
            y.append(average)
        return x, y
    def read_file(self, filename):
        asc_file = open(filename, 'r')
        lines = asc_file.readlines()
        self.data = np.zeros((256, 1024))
        self.min_w = float(lines[0].split()[0].replace(',', '.'))
        self.max_w = float(lines[-1].split()[0].replace(',', '.'))
        i = 0
        for line in lines:
            line_data = line.split()
            for j in range(1, len(line_data) - 1):
                val = float(line_data[j])
                self.max_val = max(val, self.max_val)
                self.data[j][i] = val
            i += 1
        asc_file.close()

WINDOW_TITLE = 'Extractor'
ASC_DATA = AscData()
ASC_BCK = AscData()
ASC_REF = AscData()
IMG = np.array([])
VMIN = 0
VMAX = 0
SELECTION = []

def update_num_of_select():
    global SELECTION
    global LABEL_2
    LABEL_2.config(text = "Num Of Selections: "+ str(len(SELECTION)))

def choose_file():
    global ASC_DATA
    global IMG
    global LABEL
    global SELECTION
    filename = fd.askopenfilename()
    if filename != () and filename != "":
        print("Reading...", end='', flush=True)
        ASC_DATA.read_file(filename)
        IMG = ASC_DATA.to_image()
        LABEL.config(text = filename)
        SELECTION.clear()
        update_num_of_select()
        print("done")
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
        global VMIN
        global VMAX
        global LABEL_2
        if VMIN == VMAX:
            print("Empty Selection!")
            return
        SELECTION.append((VMIN, VMAX))
        print("Selection added: (", VMIN, ",", VMAX, ")")
        update_num_of_select()
    
    rs = SpanSelector(ax, on_select, 'vertical', props=dict(facecolor='yellow', alpha=0.1), interactive=True)
    axeselect = plt.axes([0.8, 0.0, 0.2, 0.075])
    bselect = Button(axeselect, 'Add Selection')
    bselect.on_clicked(add_selection)
    plt.show()
    plt.delaxes(ax)

def clear_selections():
    global SELECTION
    global LABEL_2
    SELECTION.clear()
    print("Selections cleared.")
    update_num_of_select()

def extract():
    global SELECTION
    global ASC_DATA
    global ASC_BCK
    global ASC_REF
    if len(SELECTION) == 0:
        print("No selection.")
        return
    figs = []
    count = 1
    for v in SELECTION:
        print("Extraction of selection", v, "...", end='', flush=True)
        figs.append(plt.figure(count))
        x, y = ASC_DATA.extract_curve(v[0], v[1], ASC_BCK, ASC_REF)
        plt.title("(" + str(v[0]) + ", " + str(v[1]) + ")")
        plt.plot(x, y)
        print(" done")
        count += 1
    print("Showing...", end='', flush=True)
    plt.show()
    print("done")

def choose_bck_ref():
    global ASC_BCK
    global ASC_REF
    filename_bg = fd.askopenfilename()
    if filename_bg != () and filename_bg != "":
        print("Reading background...", end='', flush=True)
        ASC_BCK.read_file(filename_bg)
        print("done")
        print("Selected background: ", filename_bg)
    filename_ref = fd.askopenfilename()
    if filename_ref != () and filename_ref != "":
        print("Reading reference...", end='', flush=True)
        ASC_REF.read_file(filename_ref)
        print("done")
        print("Selected reference: ", filename_ref)

root = tk.Tk(className=WINDOW_TITLE)
root.wm_title(WINDOW_TITLE)
root.geometry('400x210+50+50')
LABEL = tk.Label(root, text='No file selected')
LABEL.pack()
LABEL_2 = tk.Label(root, text='No Selection')
LABEL_2.pack()
update_num_of_select()
tk.Button(root, text='Choose File', command=choose_file, justify=tk.CENTER).pack()
tk.Button(root, text='Select Areas', command=select_areas, justify=tk.CENTER).pack()
tk.Button(root, text='Clear Selections', command=clear_selections, justify=tk.CENTER).pack()
tk.Button(root, text='Extract', command=extract, justify=tk.CENTER).pack()
tk.Button(root, text='Choose Bck-Ref', command=choose_bck_ref, justify=tk.CENTER).pack()
root.mainloop()

print("Program well terminated!")