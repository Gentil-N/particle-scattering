import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib.widgets  import SpanSelector
from matplotlib.widgets  import Button

img = mpimg.imread('./res/test.webp')

#xdata = np.linspace(0,9*np.pi, num=301)
#ydata = np.sin(xdata)

fig, ax = plt.subplots()
plt.imshow(img)
#line, = ax.plot(xdata, ydata)

def onselect(vmin, vmax):
    print(vmin, vmax)

rs = SpanSelector(ax, onselect, 'vertical', props=dict(facecolor='grey', alpha=0.35), interactive=True)
axcut = plt.axes([0.9, 0.0, 0.1, 0.075])
bcut = Button(axcut, 'Extract', color='grey', hovercolor='cyan')
bcut.on_clicked(exctract)

plt.show()
print("haha")