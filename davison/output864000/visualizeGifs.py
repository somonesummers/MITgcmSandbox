import imageio
import numpy as np

startStep = 864
endStep = 8640
stepSize = 864

toPlot = ['U','W','T','S']

for q in toPlot:
    filenames = ['']
    images = []
    for i in np.arange(stepSize, endStep+1, stepSize):
        if i == stepSize:
            filenames = ['figs/%s%05i.png' % (q,i)]
        else:
            filenames.append('figs/%s%05i.png' % (q,i))

    print(filenames)

    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave('figs/45%sGif.gif' % q, images)