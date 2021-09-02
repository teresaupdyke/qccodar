#!/usr/bin/env python
#
# Last modified: Time-stamp: <2015-10-12 16:01:53 haines>
""" new classes of slider 
"""

import sys
import matplotlib.widgets
import numpy
from functools import reduce

class DiscreteSlider(matplotlib.widgets.Slider):
    """A matplotlib slider widget with discrete steps."""
    def __init__(self, *args, **kwargs):
        """Identical to Slider.__init__, except for the "increment" kwarg.
        "increment" specifies the step size that the slider will be discritized
        to."""
        self.inc = kwargs.pop('increment', 0.5)
        matplotlib.widgets.Slider.__init__(self, *args, **kwargs)

    def set_val(self, val):
        discrete_val = int(val / self.inc) * self.inc
        # We can't just call Slider.set_val(self, discrete_val), because this 
        # will prevent the slider from updating properly (it will get stuck at
        # the first step and not "slide"). Instead, we'll keep track of the
        # the continuous value as self.val and pass in the discrete value to
        # everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: 
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: 
            return
        for cid, func in self.observers.items():
            func(discrete_val)

class IndexedSlider(matplotlib.widgets.Slider):
    """A matplotlib slider widget with discrete steps that index a list of values."""
    def __init__(self, *args, **kwargs):
        """Identical to Slider.__init__, except for the "seqvals" kwarg.
        "seqvals" specifies the sequence of values to be indexed by
        the slider.  The slider will be discritized so that self.val
        is an index to specific value in "seqvals".
        """
        seqvals = kwargs.pop('seqvals', list(range(1,11)))
        if type(seqvals) is type(numpy.array([])):
            self.seqvals = seqvals.tolist()
        elif type(seqvals) is type([]):
            self.seqvals = seqvals
        else:
            print('Input param "seqvals" must be a list or numpy.array.')
        self.step = 1

        # fix args to pass index values for valmin and valmax
        args = list(args)
        if len(args) <= 2:
            args.append(0)
            args.append(len(self.seqvals))
        elif len(args) == 4:
            args[2] = 0
            args[3] = len(self.seqvals)

        # fix kwargs['valinit'] to set initial index 
        valinit = kwargs.pop('valinit', None)
        if not valinit:
            self.valinit = 0
        elif self.seqvals.count(valinit):
            # set valinit as index of seqvalinit in seqvals
            self.valinit = self.seqvals.index(valinit)
        else:
            # look around for closest value by finding index of min(absdiff) 
            min = (lambda x, y: x if x < y else y)
            absdiff = [abs(x-valinit) for x in seqvals]
            idx = absdiff.index(reduce(min, absdiff))
            seqvalinit = self.seqvals[idx]
            # set valinit as index of seqvalinit in seqvals
            self.valinit = self.seqvals.index(seqvalinit)
        # put the value back in kwargs to inherit for init
        kwargs['valinit'] = self.valinit

        matplotlib.widgets.Slider.__init__(self, *args, **kwargs)
        self.valtext.set_text(self.valfmt % seqvals[self.valinit])

    def set_val(self, val):
        discrete_val = int(val / self.step) * self.step
        # keep track of the the continuous value as self.val and pass
        # in the discrete value to everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % self.seqvals[discrete_val])
        if self.drawon: 
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: 
            return
        for cid, func in self.observers.items():
            func(discrete_val)


def _test():
    import matplotlib.pyplot as plt
    import matplotlib.pylab as pylab
    
    def discrete_slider_update(val):
        dot.set_xdata([val])
        dot.set_ydata([val])
        fig.canvas.draw()
        
    def indexed_slider_update(val):
        dot.set_xdata([val])
        dot.set_ydata([isbear.seqvals[val]])
        fig.canvas.draw()
        
    fig, ax = plt.subplots()
        
    # seqvals = range(10,40,1)
    # seqvals = range(355,359,1)
    # seqvals.extend(range(0,147,1))
    seqvals = list(range(50,0,-1))
    xseq = list(range(0, len(seqvals), 1))
    
    # trying the discrete slider first
    dsax = fig.add_axes([0.2, 0.2, 0.6, 0.03])
    dsbear = DiscreteSlider(dsax, 'Discrete', 0, 10, increment=1)
    dsbear.on_changed(discrete_slider_update)

    # now try the indexed slider
    isax = fig.add_axes([0.2, 0.1, 0.6, 0.03])
    # isbear = IndexedSlider(isax, 'Indexed')
    isbear = IndexedSlider(isax, 'Indexed', 50, 0, valinit=25.5, seqvals=seqvals)
    isbear.on_changed(indexed_slider_update)
    
    ax.plot(xseq, seqvals, 'ro')
    dot, = ax.plot(isbear.valinit, seqvals[isbear.valinit], 'bo', markersize=18)
