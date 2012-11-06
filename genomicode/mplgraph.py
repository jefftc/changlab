#mplgraph.py

"""
Functions:
boxplot
scatter
lineplot
barplot
"""

import numpy
import pylab
from matplotlib import axis
import matplotlib.pyplot as plt


def boxplot(*args,**keywds):
    """Return the pylab figure object.
    data       List of list for the data.
    xlabel      what label to put on x axis.
    ylabel      what label to put on y axis.
    title       what label to put on title.
    box_label   what to put in xticklabel.
    box_color   color for the boxes.<"black','red'>,etc.
    whisker_color    color for the whisker.<"black','red'>,etc.
    tick_size   what font size of the tick label.
    left        the left margin distance
    right       the right margin distance
    top          the top margin distance
    bottom      the bottom margin distance
    """
    assert len(args) == 1, "Specify data"
    data,=args
    data=[numpy.array(i) for i in data]
    xlabel = keywds.get("xlabel",None)
    ylabel = keywds.get("ylabel",None)
    title = keywds.get("title",None)
    box_label = keywds.get("box_label",None)
    tick_size = keywds.get("tick_size",10)
    left = keywds.get("left",0.075)
    right = keywds.get("right",0.95)
    top = keywds.get("top",0.9)
    bottom = keywds.get("bottom",0.25)
    box_color = keywds.get("box_color",'blue')
    whisker_color = keywds.get("whisker_color",'black')
    assert data,'No data provided for the box plot.'
    
    #check the inputs
    if box_label:
        assert len(box_label)==len(data)
    
    fig=pylab.figure()
    bp = pylab.boxplot(data, sym='o', vert=1, whis=1.5)
    if xlabel:
        pylab.xlabel(xlabel)
    if ylabel:
        pylab.ylabel(ylabel)
    if title:
        pylab.title(title)
    if box_label:
        name = box_label
        if len(data)<=12:
            label = box_label
        else:
            label = ['']*len(data)
            index = [int(round(len(data)/12.0*i)) for i in range(12)]
            for i in range(12):
                label[index[i]] = box_label[index[i]]
        ax = pylab.gca()
        ax.set_xticks(range(1,len(data)+1))
        ax.set_xticklabels(tuple(label),rotation='vertical',fontsize=tick_size)
    pylab.setp(bp['boxes'], color=box_color)
    pylab.setp(bp['whiskers'], color=whisker_color)
    pylab.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    return fig

##def scatter(*args,**keywds):
##    """Return the figure object.
##    x1,x2       List of list for the two dimension data.
##    xlabel      what label to put on x axis.
##    ylabel      what label to put on y axis.
##    title       what label to put on title.
##    left        the left margin distance
##    right       the right margin distance
##    top          the top margin distance
##    bottom      the bottom margin distance
##    label       a list of label for each point
##    color       list of color for each point or a single color for all the points
##    legend      list of text for legend
##    """
##    assert len(args) == 2, "Input data should be two dimension"
##    [x1, x2] = args
##    xlabel = keywds.get("xlabel",None)
##    ylabel = keywds.get("ylabel",None)
##    title = keywds.get("title",None)
##    labels = keywds.get('label',None)
##    legend = keywds.get('legend',None)
##    color = keywds.get('color','b')
##    left = keywds.get("left",0.17)
##    right = keywds.get("right",0.88)
##    top = keywds.get("top",0.86)
##    bottom = keywds.get("bottom",0.13)
##    assert [x1,x2],'No data provided for the box plot.'
##    #check the inputs
##    if labels:
##        assert len(labels)==len(x1)
##    fig=pylab.figure()
##    pylab.scatter(x1,x2,marker = 'o',s=50,c=color)
##    if xlabel:
##        pylab.xlabel(xlabel)
##    if ylabel:
##        pylab.ylabel(ylabel)
##    if title:
##        pylab.title(title)
##    pylab.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
##    if labels:
##        txt_height = 0.08*(pylab.ylim()[1] - pylab.ylim()[0])
##        txt_width = 0.04*(pylab.xlim()[1] - pylab.xlim()[0])
##        text_positions = get_text_positions(x1, x2, txt_width, txt_height)
##        for label, x, y,t in zip(labels, x1, x2,text_positions):
##            pylab.annotate(
##        label, 
##        xy = (x, y),xytext = (x-txt_width, t*1.05),
##        textcoords = 'data',size=6,
##        #bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
##        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
##    return fig
def scatter(*args,**keywds):
    """Return the figure object.
    x1,x2       List of list for the two dimension data.
    xlabel      what label to put on x axis.
    ylabel      what label to put on y axis.
    title       what label to put on title.
    left        the left margin distance
    right       the right margin distance
    top          the top margin distance
    bottom      the bottom margin distance
    label       a list of label for each point
    color       list of color for each point or a single color for all the points
    legend      list of text for legend for each point
    """
    assert len(args) == 2, "Input data should be two dimension"
    [x1, x2] = args
    xlabel = keywds.get("xlabel",None)
    ylabel = keywds.get("ylabel",None)
    title = keywds.get("title",None)
    labels = keywds.get('label',None)
    legend = keywds.get('legend',None)
    color = keywds.get('color','b')
    left = keywds.get("left",0.17)
    right = keywds.get("right",0.88)
    top = keywds.get("top",0.86)
    bottom = keywds.get("bottom",0.13)
    assert [x1,x2],'No data provided for the box plot.'
    #check the inputs
    if labels:
        assert len(labels)==len(x1)
    if legend:
        assert len(legend)==len(x1)
    if len(color)>1:
        assert len(color)==len(x1)
    fig=pylab.figure()
    if len(color)>1 and legend:
        old_legend = numpy.array(legend)
        legend = reduce(lambda x, y: x if y in x else x + [y], legend, [])
        color = reduce(lambda x, y: x if y in x else x + [y], color, [])
        for l,t in zip(legend,color):
            pylab.scatter(numpy.array(x1)[old_legend==l],numpy.array(x2)[old_legend==l],label=l,marker = 'o',s=50,c=t)
        pylab.legend()
    elif not legend:  
        pylab.scatter(x1,x2,marker = 'o',s=50,c=color)
    if xlabel:
        pylab.xlabel(xlabel)
    if ylabel:
        pylab.ylabel(ylabel)
    if title:
        pylab.title(title)
    pylab.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    if labels:
        txt_height = 0.08*(pylab.ylim()[1] - pylab.ylim()[0])
        txt_width = 0.04*(pylab.xlim()[1] - pylab.xlim()[0])
        text_positions = get_text_positions(x1, x2, txt_width, txt_height)
        for label, x, y,t in zip(labels, x1, x2,text_positions):
            pylab.annotate(
        label, 
        xy = (x, y),xytext = (x-txt_width, t*1.05),
        textcoords = 'data',size=6,
        #bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    return fig

def get_text_positions(x_data, y_data, txt_width, txt_height):
    a = zip(y_data, x_data)
    text_positions = y_data[:]
    for index, (y, x) in enumerate(a):
        local_text_positions = [i for i in a if i[0] > (y - txt_height) 
                            and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
        if local_text_positions:
            sorted_ltp = sorted(local_text_positions)
            if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
                differ = numpy.diff(sorted_ltp, axis=0)
                a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
                text_positions[index] = sorted_ltp[-1][0] + txt_height
                for k, (j, m) in enumerate(differ):
                    #j is the vertical distance between words
                    if j > txt_height * 2: #if True then room to fit a word in
                        a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
                        text_positions[index] = sorted_ltp[k][0] + txt_height
                        break
    return text_positions

def lineplot(*args,**keywds):
    """Return the pylab figure object.
    line1        List of the (x, y) coordinates of the points.
    line2        ...
    xlabel      what label to put on x axis.
    ylabel      what label to put on y axis.
    title       what label to put on title.
    box_label   what to put in xticklabel.
    tick_size   what font size of the tick label.
    left        the left margin distance
    right       the right margin distance
    top         the top margin distance
    bottom      the bottom margin distance
    legend      list of legend, 1 for each line
    color       List of colors, 1 for each line
    """
    assert args, "No lines given."
    lines = args
    for line in lines:
        for x in line:
            assert len(x) in [2, 3], str(x)
    X=[[x[0] for x in line] for line in lines]
    Y=[[x[1] for x in line] for line in lines]
    xlabel = keywds.get("xlabel",None)
    ylabel = keywds.get("ylabel",None)
    title = keywds.get("title",None)
    box_label = keywds.get("box_label",None)
    tick_size = keywds.get("tick_size",10)
    left = keywds.get("left",0.075)
    right = keywds.get("right",0.95)
    top = keywds.get("top",0.9)
    bottom = keywds.get("bottom",0.25)
    color = keywds.get("color",'b')
    legend = keywds.get("legend",None)
    ylim_min = keywds.get("ylim_min",None)
    assert lines,'No data provided for the line plot.'
    #check the inputs
    if box_label:
        assert len(box_label)==len(lines[0])
    fig=pylab.figure()
    if isinstance(color,list):
        for x,y,c in zip(X,Y,color):
            bp = pylab.plot(x,y,color=c)
    else:
        for x,y in zip(X,Y):
            bp = pylab.plot(x,y)
    if xlabel:
        pylab.xlabel(xlabel)
    if ylabel:
        pylab.ylabel(ylabel)
    if title:
        pylab.title(title)
    if box_label:
        if len(lines[0])<=12:
            label = box_label
        else:
            label = ['']*len(lines[0])
            index = [int(round(len(lines[0])/12.0*i)) for i in range(12)]
            for i in range(12):
                label[index[i]] = box_label[index[i]]
        ax = pylab.gca()
        ax.set_xticks(x)
        ax.set_xticklabels(tuple(label),rotation='vertical',fontsize=tick_size)
    if legend:
        pylab.legend(legend, 'best', shadow=True, fancybox=True)
    if ylim_min is not None:
        b=max(max(Y))
        pylab.ylim((ylim_min,b+1))
    pylab.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    return fig

def barplot(*args,**keywds):
    """Return the pylab figure object.
    mean       List of means for the data.
    std        List of standard variation.
    xlabel      what label to put on x axis.
    ylabel      what label to put on y axis.
    title       what label to put on title.
    box_label   what to put in xticklabel.
    tick_size   what font size of the tick label.
    left        the left margin distance.
    right       the right margin distance.
    top          the top margin distance.
    bottom      the bottom margin distance.
    xtick_rotation rotation the box_label in vertical.
    ylim       (min,max) of the limit of y_axis.
    ytick_pos   List of position to put yticks.
    yticks      List of ticks to put in y_axis.
    """
    assert len(args) >=1, "Specify data"
    if len(args)==1:
        mean, = args
        std = [0]*len(mean)
    elif len(args)==2:
        mean,std = args
    xlabel = keywds.get("xlabel",None)
    ylabel = keywds.get("ylabel",None)
    title = keywds.get("title",None)
    box_label = keywds.get("box_label",None)
    tick_size = keywds.get("tick_size",10)
    left = keywds.get("left",0.085)
    right = keywds.get("right",0.95)
    top = keywds.get("top",0.9)
    bottom = keywds.get("bottom",0.35)
    xtick_rotation  = keywds.get('xtick_rotation',None)
    assert mean,'No data provided for the bar plot.'
    #check the inputs
    if box_label:
        assert len(box_label)==len(mean)
    ind = numpy.arange(len(mean))
    width = 0.35 
    fig=pylab.figure()
    bp = pylab.bar(ind,mean,width,yerr=std,color='y')
    
    if xlabel:
        pylab.xlabel(xlabel)
    if ylabel:
        pylab.ylabel(ylabel)
    if title:
        pylab.title(title)
    if box_label:
        if len(mean)<=12:
            label = box_label
        else:
            label = ['']*len(mean)
            index = [int(round(len(mean)/12.0*i)) for i in range(12)]
            for i in range(12):
                label[index[i]] = box_label[index[i]]
        
        pylab.xticks(ind+width/2.,label,rotation=xtick_rotation)
    if keywds.get('ylim'):
        pylab.ylim(keywds.get('ylim'))
    if keywds.get('ytick_pos') and keywds.get('yticks'):
        pylab.yticks(keywds.get('ytick_pos'),keywds.get('yticks'))
        
    pylab.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    return fig
