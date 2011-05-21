"""

Functions:
image
rectangle
line
text
write

get_text_size
fit_fontsize_to_height

"""
import sys, os

def image(width, height):
    # XXX add camera, etc.
    raise NotImplementedError

def rectangle(image, x, y, width, height, color, outline=None):
    # XXX
    raise NotImplementedError
    
def line(image, x, y, width, height, color):
    raise NotImplementedError

def polygon(image, coords, color):
    raise NotImplementedError
    
def text(image, x, y, text, fontsize, color):
    raise NotImplementedError

def text90(image, x, y, text, fontsize, color):
    raise NotImplementedError

def write(image, handle):
    raise NotImplementedError

def get_text_size(text, fontsize):
    raise NotImplementedError

def fit_fontsize_to_height(height):
    raise NotImplementedError
