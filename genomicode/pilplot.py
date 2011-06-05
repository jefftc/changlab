"""

The Python Imaging Library uses a coordinate system with (0, 0) in the
upper left corner.

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
    from PIL import Image
    return Image.new("RGB", (width, height), (255, 255, 255))

def rectangle(image, x, y, width, height, color, outline=None):
    # color is (r, g, b) where r, g, b from 0-255.
    # Is the fill color.  outline is the colorfor the outline.
    from PIL import ImageDraw
    draw = ImageDraw.Draw(image)
    # width and height are exclusive, but PIL is inclusive, so
    # subtract 1 to compensate.
    # outline???
    draw.rectangle((x, y, x+width-1, y+height-1), outline=outline, fill=color)
    del draw
    
def line(image, x, y, width, height, color):
    from PIL import ImageDraw
    draw = ImageDraw.Draw(image)
    draw.line((x, y, x+width-1, y+height-1), fill=color)
    del draw

def polygon(image, coords, color):
    # Coordinates is an array of (x, y).
    from PIL import ImageDraw
    draw = ImageDraw.Draw(image)
    draw.polygon(coords, fill=color)
    del draw
    
def text(image, x, y, text, fontsize, color):
    from PIL import ImageDraw
    font = _load_font(fontsize)
    draw = ImageDraw.Draw(image)
    draw.text((x, y), text, font=font, fill=color)
    del draw

def text90(image, x, y, text, fontsize, color):
    from PIL import Image, ImageDraw
    
    bgcolor = (255, 255, 255)
    font = _load_font(fontsize)
    
    width, height = get_text_size(text, fontsize)
    tim = Image.new("RGB", (width, height), bgcolor)
    draw = ImageDraw.Draw(tim)
    draw.text((0, 0), text, font=font, fill=color)
    del draw

    # BUG: Should paste with a mask, so don't overwrite underlying
    # stuff.
    tim = tim.rotate(90)
    image.paste(tim, (x, y))

def write(image, handle):
    image.save(handle, "png")

def _load_font(size):
    from PIL import ImageFont
    this_module = sys.modules[__name__]
    path = os.path.split(this_module.__file__)[0]
    filename = os.path.join(path, "MS PGothic.ttf")
    assert os.path.exists(filename), "Missing: %s" % filename
    return ImageFont.truetype(filename, size)

def get_text_size(text, fontsize):
    from PIL import Image, ImageDraw, ImageFont
    font = _load_font(fontsize)
    image = Image.new("RGB", (1, 1), (255, 255, 255))
    draw = ImageDraw.Draw(image)
    width, height = draw.textsize(text, font=font)
    del draw
    return width, height

def fit_fontsize_to_height(height):
    # Not needed.  The fontsize is the same as the height!
    #min_size, max_size = 1, 50
    #while min_size <= max_size:
    #    fontsize = (min_size + max_size) / 2
    #    w, h = self._get_text_size("test", fontsize)
    #    if h < height:
    #        min_size = fontsize + 1
    #    elif h > height:
    #        max_size = fontsize - 1
    #    else:
    #        break
    #else:
    #    raise AssertionError, "I could not find the font size."
    #return fontsize
    return height
