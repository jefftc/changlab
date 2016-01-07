"""

Colors are specified as (R, G, B) tuples with values from from 0.0 to 1.0.

red_shade
green_shade
white_shade
rg_array_colors
by_array_colors
red_green_soft
red_blue_soft
rgb_colors          From red, green, blue colorwheel.
matlab_colors
bild_colors
broad_colors
genespring_colors
yahoo_weather_colors

brewer_prgn_div
brewer_rdbu_div
brewer_rdylbu_div
brewer_rdylgn_div
brewer_spectral_div



choose_contrasting_color
choose_contrasting_bw

rgb2hex
rgb2hsl
hsl2rgb

"""


def _matrix2color(matrix, pos):
    # Make sure matrix is sorted from low to high.
    matrix = sorted(matrix)
    
    # pos is [0, 1]
    # Return (R, G, B) where R, G, and B are percentages from 0 to 1.
    breaks = [x[0] for x in matrix]
    
    # Set i1 to the index of the biggest breaks that is <= pos.  i1
    # can range from [0, len(breaks))
    i1 = len([x for x in breaks if pos >= x]) - 1
    assert i1 >= 0
    x = matrix[i1][1:]
    if i1 < len(matrix)-1:
        # pos is [breaks[i1] , breaks[i1+1]).  Interpolate with the
        # next index.
        i2 = i1 + 1

        delta = float(pos-breaks[i1]) / (breaks[i2]-breaks[i1])

        r1, g1, b1 = matrix[i1][1:]
        r2, g2, b2 = matrix[i2][1:]
        r = r1 + delta*(r2-r1)
        g = g1 + delta*(g2-g1)
        b = b1 + delta*(b2-b1)
        x = r, g, b
    return [x/255.0 for x in x]

def _colors(color_matrix, n):
    assert n > 0, "Need at least 1 color."
    if n == 1:
        x = [_matrix2color(color_matrix, 0.5)]
    else:
        x = [_matrix2color(color_matrix, float(i)/(n-1)) for i in range(n)]
    return x

def red_shade(n):
    color_matrix = [
        (0.0,   0, 0, 0),
        (1.0, 255, 0, 0),
        ]
    x = _colors(color_matrix, n)
    return x

def green_shade(n):
    color_matrix = [
        (0.0,   0, 0, 0),
        (1.0, 0, 255, 0),
        ]
    x = _colors(color_matrix, n)
    return x

def white_shade(n):
    color_matrix = [
        (0.0,   0, 0, 0),
        (1.0, 255, 255, 255),
        ]
    x = _colors(color_matrix, n)
    return x

def rg_array_colors(n):
    color_matrix = [
        (0.0,   0, 255, 0),
        (0.5,   0,   0, 0),
        (1.0, 255,   0, 0),
        ]
    x = _colors(color_matrix, n)
    return x

def by_array_colors(n):
    color_matrix = [
        (0.0,   0,   0, 255),
        (0.5,   0,   0,   0),
        (1.0, 255, 255,   0),
        ]
    x = _colors(color_matrix, n)
    return x

def red_green_soft(n):
    color_matrix = [
        (0.00,   52, 112,  11),
        (0.25,   95, 173,  93),
        (0.50,  245, 245, 245),
        #(0.75,  229,  89,  95),
        #(1.00,  249,  39,  39),
        (0.75,  229,  115,  105),
        (1.00,  128,  52,  29),
        ]
    x = _colors(color_matrix, n)
    return x

def red_blue_soft(n):
    color_matrix = [
        (1.00,   86,  15,  24),
        (0.75,  213,  96,  76),
        (0.50,  245, 245, 245),
        (0.25,   79, 148, 194),
        (0.00,    8,  34,  79),
        ]
    x = _colors(color_matrix, n)
    return x

def rgb_colors(n):
    color_matrix = [
        (0.00,   0,   0, 255),
        (0.25,   0, 255, 255),
        (0.50,   0, 255,   0),
        (0.75, 255, 255,   0),
        (1.00, 255,   0,   0),
        ]
    x = _colors(color_matrix, n)
    return x

def matlab_colors(n):
    color_matrix = [
        (0.000,   0,   0, 143),
        (0.125,   0,   0, 255),
        (0.250,   0, 127, 255),
        (0.375,   0, 255, 255),
        (0.500, 127, 255, 127),
        (0.625, 255, 255,   0),
        (0.750, 255, 127,   0),
        (0.875, 255,   0,   0),
        (1.000, 127,   0,   0),
        ]
    x = _colors(color_matrix, n)
    return x

def bild_colors(n):
    color_matrix = [
        (0.000,  49,  50, 114),
        (0.050,  61,  69, 137),
        (0.100,  62,  84, 154),
        (0.150,  67,  89, 160),
        (0.200,  85, 108, 176),
        (0.250, 115, 145, 201),
        (0.300, 160, 205, 240),
        (0.350, 180, 220, 243),
        (0.400, 169, 216, 211),
        (0.450, 160, 208, 164),
        (0.500, 179, 213, 112),
        (0.550, 203, 220,  61),
        (0.600, 232, 231,  61),
        (0.650, 255, 234,  47),
        (0.700, 250, 180,  50),
        (0.750, 243, 136,  54),
        (0.800, 231,  80,  61),
        (0.850, 218,  54,  55),
        (0.900, 204,  55,  59),
        (0.950, 160,  52,  52),
        (1.000, 114,  39,  44),
        ]
    x = _colors(color_matrix, n)
    return x

def broad_colors(n):
    # Default color scheme in GenePattern HeatMapImage module.
    color_matrix = [
        (0.000,  69,   0, 173),
        (0.091,  39,   0, 209),
        (0.182, 107,  88, 239),
        (0.273, 136, 136, 255),
        (0.364, 199, 193, 255),
        (0.455, 213, 213, 255),
        (0.545, 255, 192, 229),
        (0.636, 255, 137, 137),
        (0.727, 255, 112, 128),
        (0.818, 255,  90,  90),
        (0.909, 239,  64,  64),
        (1.000, 214,  12,   0),
        ]
    x = _colors(color_matrix, n)
    return x


def genespring_colors(n):
    color_matrix = [
        (0.0,   0,   0, 255),
        (0.5, 255, 255,   0),
        (1.0, 255,   0,   0),
        ]
    x = _colors(color_matrix, n)
    return x

def yahoo_weather_colors(n):
    color_matrix = [
        (0.0, 255, 255, 255),
        (0.1, 204, 255, 255),
        (0.2, 153, 255, 255),
        (0.3, 102, 204, 255),
        (0.4,  84, 169, 255),
        (0.5, 204, 255, 103),
        (0.6, 255, 255, 103),
        (0.7, 255, 204, 102),
        (0.8, 255, 153, 102),
        (0.9, 204, 102, 102),
        (1.0, 209,  73,  73),
        ]
    x = _colors(color_matrix, n)
    return x


def brewer_prgn_div(n):
    color_matrix = [
        (1.00, 61, 1, 80),
        (0.90, 119, 44, 135),
        (0.80, 154, 112, 171),
        (0.70, 195, 166, 208),
        (0.60, 232, 212, 233),
        (0.50, 248, 248, 248),
        (0.40, 217, 241, 211),
        (0.30, 169, 221, 162),
        (0.20, 79, 174, 106),
        (0.10, 26, 128, 63),
        (0.00, 0, 69, 30),
        ]
    x = _colors(color_matrix, n)
    return x
    

def brewer_rdbu_div(n):
    color_matrix = [
        (1.00, 103, 0, 28),
        (0.90, 177, 20, 42),
        (0.80, 215, 96, 78),
        (0.70, 245, 165, 130),
        (0.60, 253, 220, 200),
        (0.50, 248, 248, 248),
        (0.40, 210, 230, 240),
        (0.30, 146, 198, 223),
        (0.20, 67, 147, 197),
        (0.10, 33, 102, 173),
        (0.00, 6, 46, 95),
        ]
    x = _colors(color_matrix, n)
    return x
    

def brewer_rdylbu_div(n):
    color_matrix = [
        (1.00, 154, 1, 64),
        (0.90, 216, 47, 39),
        (0.80, 245, 107, 70),
        (0.70, 254, 175, 100),
        (0.60, 254, 226, 144),
        (0.50, 254, 255, 182),
        (0.40, 227, 245, 249),
        (0.30, 171, 217, 235),
        (0.20, 116, 174, 209),
        (0.10, 69, 118, 181),
        (0.00, 44, 52, 147),
        ]
    x = _colors(color_matrix, n)
    return x


def brewer_rdylgn_div(n):
    color_matrix = [
        (1.00, 157, 1, 55),
        (0.90, 216, 47, 39),
        (0.80, 245, 107, 70),
        (0.70, 254, 175, 100),
        (0.60, 254, 225, 139),
        (0.50, 254, 255, 192),
        (0.40, 218, 240, 140),
        (0.30, 164, 217, 106),
        (0.20, 103, 190, 99),
        (0.10, 26, 152, 76),
        (0.00, 10, 111, 54),
        ]
    x = _colors(color_matrix, n)
    return x


def brewer_spectral_div(n):
    color_matrix = [
        (1.00, 150, 0, 69),
        (0.90, 213, 62, 80),
        (0.80, 245, 107, 70),
        (0.70, 254, 175, 100),
        (0.60, 254, 225, 139),
        (0.50, 254, 255, 192),
        (0.40, 231, 246, 152),
        (0.30, 169, 221, 163),
        (0.20, 102, 195, 165),
        (0.10, 61, 143, 194),
        (0.00, 102, 81, 164),
        ]
    x = _colors(color_matrix, n)
    return x


def hex2rgb(h):
    assert h.startswith("0x") or h.startswith("#")
    if h.startswith("0x"):
        h = h[2:]
    if h.startswith("#"):
        h = h[1:]
    assert len(h) == 6
    r, g, b = h[:2], h[2:4], h[4:]
    r, g, b = int(r, 16), int(g, 16), int(b, 16)
    r, g, b = r/255.0, g/255.0, b/255.0
    return r, g, b


def rgb2hex(c):
    # red, green, blue values from 0-1.
    r, g, b = c
    NUMERIC_TYPE = [type(0.0), type(0)]
    assert type(r) in NUMERIC_TYPE
    assert type(g) in NUMERIC_TYPE
    assert type(b) in NUMERIC_TYPE
    r, g, b = int(r*255), int(g*255), int(b*255)
    assert r >= 0 and r < 256
    assert g >= 0 and g < 256
    assert b >= 0 and b < 256
    r = hex(r)[2:].upper()
    g = hex(g)[2:].upper()
    b = hex(b)[2:].upper()
    if len(r) == 1:
        r = "0%s" % r
    if len(g) == 1:
        g = "0%s" % g
    if len(b) == 1:
        b = "0%s" % b
    x = "0x%s%s%s" % (r, g, b)
    assert len(x) == 8
    return x

def choose_contrasting_color(col):
    r, g, b = col

    assert r >= 0 and r <= 1
    assert g >= 0 and g <= 1
    assert b >= 0 and b <= 1

    h, s, l = rgb2hsl(col)
    l = l + 0.50
    if l > 1.0:
        l -= 1
    col = hsl2rgb((h, s, l))
    return col

def choose_contrasting_bw(col):
    r, g, b = col
    assert r >= 0 and r <= 1
    assert g >= 0 and g <= 1
    assert b >= 0 and b <= 1

    Y = 0.2126*(r**2.2) + 0.7151*(g**2.2) + 0.0721*(b**2.2)
    contrast = 0, 0, 0       # black
    if Y <= 0.18:
        contrast = 1, 1, 1   # white
    return contrast

def rgb2hsl(col):
    r, g, b = col
    assert r >= 0 and r <= 1
    assert g >= 0 and g <= 1
    assert b >= 0 and b <= 1
    r, g, b = float(r), float(g), float(b)

    maxcolor = max(r, g, b)
    mincolor = min(r, g, b)
    L = (maxcolor+mincolor)/2.0
    if maxcolor == mincolor:
        S = H = 0
        return H, S, L
        
    if L < 0.5:
        S = (maxcolor-mincolor)/(maxcolor+mincolor)
    else:
        S = (maxcolor-mincolor)/(2.0-maxcolor-mincolor)

    if r == maxcolor:
        H = (g-b)/(maxcolor-mincolor)
    elif g == maxcolor:
        H = 2.0+(b-r)/(maxcolor-mincolor)
    else:
        H = 4.0+(r-g)/(maxcolor-mincolor)

    # Scale H to 0-360.
    H = H*60.0
    if H < 0:
        H += 360
    return H, S, L

def _hsl2rgb_h(H, temp1, temp2, temp3):
    # What is H for?
    if temp3 < 0:
        temp3 += 1
    elif temp3 > 1:
        temp3 -= 1

    if 6.0*temp3 < 1:
        c = temp1 + (temp2-temp1)*6.0*temp3
    elif 2.0*temp3 < 1:
        c = temp2
    elif 3.0*temp3 < 2:
        c = temp1+(temp2-temp1)*((2.0/3.0)-temp3)*6.0
    else:
        c = temp1
    return c

def hsl2rgb(col):
    H, S, L = col

    if S == 0:
        return L, L, L

    if L < 0.5:
        temp2 = L*(1.0+S)
    else:
        temp2 = L+S - L*S
    temp1 = 2.0*L - temp2

    H = H/360.0

    R = _hsl2rgb_h(H, temp1, temp2, H+1.0/3.0)
    G = _hsl2rgb_h(H, temp1, temp2, H)
    B = _hsl2rgb_h(H, temp1, temp2, H-1.0/3.0)
    return R, G, B



BREWER_QUALITATIVE_SET1 = [
    hex2rgb("#D23329"),  # red
    hex2rgb("#467CB4"),  # blue
    hex2rgb("#69AD59"),  # green
    hex2rgb("#8F559D"),  # purple
    hex2rgb("#E97F32"),  # orange
    hex2rgb("#FBFD60"),  # yellow
    hex2rgb("#975C2C"),  # brown
    hex2rgb("#E383BA"),  # pink
    ]

