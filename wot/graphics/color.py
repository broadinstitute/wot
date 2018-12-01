import matplotlib


def hexstring_of_rgba(c):
    return '#%02x%02x%02x%02x' % tuple(int(255 * x) for x in c)


def hexstring_of_hsv(c):
    return make_transparent(matplotlib.colors.hsv_to_rgb(c), 1)


def make_transparent(color, alpha_channel):
    return hexstring_of_rgba(matplotlib.colors.to_rgba(color, alpha_channel))


# Color generation

def color_heatmap(x):
    """Get the heatmap color of a value in [0, 1] : the higher the hotter"""
    return hexstring_of_hsv((.7 * (1 - x), 1.0, .92))


def color_mix(a, b, p):
    return tuple(a[i] * p + b[i] * (1 - p) for i in range(len(a)))


def color_linear_gradient(f, t, steps):
    return [color_mix(f, t, i / (steps - 1)) for i in range(steps)]
