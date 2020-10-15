def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip("#")
    lv = len(value)
    return tuple(int(value[i : i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb):
    """Return color as #rrggbb for the given color values."""
    return "#%02x%02x%02x" % (rgb[0], rgb[1], rgb[2])


def color_blind(n=None, alpha=None, black=False, return_cmap=False):
    """Returns a color blind friendly color palette"""
    from earthtools import objects

    # set metadata
    name = "color_blind"
    palette = ["#0072B2", "#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#E69F00", "#D55E00", "#000000"]

    # remove the black value if set
    if not black:
        palette.remove("#000000")

    # create the baseline object
    cp = objects.color(palette=palette, n=n, name=name)

    # add alpha values if set
    if alpha:
        cp.alpha(alpha)

    # default to return the palette, otherwise to return the color map
    if return_cmap:
        return cp.cmap
    else:
        return cp.palette


def stanford(n=None, alpha=None, return_cmap=False):
    """Returns the stanford base color palette"""
    from earthtools import objects

    # set metadata
    name = "stanford_primary"
    palette = ["#8C1515", "#4D4F53", "#000000"]

    # create the baseline object
    cp = objects.color(palette=palette, n=n, name=name)

    # add alpha values if set
    if alpha:
        cp.alpha(alpha)

    # default to return the palette, otherwise to return the color map
    if return_cmap:
        return cp.cmap
    else:
        return cp.palette


def stanford_bright(n=None, alpha=None, return_cmap=False):
    """Returns the stanford bright color palette"""
    from earthtools import objects

    # set metadata
    name = "stanford_bright"
    palette = ["#0098DB", "#007C92", "#009B76", "#B26F16", "#E98300", "#EAAB00"]

    # create the baseline object
    cp = objects.color(palette=palette, n=n, name=name)

    # add alpha values if set
    if alpha:
        cp.alpha(alpha)

    # default to return the palette, otherwise to return the color map
    if return_cmap:
        return cp.cmap
    else:
        return cp.palette


def stanford_dark(n=None, alpha=None, return_cmap=False):
    """Returns the stanford dark color palette"""
    from earthtools import objects

    # set metadata
    name = "stanford_dark"
    palette = ["#5E3032", "#8D3C1E", "#00505C", "#175E54", "#53284F"]

    # create the baseline object
    cp = objects.color(palette=palette, n=n, name=name)

    # add alpha values if set
    if alpha:
        cp.alpha(alpha)

    # default to return the palette, otherwise to return the color map
    if return_cmap:
        return cp.cmap
    else:
        return cp.palette


def stanford_neutral(n=None, alpha=None, return_cmap=False):
    """Returns the stanford neutral color palette"""
    from earthtools import objects

    # set metadata
    name = "stanford_neutral"
    palette = ["#D2C295", "#B3995D", "#928B81", "#B6B1A9", "#DAD7CB", "#C7D1C5", "#4D5357", "#5F574F"]

    # create the baseline object
    cp = objects.color(palette=palette, n=n, name=name)

    # add alpha values if set
    if alpha:
        cp.alpha(alpha)

    # default to return the palette, otherwise to return the color map
    if return_cmap:
        return cp.cmap
    else:
        return cp.palette
