def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(rgb):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (rgb[0], rgb[1], rgb[2])

def color_blind(n = None, alpha = None):
    """Returns a color blind friendly color palette
    """
    from aei import objects
    base = ["#E69F00", "#009E73", "#0072B2", "#D55E00", "#F0E442", 
        "#56B4E9", "#CC79A7", "#000000"]
    
    # create the baseline object
    cp = objects.color(base)
    
    # interploate the colors if set
    if n:
        cp.interpolate(n)
        
    # add alpha values if set
    if alpha:
        cp.add_alpha(alpha)
        
    return cp.palette
    
def stanford(n = None, alpha = None):
    """Returns the stanford base color palette
    """
    from aei import objects
    base = ["#8C1515", "#4D4F53", "#000000"]
    
    # create the baseline object
    cp = objects.color(base)
    
    # interploate the colors if set
    if n:
        cp.interpolate(n)
        
    # add alpha values if set
    if alpha:
        cp.add_alpha(alpha)
        
    return cp.palette
        
def stanford_bright(n = None, alpha = None):
    """Returns the stanford bright color palette
    """
    from aei import objects
    base = ["#0098DB", "#007C92", "#009B76", "#B26F16", "#E98300", "#EAAB00"]
    
    # create the baseline object
    cp = objects.color(base)
    
    # interploate the colors if set
    if n:
        cp.interpolate(n)
        
    # add alpha values if set
    if alpha:
        cp.add_alpha(alpha)
        
    return cp.palette
        
def stanford_dark(n = None, alpha = None):
    """Returns the stanford dark color palette
    """
    from aei import objects
    base = ["#5E3032", "#8D3C1E", "#00505C", "#175E54", "53284F"]
    
    # create the baseline object
    cp = objects.color(base)
    
    # interploate the colors if set
    if n:
        cp.interpolate(n)
        
    # add alpha values if set
    if alpha:
        cp.add_alpha(alpha)
        
    return cp.palette
        
def stanford_neutral(n = None, alpha = None):
    """Returns the stanford neutral color palette
    """
    from aei import objects
    base = ["#D2C295", "#B3995D", "#928B81", "#B6B1A9", 
        "#DAD7CB", "#C7D1C5", "#4D5357", "#5F574F"]
    
    # create the baseline object
    cp = objects.color(base)
    
    # interploate the colors if set
    if n:
        cp.interpolate(n)
        
    # add alpha values if set
    if alpha:
        cp.add_alpha(alpha)
        
    return cp.palette