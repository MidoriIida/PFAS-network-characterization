for name in list(globals().keys()):
    if not name.startswith('_'):
        del globals()[name]
