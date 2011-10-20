def testattr(**kwargs):
    """Add attributes to a test function/method/class.
    
    This function is needed to be able to add
      @attr(slow = True)
    for functions.
    
    """
    def wrap(func):
        func.__dict__.update(kwargs)
        return func
    return wrap
