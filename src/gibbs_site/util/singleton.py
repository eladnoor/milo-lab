def Singleton(cls):
    """Decorator that makes a class a singleton.
    
    TODO(flamholz): move somewhere more central.
    """
    instance_container = []
    def getinstance():
        if not len(instance_container):
            instance_container.append(cls())
        return instance_container[0]
    return getinstance