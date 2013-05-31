class FatalError:
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class RerunnableError:
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
