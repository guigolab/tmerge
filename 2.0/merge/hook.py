class Hook:
    def __init__(self):
        self.tapped_funcs = []

    def exec(self, *args):
        # Execute all functions that are tapped to this hook
        for tapped in self.tapped_funcs:
            tapped(*args)

    def tap(self, func):
        # Add a function to be executed when this hook is called
        self.tapped_funcs.append(func)