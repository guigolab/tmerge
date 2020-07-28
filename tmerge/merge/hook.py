class Hook:
    """
    Defines a "hook" that is used to hook into tmerge lifecycle events.
    """
    def __init__(self):
        self.tapped_funcs = []

    def exec(self, *args):
        """
        Executes all of the functions that are tapped into the hook instance.
        """
        # Execute all functions that are tapped to this hook
        for tapped in self.tapped_funcs:
            tapped(*args)

    def tap(self, func):
        """
        Adds a function to the hook instance.
        """
        # Add a function to be executed when this hook is called
        self.tapped_funcs.append(func)