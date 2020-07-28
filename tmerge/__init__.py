from pkg_resources import iter_entry_points
from .merge import Merge as _Merge
from . import plugins

_discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point
    in iter_entry_points('tmerge.plugins')
}

def merge(**kwargs):
    input_path = kwargs.get("input_path")
    output_path = kwargs.get("output_path")
    tolerance = kwargs.get("tolerance", 0)
    processes = kwargs.get("processes", None)
    plugins = kwargs.get("plugins", [])
    
    merger = _Merge(input_path, output_path, tolerance, processes)

    for plugin in _discovered_plugins.values():
        plugin(merger.hooks, **kwargs)

    for plugin in plugins:
        plugin(merger.hooks, **kwargs)

    merger.merge()

