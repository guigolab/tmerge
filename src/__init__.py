from .merge import Merge
from pkg_resources import iter_entry_points

discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point
    in iter_entry_points('tmerge.plugins')
}

def tmerge(input_path, output_path, tolerance, processes, **kwargs):
    merger = Merge(input_path, output_path, tolerance, processes)

    for plugin in discovered_plugins.values():
        plugin(merger, input_path, output_path, tolerance, processes, **kwargs)

    merger.merge()