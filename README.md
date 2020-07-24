# tmerge [![Build Status](https://travis-ci.com/jacobwindsor/tmerge.svg?branch=master)](https://travis-ci.com/jacobwindsor/tmerge)
Build highly accurate full-length transcripts from third generation sequencing alignments.

tmerge compares transcript structures (or read-to-genome alignments) present in the input and attempts to reduce transcript redundancy, i.e., merge compatible input transcripts into non-redundant transcript models.

tmerge is fast and can typically process several millions of aligned long reads in a few minutes.

![tmerge](/images/tmerge2_basic.png)

## Installation
```
pip install tmerge
```

It is recommended to install tmerge within a [virtual environment](https://docs.python.org/3/tutorial/venv.html)


## Usage
tmerge offers both a CLI and a Python module. The CLI is built upon the Python module and includes several built in "plugins" (see below for description of plugins).

### CLI
Once you have installed tmerge via pip, the CLI will be available on your PATH. If you have installed tmerge into a virtual environment, you will need to activate that virtual environment to run the CLI.

Run `tmerge --help` for a description of the options.

### Python module
You may import tmerge as a Python module and call the `merge` function to run it. `merge` takes 2 mandatory kwargs (`input_path` and `output_path`) and two optional kwargs (`tolerance` and `processes`). Any additional keyword arguments will be sent to any registered plugins (see below). 

```python
from tmerge import merge
input_path = "my.gff"
output_path = "output.gff"

merge(input_path=input_path, output_path=output_path) # Will block until completion

# Can now do more things with the output file
with open(output_path, "a") as f:
    f.write("# A comment \n")
    f.flush()

```

## Plugins
Plugins allow you to "hook" into tmerge's lifecycle events and allow you to view, edit or remove the transcripts passing through tmerge and adapt it to your lab's specific needs. For example, adding Hi-Seq support.

This section explains how to write plugins and register them to tmerge. 

### Transcript Models and Contigs
Before writing a plugin, it is important to understand the concept of Transcript Models and Contigs. Transcripts are represented in tmerge as `TranscriptModel` objects, at first these are the transcripts defined in the input file but are altered throughout the lifecycle of tmerge, either having other transcript models merged into them or removed entirely. Contigs are lists of overlapping transcript models. Merging of transcript models is only performed within a contig and not between contigs.

### Plugin class
A plugin is a simple class that registers itself to one or more "hooks" in it's init method. It receives the `hooks` dict as it's first argument followed by all of the kwargs that are passed to `tmerge.merge`.

```python
class Counter:
    def __init__(self, hooks,**kwargs):
        self.count = 0
        hooks["transcript_added"].tap(self.add_one)

    def add_one(self, *args):
        self.count += 1

    def print(self, *args);
        print(f"There are {self.count} transcripts")
```

### Editing transcripts
Some of the hooks send transcripts to the hooked-in function (see table below). You can edit or remove any of these transcripts and changes will be reflected in the output merged file. Further, any key/value pairs added to the `meta` dict will be appended to the "attributes" column of the output merged GFF.

```python
class MyPointlessPlugin:
    def __init__(self, hooks, extra_attribute, bad_id, **kwargs):
        self.extra_attribute = extra_attribute
        self.bad_id = bad_id

        hooks["transcript_added"].tap(self.add_meta)
        hooks["contig_merged"].tap(self.remove_if_matches)

    def add_meta(self, transcript, *args):
        # When tmerge.merge(input_path=output, output_path=output, extra_attribute="Pointless") is ran 'extra: "Pointless"' will be added to the attributes column for every transcript
        self.transcript.meta["extra"] = self.extra_attribute

    def remove_if_matches(self, contig, *args):
        # Running tmerge.merge(input_path=output, output_path=output, bad_id="bad") will remove any transcript with the id of "bad" from the result
        for transcript in contig:
            if transcript.id = self.bad_id:
                transcript.remove() # Flags a transcript for removal
```

### Registering plugins
#### Simple list
The easiest way to provide tmerge with plugins is to pass the `plugins` kwarg to `tmerge.merge`.

```python
from myplugins import MySimplePlugin, MyAdvancedPlugin
from tmerge import merge

merge(
    input_path="input.gff",
    output_path="output.gff",
    plugins=[
        MySimplePlugin,
        MyAdvancedPlugin
    ]
)
```

#### Dynamic Plugin Discovery
If you're already using setup_tools in your project, then you can use [dynamic plugin discovery](https://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins) to easily drop in plugins to tmerge.

In your project's `setup.py` add your plugin to the `tmerge.plugins` group:

```python
# setup.py
setup(
    ...
    entry_points={
        "tmerge.plugins": "plugin_name = my.plugin.module.MyPlugin"
    }
)
```

This will automatically register your plugin with tmerge and the plugin will be executed with `tmerge.merge`.

### Lifecycle events
![tmerge2 hookes](/images/tmerge2_hooks.png)


| Hook Name          | When?                                                                                        | Arguments sent to hooked-in functions                                                                                                                               |
|--------------------|----------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chromosome_parsed  | When one chromosome is parsed from the input                                                 | `chromosome` (list of `TranscriptModel`s)                                                                                                                           |
| transcript_added   | When a transcript is added to a contig                                                       | `transcript` (`TranscriptModel`)                                                                                                                                    |
| contig_built       | When one contig (group of overlapping transcripts) is built                                  | `contig` (list of `TranscriptModel`s)                                                                                                                               |
| transcripts_merged | When one transcript is merged into another                                                   | `host_transcript` (`TranscriptModel`), `merged_transcript` (`TranscriptModel`)  `host_transcript` is the transcript that has had `merged_transcript` merged into it |
| contig_merged      | When one contig is fully merged                                                              | `contig` (list of `TranscriptModel`s)                                                                                                                               |
| contig_complete    | Contig has been fully merged, transcript flagged for removal removed, and queued for writing | `contig` (list of `TranscriptModel`s)                                                                                                                               |
| merging_complete   | All transcripts have been merged                                                             | None                                                                                                                                                                |
| pre_sort           | Just before the merged output is sorted                                                      | None                                                                                                                                                                |
| post_sort          | Just after the merged output is sorted                                                       | None                                                                                                                                                                |
| complete           | Everything complete                                                                          | None                                                                                                                                                                |
### Examples
See the `plugins/` folder for examples of various plugins.


## Authors
Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

Jacob Windsor, CRG, Barcelona, contact me@jcbwndsr.com
