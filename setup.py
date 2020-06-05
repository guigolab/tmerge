from setuptools import setup, find_packages
import os

with open(os.path.join("../", os.path.abspath(os.path.dirname(__file__)), 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="tmerge",
    version="2.0.0",
    description="Merge transcriptome read-to-genome alignments into non-redundant transcript models.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/jacobwindsor/tmerge",
    author="Jacob Windsor",
    author_email="me@jcbwndsr.com",
    license="MIT",
    keywords="tmerge merge transcriptomics stringtie2 rna-seq rnaseq long-read-rna long-read transcriptome",
    packages=find_packages(include=["src", "plugins"]),
    install_requires=["pyfaidx", "typing"],
    python_requires='>=3',
    entry_points={
        "console_scripts": [
            "tmerge = cli:main"
        ]
    }
)