from setuptools import setup

setup(
    name="rbseq_workflow",
    version="0.0.1",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("RBSeq analysis workflow"),
    license="LICENSE",
    keywords="rbseq",
    #url="",
    install_requires=[
        'click', 'pyaml'],
    packages=['rbseq_workflow'],
    entry_points={
        'console_scripts': ['rbseq=rbseq_workflow.main:main'],
    }
)

