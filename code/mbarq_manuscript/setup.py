from setuptools import setup

setup(
    name="mbarq_manuscript_workflow",
    version="0.0.1",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("mBARq manuscript workflow"),
    license="LICENSE",
    keywords="mbarq_man",
    #url="",
    install_requires=[
        'click', 'pyaml'],
    packages=['mbarq_manuscript_workflow'],
    entry_points={
        'console_scripts': ['mflow=mbarq_manuscript_workflow.main:main'],
    }
)

