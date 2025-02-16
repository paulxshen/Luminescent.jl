from setuptools import setup, find_packages

setup(
    name='luminescent',  # Your package name
    version='0.5.1',  # Your package version
    description='A description of your package',
    author='Paul Shen',
    author_email='pxshen@alumni.stanford.edu',
    packages=find_packages(),  # Automatically find your package(s)
    install_requires=[
        "gdsfactory",
        "pillow",
        "pymeshfix",
        "electromagneticpython",
        "sortedcontainers",
        'scikit-rf',
        # 'python-opencv',
    ],
)
# cd C:\Users\pxshe\OneDrive\Desktop\beans\Luminescent.jl\luminescent
# python -m build
# twine upload dist/*

# pip install gdsfactory pillow pymeshfix electromagneticpython sortedcontainers scikit-rf
