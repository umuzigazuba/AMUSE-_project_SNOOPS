from setuptools import find_packages, setup

setup(
    name = 'src',
    packages = find_packages(),
    version = '0.1.0',
    description = 'The formation of a second generation of stars in a globular cluster after a collision with a molecular cloud',
    author = 'Yiqi Wu, Kostas Tsalapatas, Erin Umuzigazuba',
    license = 'MIT',
    url = 'https://github.com/umuzigazuba/AMUSE_project_YESSIR/tree/Yiqi',
    download_url = 'https://github.com/umuzigazuba/AMUSE_project_YESSIR.git',
    install_requires=[
        'amuse-framework',
        'amuse-seba',
        'amuse-sse',
        'amuse-bhtree',
        'amuse-fi',
        'natsort',
        'numpy',
        'matplotlib',
        'plotly'
    ],
)
