from setuptools import setup, find_packages
import sys

__version__ = '2.0.1'

install_requires = [
    'docopt', 
    'numpy',
    'geopy',
    ]

qcviz_requires=[
    'matplotlib',
    'ipython'
    ]

tests_requires= [
    'nose',
    ]

setup(name='qccodar3',
      version=__version__,
      description="Apply quality controls to improve CODAR output",
      long_description="",
      classifiers=[
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    ],
      keywords='qc codar radials currents ocean science data',
      author='Sara Haines, Teresa Updyke',
      author_email='sarahaines@unc.edu',
      maintainer='Teresa Updyke',
      maintainer_email='tupdyke@odu.edu',

      url='http://nccoos.unc.edu',
      license='GNU',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      include_package_data=True,
      zip_safe=False,
      namespace_packages=["qccodar3", "qccodar3.qcviz", "qccodar3.test"],
      install_requires=install_requires,
      extras_require={
        'tests' : tests_requires,
        'qcviz' : qcviz_requires
        },
      test_suite="qccodar3.test",
      entry_points="""
        [console_scripts]
        qccodar3 = qccodar3.app:main
      """,
      data_files=[('config', ['src/qccodar3/config/qccodar.plist']),
                  ('file_formats',['src/qccodar3/file_formats/radialshort_LLUV_RDL7.ruv'])],
      #scripts=['qccodar3/qcviz/qcviz.py']
      )