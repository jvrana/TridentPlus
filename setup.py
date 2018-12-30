from distutils.core import setup


tests_require = [
    'pytest',
]

install_requires = [
    'pandas',
    'tqdm',
    'dill',
    'intermine',
    'primer3-py',
    'jdna',
    'pydent',
    'pyblast',
    'primer3plus'
]

setup(name='tridentplus',
      packages=['tridentplus', 'tridentplus.utils'],
      package_data={
          'tridentplus': ['data/*']
      },
      version='0.0.1',
      install_requires=install_requires,
      tests_require=tests_require,
     )
