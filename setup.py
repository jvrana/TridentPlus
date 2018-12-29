from distutils.core import setup

setup(name='tridentplus',
      packages=['tridentplus', 'tridentplus.utils'],
      package_data={
          'tridentplus': ['data/*']
      },
      version='0.0.1'
     )
