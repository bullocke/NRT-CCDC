from setuptools import setup

setup(name='nrtm',
      version='0.1',
      description='Near-real time forest monitoring',
      url='http://bullocke.github.io/NRT',
      author='Eric Bullock',
      author_email='bullocke@bu.edu',
      license='MIT',
      packages=['nrtm'],
      zip_safe=False,
      entry_points='''
      [console_scripts]
      nrtm=nrtm.cli.main:cli
      [nrtm.nrtm_commands]
      monitor=nrtm.cli.monitor:monitor
      ''',
      )
