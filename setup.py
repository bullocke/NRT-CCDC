from setuptools import setup

setup(name='cfm',
      version='0.1',
      description='Near-real time forest monitoring',
      url='http://bullocke.github.io/NRT',
      author='Eric Bullock',
      author_email='bullocke@bu.edu',
      license='MIT',
      packages=['cfm'],
      zip_safe=False,
      entry_points='''
      [console_scripts]
      cfm=cfm.cli.main:cli
      [cfm.cfm_commands]
      monitor=cfm.cli.monitor:monitor
      monitor_map=cfm.cli.monitor_map:monitor_map
      ''',
      )
