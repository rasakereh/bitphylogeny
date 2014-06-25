from distutils.core import setup
 
setup(
      name='BitPhylogeny',
      version='0.1',
      description='Bayesian inference of intra tumor phylogenies.',
      author='Ke Yuan',
      author_email='ke.yuan.09@gmail.com',
      url='https://bitbucket.org/ke_yuan/bitphylogeny',
      package_dir = {'': 'lib'},    
      packages=['bitphylogeny', 'bitphylogeny.post_process'],
      scripts=['bin/BitPhylogeny']
     )
