from setuptools import setup, find_packages

setup(  name                = 'diffuse',
        version             = '0.1',
        description         = 'A Python package for extracting and simulating diffuse scattering of X-rays in protein crystals',
        author              = 'J A van der Horn',
        author_email        = 'jitsie@gmail.com',
        url                 = '',
        license             = 'GPLv3',
        install_requires    = ( ),
        package_dir         = {'diffuse':'src'},
        packages            = find_packages(where='diffuse'),
        include_package_data= True,
        classifiers         = ['Programming Language :: Python :: 3'],
        )
