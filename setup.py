from setuptools import setup, find_packages, findall

setup(  name                = 'diffuse',
        version             = '0.1',
        description         = 'A Python package for simulating diffuse scattering of X-rays by protein crystals',
        author              = 'J A van der Horn',
        author_email        = 'jitsie@gmail.com',
        url                 = '',
        license             = '',
        install_requires    = ( ),
        package_dir         = {'':'scud'},
        packages            = find_packages(where   ='diffuse-python',
                                            exclude = ( )
                                            ),
        include_package_data= True,
        classifiers         = [ ],
        )
