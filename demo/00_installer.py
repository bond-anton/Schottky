from __future__ import division, print_function

from BDProjects.Client import Connector, Installer

installer = Installer(Connector(config_file_name='config.ini'), overwrite=True)
