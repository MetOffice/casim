#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

'''
This module contains a Fab-based build script for CASIM.
'''

import logging
from pathlib import Path

from fab.fab_base.fab_base import FabBase
from fab.api import Category, grab_files


class FabCasim(FabBase):
    '''
    A class to build CASIM using Fab as base class.

    :param str name: name of the build.
    '''

    def __init__(self, name):
        super().__init__(name, link_target="static-library")
        # Store the root directory of MONC:
        self._root = Path(__file__).resolve().parents[1]

        # Casim uses UM's vectorlib, and that seems to need 8-byte-integers.
        # So add the flag to the compiler
        fc = self.config.tool_box.get_tool(Category.FORTRAN_COMPILER)
        fc.add_flags(fc["default-8-byte-integer"], "base")

    def grab_files_step(self) -> None:
        '''
        Extracts all the required source files from the repositories.
        It then sets the include path (since include files are not
        copied into the build tree).

        :raises RuntimeError: the expected `rose-meta/um-atmos` file
            does not exist, indicating an invalid directory structure.
        '''
        grab_files(self.config, self._root / "src")

    def define_preprocessor_flags_step(self) -> None:
        '''
        Defines the preprocessor flags.
        '''
        super().define_preprocessor_flags_step()

        flags = ['-DDEF_MODEL=MODEL_MONC', '-DMODEL_MONC=4']

        self.add_preprocessor_flags(flags)


# ==========================================================================
if __name__ == "__main__":

    # Enable full fab logging for now:
    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s: %(name)s: %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    fab_casim = FabCasim("casim")
    fab_casim.build()
