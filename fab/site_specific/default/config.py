#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

'''
This module contains the default Baf configuration class.
'''

import argparse
from typing import cast, List, Optional

from fab.api import (BuildConfig, Category, Compiler, ProfileFlags,
                     ToolRepository)


class Config:
    '''
    This class is the default Configuration object for Fab builds.
    It provides several callbacks which will be called from the build
    scripts to allow site-specific customisations.
    '''

    def __init__(self) -> None:
        self._args: argparse.Namespace

    @property
    def args(self) -> argparse.Namespace:
        '''
        :returns: the command line options specified by the user.
        '''
        return self._args

    def get_valid_profiles(self) -> List[str]:
        '''
        Determines the list of all allowed compiler profiles. The first
        entry in this list is the default profile to be used. This method
        can be overwritten by site configs to add or modify the supported
        profiles.

        :returns: List of all supported compiler profiles.
        '''
        return ["debug", "normal"]

    def handle_command_line_options(self, args: argparse.Namespace) -> None:
        '''
        Additional callback function executed once all command line
        options have been added. This is for example used to add
        Vernier profiling flags, which are site-specific.

        :param argparse.Namespace args: the command line options added in
        the site configs
        '''
        # Keep a copy of the args, so they can be used when
        # initialising compilers
        self._args = args

    def update_toolbox(self, build_config: BuildConfig) -> None:
        '''
        Set the default compiler flags for the various compiler
        that are supported.

        :param build_config: the Fab build configuration instance
        '''
        # Define a base profile, which contains the common
        # compilation flags. This 'base' is not accessible to
        # the user, so it's not part of the profile list.
        ProfileFlags.define_profile("base")
        for profile in self.get_valid_profiles():
            ProfileFlags.define_profile(profile, inherit_from="base")

        # Note that CASIM needs UM .mod files and DrHook in order to compile.
        # E.g. you will have to add:
        # gfortran.add_flags(
        #      ['-I', "PATH-TO-UM/build_output/"], "base")
        # gfortran.add_flags(
        #      ['-I', "PATH-TO-UM/build_output/"], "base")

        self.setup_intel_classic(build_config)
        self.setup_intel_llvm(build_config)
        self.setup_gnu(build_config)
        self.setup_nvidia(build_config)
        self.setup_cray(build_config)

    @staticmethod
    def get_compiler(name: str) -> Optional[Compiler]:
        """
        Searched for a compiler with the specified name. If the
        compiler is not found, it will try to use mpif90 (in case that
        e.g. gfortran is not in PATH, but mpif90-gfortran is).

        :returns: the Fab compiler object, or None if the compiler is
            not available.
        """
        tr = ToolRepository()
        compiler = tr.get_tool(Category.FORTRAN_COMPILER, name)
        if not compiler.is_available:
            try:
                compiler = tr.get_tool(Category.FORTRAN_COMPILER,
                                       f"mpif90-{name}")
                if not compiler.is_available:
                    return None
            except KeyError:
                # mpif90-NAME does not exist.
                return None
        return cast(Compiler, compiler)

    def setup_cray(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Cray compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        ftn = Config.get_compiler("crayftn-ftn")
        if not ftn:
            return

        ftn.add_flags(["-e", "m"], "base")
        ftn.add_flags(["-g",
                       "-Ktrap=divz,inv,ovf",    # floating point checking
                       "-R", "bcdps",  # bounds, array shape, collapse,
                                       # pointer, string checking
                       "-O0"],         # No optimisation
                      "debug")
        ftn.add_flags(["-O3"], "normal")

    def setup_gnu(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Gnu compiler and linker flags.
        For now call an external function, since it is expected that
        this configuration can be very lengthy (once we support
        compiler modes).

        :param build_config: the Fab build configuration instance
        '''
        gfortran = Config.get_compiler("gfortran")
        if not gfortran:
            return

        gfortran.add_flags(["-fbacktrace"], "base")
        gfortran.add_flags(["-g", "-fcheck=all",
                            "-ffpe-trap=invalid,zero,overflow"], "debug")
        gfortran.add_flags(["-O3"], "normal")

    def setup_intel_classic(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Intel classic compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        ifort = Config.get_compiler("ifort")
        if not ifort:
            return

        ifort.add_flags(["-g", "-traceback"], "base")
        ifort.add_flags(["-check bounds,uninit", "-no-vec"], "debug")
        ifort.add_flags(["-O2"], "normal")

    def setup_intel_llvm(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Intel LLVM compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        ifx = Config.get_compiler("ifx")
        if not ifx:
            return

        ifx.add_flags(["-g", "-traceback"], "base")
        ifx.add_flags(["-check bounds,uninit", "-no-vec"], "debug")
        ifx.add_flags(["-O2"], "normal")

    def setup_nvidia(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Nvidia compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        nvfortran = Config.get_compiler("nvfortran")
        if not nvfortran:
            return

        nvfortran.add_flags(["-g", "-traceback", "-O0"], "base")
        nvfortran.add_flags(["-fp-model=strict"], "debug")
        nvfortran.add_flags(["-O4"], "normal")
