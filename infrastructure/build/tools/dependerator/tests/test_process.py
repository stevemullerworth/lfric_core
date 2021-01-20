#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from pathlib import Path
import pytest

from dependerator.database import (FileDependencies,
                                   FortranDependencies,
                                   SQLiteDatabase)
from dependerator.process import FortranProcessor


class TestNamelistDescriptionAnalyser():
    @pytest.fixture
    def databases(self, tmp_path_factory):
        filename = tmp_path_factory.mktemp('db-', True) / 'fortran.db'
        database = SQLiteDatabase(filename)
        return FortranDependencies(database), FileDependencies(database)

    def test_compile_dependencies(self, databases):
        self._populate_database(databases[0])

        uut = FortranProcessor(databases[0], "objects", "modules")
        uut.determineCompileDependencies(databases[1])

        assert [(Path('objects/bits/bar.o'), [Path('modules/bits/baz.mod')]),
                (Path('objects/bobs/qux.o'), [Path('modules/bits/baz.mod')]),
                (Path('objects/foo.o'), [Path('modules/bits/bar.mod')])] \
            == list(databases[1].getDependencies())

    def test_link_dependencies(self, databases):
        self._populate_database(databases[0])

        import sys
        uut = FortranProcessor(databases[0], "objects", "modules")
        result = list(uut.determineLinkDependencies())

        assert [(u'objects/foo', ['objects/bits/bar.o',
                                  'objects/bits/baz.o',
                                  'objects/bobs/qux.o',
                                  'objects/foo.o'])] == result

    def _populate_database( self, database: FortranDependencies ):
        database.addProgram( 'foo', 'foo.f90' )
        database.addModule( 'bar', 'bits/bar.f90' )
        database.addModule( 'baz', 'bits/baz.f90' )
        database.addModule( 'qux', 'bobs/qux.f90' )

        database.addModuleCompileDependency( 'foo', 'bar' )
        database.addModuleCompileDependency( 'bar', 'baz' )
        database.addModuleCompileDependency( 'qux', 'baz' )

        database.addModuleLinkDependency( 'foo', 'bar' )
        database.addModuleLinkDependency( 'bar', 'baz' )
        database.addModuleLinkDependency( 'baz', 'qux' )

if __name__ == '__main__':
    unittest.main()
