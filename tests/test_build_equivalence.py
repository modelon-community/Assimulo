#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2026 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Bit-for-bit equivalence of this build with a reference build.

    python tools/dump_examples.py reference.npz          # with the reference build installed
    ASSIMULO_REFERENCE_RESULTS=reference.npz pytest tests/test_build_equivalence.py

Skipped when no reference file is given. The reference is specific to the machine and the
compilers (the point is to compare two builds of the same sources on the same machine, e.g.
the numpy.distutils build against the meson-python build).
"""
import os
import sys

import pytest

TOOLS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools")


@pytest.mark.skipif(not os.environ.get("ASSIMULO_REFERENCE_RESULTS"),
                    reason="set ASSIMULO_REFERENCE_RESULTS to a dump made with tools/dump_examples.py")
def test_examples_bit_identical_to_reference(tmp_path):
    sys.path.insert(0, TOOLS)
    import dump_examples
    out, skipped = dump_examples.run_all()
    assert not skipped, skipped
    current = str(tmp_path / "current.npz")
    dump_examples.save(current, out, skipped)
    assert dump_examples.compare(os.environ["ASSIMULO_REFERENCE_RESULTS"], current)
