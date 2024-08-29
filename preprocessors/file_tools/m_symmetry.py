
"""
this file is part of the phonon explorer package!
author: tyler c. sterling & dmitry reznik.
email: ty.sterling@colorado.edu
affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
date: 08/29/2024
description:
    tools to check symmetry args being passed to mantid
"""

from file_tools.m_file_utils import *

# --------------------------------------------------------------------------------------------------

def _do_space_group_stuff(space_group):

    """
    print symmetry operations from mantid
    """

    print('\nfound valid space-group arg.')
    symm_ops = space_group.getSymmetryOperationStrings()

    print('\nsymmetry-operations:')
    for symm_op in symm_ops:
        print(' ',symm_op)

    print('')

# --------------------------------------------------------------------------------------------------

def _do_symm_ops_stuff(symm_ops):

    """
    print symmetry operations from mantid
    """

    print('\nfound valid symmetry-operation arg.')
    print('\nsymmetry-operations:')
    for symm_op in symm_ops:
        print(' ',symm_op.getIdentifier())

    print('')

# --------------------------------------------------------------------------------------------------

def check_symmetry_args(SymmetryArgs):

    """
    do some pre-processing of symmetry args to avoid wasting time loading nxspe data etc.
    crash now if args are likely to cause mantid to crash later.
    """

    # do nothing if default
    if SymmetryArgs is None:
        return

    print('\n*** checking symmetry args ***')
    print('\nSymmetryArgs:',SymmetryArgs)

    # otherwise, use mantid api to check args
    from mantid.geometry import SpaceGroupFactory, SymmetryOperationFactory
    
    # try to pass args along to relevant class
    try:
        space_group = SpaceGroupFactory.createSpaceGroup(SymmetryArgs)
        _do_space_group_stuff(space_group)
        return
    except Exception as ex:
        print('\nnot a proper space-group arg, here is the exception:')
        print(' ',ex)

    try:
        symm_ops = SymmetryOperationFactory.createSymOps(SymmetryArgs)
        _do_symm_ops_stuff(symm_ops)
        return
    except Exception as ex:
        print('\nnot a proper symmetry-operation arg, here is the exception:')
        print(' ',ex)

    msg = '\nwrong SymmetryArgs. give either space-group name, or list of symmetry operations\n'
    msg += '  e.g. SymmetryArgs=\'x,y,z;-x,y,z;x,-y,z;x,y,-z\' \n'
    msg += 'or one of the following space-group names:'
    print(msg)
    
    symbols = SpaceGroupFactory.getAllSpaceGroupSymbols()
    for symbol in symbols:
        print(' ',symbol)

    crash()

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    
    SymmetryArgs = 'P 4'
    check_symmetry_args(SymmetryArgs)
 

