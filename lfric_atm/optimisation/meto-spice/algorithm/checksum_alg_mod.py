##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''PSyclone transformation script for the Dynamo0p3 API to apply
colouring and redundant computation to the level1 halo for
setval_c generically.

The inner product used by the checksum algorithm produces different
results on different numbers of threads if OpenMP is applied. It means
the checksums are different even if the model result is the
same. Therefore, we do not apply the OpenMP transform to the
algorithm.

'''
from psyclone.transformations import Dynamo0p3ColourTrans, \
                                     Dynamo0p3RedundantComputationTrans

from psyclone.domain.lfric import LFRicConstants


def trans(psy):
    '''Applies PSyclone colouring and redundant computation
    transformations.

    '''
    ctrans = Dynamo0p3ColourTrans()
    rtrans = Dynamo0p3RedundantComputationTrans()
    const = LFRicConstants()

    setval_count = 0
    # Loop over all of the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:

        print("Transforming invoke '{0}' ...".format(invoke.name))
        schedule = invoke.schedule

        # Make setval_c compute redundantly to the level 1 halo if it
        # is in its own loop.
        for loop in schedule.loops():
            if loop.iteration_space == "dof":
                if len(loop.kernels()) != 1:
                    raise Exception(
                        "Expecting loop to contain 1 call but found '{0}'".
                        format(len(loop.kernels())))
                if loop.kernels()[0].name in ["setval_c"]:
                    setval_count += 1
                    rtrans.apply(loop, options={"depth": 1})

        # Colour loops over cells unless they are on discontinuous
        # spaces or over dofs
        for loop in schedule.loops():
            if loop.iteration_space == "cell_column" \
                and loop.field_space.orig_name \
                    not in const.VALID_DISCONTINUOUS_NAMES:
                ctrans.apply(loop)

        # Take a look at what we've done
        print("Found {0} setval calls".format(setval_count))
        schedule.view()

    return psy
