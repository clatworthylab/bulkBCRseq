#!/usr/bin/env python
import sh


def test_script(option=None, metadata=None, bsub=True, verbose=True, execute=False):
    """
    Test to print the bsub commands

    Parameters
    ----------
    option : int
        options for running the main script
    metadata : str, optional
        path to metadata.txt file
    bsub : bool
        whether or not use bsub command
    verbose : bool
        whether or not to print command
    execute : bool
        where or not to run command

    """

    if option is None:
        opt = 1
    else:
        opt = option

    if metadata is None:
        meta = 'tests/data/Sample_metadata.txt'
    else:
        meta = metadata

    if bsub:
        bsub_ = 'Y'
    else:
        bsub_ = 'N'

    if verbose:
        verbose_ = 'Y'
    else:
        verbose_ = 'N'

    if execute:
        execute_ = 'Y'
    else:
        execute_ = 'N'

    cmd = ['Processing_sequences_large_scale.py',
           meta,
           str(opt),
           bsub_,
           verbose_,
           execute_]
    try:
        # run the shell command
        sh.python(cmd)
    except sh.ErrorReturnCode as e:
        # print the error, so we know what went wrong
        print e
        # make sure the test fails
        pytest.fail(e)


if __name__ == "__main__":
    # print bsub commands
    test_script()
    test_script(2)
    test_script(3)
    test_script(4)
    # print non-bsub commands
    test_script(bsub=False)
    test_script(2, bsub=False)
    test_script(3, bsub=False)
    test_script(4, bsub=False)
