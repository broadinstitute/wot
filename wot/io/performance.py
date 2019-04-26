import os
from math import log10


def output_progress(count, total=1.0):
    """
    Prints a nice progress bar when doing long computations.

    Parameters
    ----------
    count : float
        Current progress.
    total : float
        Progress when completed. 1 by default.

    Note
    ----
    If total is more than 1, progress is displayed as "X / Y".
    If total is 1, progress is displayed as "XX.XX %"
    """
    p = count / total if total > 0 else 0
    p = min(p, 1)
    columns, _ = os.get_terminal_size(0)
    if total > 1:
        l = int(log10(max(1, total)) + 1)
        columns -= 7 + 2 * l
        done = int(columns * p)
        print('\r\033[K[' + '#' * done + ' ' * (columns - done) + ']' +
              ' {:>{}} / {:>{}}'.format(int(count), l, int(total), l),
              end='', flush=True)
    else:
        columns -= 12
        done = int(columns * p)
        print('\r\033[K[' + '#' * done + ' ' * (columns - done) + ']' +
              ' {:>6.2f} %'.format(100 * p),
              end='', flush=True)


def init_progress():
    output_progress(0)


def finalize_progress():
    output_progress(1.0)
    print('')
