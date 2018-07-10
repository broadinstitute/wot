import time

start_time = time.time()

def time_verbose(message, update=False, prefix='', suffix=''):
    """Prints a message with a timestamp. For performance measuring purposes"""
    prefix = "\r\033[K" + prefix
    if update:
        suffix += "\r\033[1A"
    print("{}{:>50} [ {:>10.3f} ]{}".format(
        prefix, message,
        time.time() - start_time,
        suffix))
