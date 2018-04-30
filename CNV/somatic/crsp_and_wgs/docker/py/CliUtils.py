import logging

def setup_logging(log_filename, is_verbose_stdout=True):
    """ Create logging for stdout and to a file.
    """
    # Create a basic logger to a file
    loggingFormat = '%(asctime)s %(levelname)s [%(name)s:%(lineno)d] %(message)s'
    logging.basicConfig(filename=log_filename, level=logging.INFO, format=loggingFormat)


    # Add a console logger to the root logger, which means that all loggers generated will have the console dump.
    #    Output on the console will be the same as what is in the log file.
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARN)
    formatter = logging.Formatter(loggingFormat)
    ch.setFormatter(formatter)

    if is_verbose_stdout:
        ch.setLevel(logging.INFO)

    logging.getLogger('').addHandler(ch)

def count_lines(filename):
    """ http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    :param filename:
    :return:
    """
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1