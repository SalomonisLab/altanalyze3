import logging


def setup_logger(logger, log_level, log_format=None):
    log_format = "%(processName)12s (%(asctime)s): %(message)s" if log_format is None else log_format
    for log_handler in logger.handlers:
        logger.removeHandler(log_handler)
    for log_filter in logger.filters:
        logger.removeFilter(log_filter)
    logging.basicConfig(level=log_level, format=log_format)