import logging
import coloredlogs

class Logger:

    def __init__(self, name, level = 10):
        fields = coloredlogs.DEFAULT_FIELD_STYLES
        fields['levelname']['color'] = 'magenta'
        logging.captureWarnings(True)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        handler = logging.StreamHandler()  # For console output
        handler.setLevel(level)

        format = '%(asctime)s %(levelname)-7s %(name)s:%(lineno)d: %(message)s'
        coloredlogs.install(level=level, logger=logger, fmt=format, datefmt='%H:%M:%S', field_styles=fields)

        self.LOGGER = logger
        self.HANDLER = handler