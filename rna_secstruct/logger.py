"""
Logging Setup Module

This module provides functions to set up logging configurations for Python applications.
It includes functionalities to set up root-level logging, application-level logging, and
to get logger instances with specific module names.

Functions:
    - setup_logging(file_name: Optional[str] = None) -> logging.Logger
        Set up the root logging configuration with optional file logging.

    - get_logger(module_name: str = "") -> logging.Logger
        Get a logger instance with the specified module name.
"""

import logging
import sys
from typing import Optional

# logging #####################################################################

APP_LOGGER_NAME = "rna-secstruct"


def setup_logging(file_name: Optional[str] = None) -> logging.Logger:
    """
    Set up logging configuration.

    Args:
        file_name (str, optional): The name of the log file. If provided, logs will be
        written to this file.

    Returns:
        None
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)  # Set the root logger level

    # Create a stream handler for output to console
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)  # Set the desired level for console output
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        root_logger.addHandler(fh)

    return root_logger


def get_logger(module_name: str = "") -> logging.Logger:
    """
    Get a logger instance with the specified module name.

    Args:
        module_name (str): The name of the module to be included in the logger name.

    Returns:
        logging.Logger: A logger instance with the specified module name.

    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
