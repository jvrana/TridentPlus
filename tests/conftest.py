import pytest
import json
import os
from pydent import AqSession
from pydent.browser import Browser

import dill


@pytest.fixture(scope='session')
def datadir():
    here = os.path.abspath(os.path.dirname(__file__))
    dir = os.path.join(here, 'data')
    return dir



@pytest.fixture(scope="session")
def config():
    """
    Returns the config dictionary for live tests.
    """
    dir = os.path.dirname(os.path.abspath(__file__))

    config_path = os.path.join(dir, "secrets", "config.json.secret")
    with open(config_path, 'rU') as f:
        config = json.load(f)
    return config


@pytest.fixture(scope="session")
def session(config):
    """
    Returns a live aquarium connection.
    """
    return AqSession(**config)


@pytest.fixture(scope='function')
def browser(session):
    """
    Returns a live trident browser
    """
    filepath = os.path.join(datadir, 'browser.pkl')
    if os.path.isfile(filepath):
        with open(filepath, 'rb') as f:
            print("Loading browser from '{}'".fromat(filepath))
            browser = dill.load(f)
    else:
        browser = Browser(session)
        with open(filepath, 'wb') as f:
            dill.dump(browser, f)
    return browser
