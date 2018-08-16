import os
import nose
import shutil
import yaml
import unittest
from breaker import cl_utils
from breaker.utKit import utKit
from breaker.stats import survey_footprint

from fundamentals import tools, times

su = tools(
    arguments={},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName="breaker"
)
arguments, settings, log, dbConn = su.setup()

# load settings
stream = file(
    "/Users/Dave/.config/breaker/breaker.yaml", 'r')
settings = yaml.load(stream)
stream.close()

# SETUP AND TEARDOWN FIXTURE FUNCTIONS FOR THE ENTIRE MODULE
moduleDirectory = os.path.dirname(__file__)
utKit = utKit(moduleDirectory)
log, dbConn, pathToInputDir, pathToOutputDir = utKit.setupModule()
utKit.tearDownModule()


class test_survey_footprint(unittest.TestCase):

    def test_survey_footprint_function(self):

        from breaker.stats import survey_footprint
        stats = survey_footprint(
            log=log,
            settings=settings,
            gwid="G184098"
        )
        stats.get()

    def test_annotator_exposure_coverage(self):

        exposures = {
            "fred": (62.26, 41.70),
            "bob": (61.26, -41.70),
            "joe": (60.36, 31.70),
            "sam": (60.45, 4.70),
            "arthur": (6.26, 41.0),
            "sarah": (10.26, 40.70),
            "jane": (160.26, -5.70)
        }

        from breaker.stats import survey_footprint
        stats = survey_footprint(
            log=log,
            settings=settings,
            gwid="G211117"
        )
        exposureAreas = stats.annotate_exposures(
            exposures=exposures,
            pointingSide=5.46
        )

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
