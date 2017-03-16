import os
import nose
import shutil
import unittest
import yaml
from breaker.utKit import utKit

from fundamentals import tools

su = tools(
    arguments={"settingsFile": None},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName="breaker"
)
arguments, settings, log, dbConn = su.setup()

# # load settings
# stream = file(
#     "/Users/Dave/.config/breaker/breaker.yaml", 'r')
# settings = yaml.load(stream)
# stream.close()

# SETUP AND TEARDOWN FIXTURE FUNCTIONS FOR THE ENTIRE MODULE
moduleDirectory = os.path.dirname(__file__)
utKit = utKit(moduleDirectory)
log, dbConn, pathToInputDir, pathToOutputDir = utKit.setupModule()
utKit.tearDownModule()

# load settings
stream = file(
    "/Users/Dave/.config/breaker/breaker.yaml", 'r')
settings = yaml.load(stream)
stream.close()

import shutil
try:
    shutil.rmtree(pathToOutputDir)
except:
    pass
# COPY INPUT TO OUTPUT DIR
shutil.copytree(pathToInputDir, pathToOutputDir)

# Recursively create missing directories
if not os.path.exists(pathToOutputDir):
    os.makedirs(pathToOutputDir)

# xt-setup-unit-testing-files-and-folders


class test_annotator(unittest.TestCase):

    def test_annotator_function(self):

        transients = {
            "fred": (62.26, 41.70),
            "bob": (61.26, -41.70),
            "joe": (60.36, 31.70),
            "sam": (60.45, 4.70),
            "arthur": (6.26, 41.0),
            "sarah": (10.26, 40.70),
            "jane": (160.26, -5.70)
        }

        from breaker.transients import annotator
        an = annotator(
            log=log,
            settings=settings,
            gwid="G211117"
        )
        transientNames, probs = an.annotate(transients)

    def test_annotator_prob_dist(self):

        from breaker.transients import annotator
        an = annotator(
            log=log,
            settings=settings,
            gwid="G211117"
        )
        transients = an._generate_probability_distribution()

    def test_annotator_function_exception(self):

        from breaker.transients import annotator
        try:
            this = annotator(
                log=log,
                settings=settings,
                fakeKey="break the code"
            )
            this.get()
            assert False
        except Exception, e:
            assert True
            print str(e)

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
