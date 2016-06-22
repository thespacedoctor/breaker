import os
import nose
import shutil
import yaml
from breaker.utKit import utKit

from fundamentals import tools

su = tools(
    arguments={"settingsFile": None},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName="breaker",
    tunnel=False
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


class test_listen():

    def test_listen_function(self):

        from breaker.gracedb import listen
        this = listen(
            log=log,
            settings=settings,
            label="EM_READY",
            farThreshold=1e-7,
            startMJD=56658.0,
            endMJD=69807.0
        )
        this.get_maps()

    def test_listen_daemon_function(self):

        from breaker.gracedb import listen
        this = listen(
            log=log,
            settings=settings,
            label="EM_READY",
            farThreshold=1e-7,
            daemon=True
        )
        this.get_maps()

    # def test_listen_function_exception(self):

    #     from breaker import listen
    #     try:
    #         this = listen(
    #             log=log,
    #             settings=settings,
    #             fakeKey="break the code"
    #         )
    #         this.get()
    #         assert False
    #     except Exception, e:
    #         assert True
    #         print str(e)

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
