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


class test_projections():

    def test_projections_function(self):

        from breaker.plots import projections
        this = projections(
            log=log,
            gwid="GW151226",
            healpixPath=pathToInputDir + "/GW151226-LALInference_skymap_final.fits",
            projection="mollweide",
            outputDirectory=pathToOutputDir
        )
        this.get()

    # def test_projections_function2(self):

    #     from breaker.plots import projections
    #     this = projections(
    #         log=log,
    #         gwid="GW151226",
    #         healpixPath=pathToInputDir + "/GW151226-LALInference_skymap_final.fits",
    #         projection="aitoff",
    #         outputDirectory=pathToOutputDir
    #     )
    #     this.get()

    # def test_projections_function3(self):

    #     from breaker.plots import projections
    #     this = projections(
    #         log=log,
    #         gwid="GW151226",
    #         healpixPath=pathToInputDir + "/GW151226-LALInference_skymap_final.fits",
    #         projection="hammer",
    #         outputDirectory=pathToOutputDir
    #     )
    #     this.get()

    # def test_projections_function4(self):

    #     from breaker.plots import projections
    #     this = projections(
    #         log=log,
    #         gwid="GW151226",
    #         healpixPath=pathToInputDir + "/GW151226-LALInference_skymap_final.fits",
    #         projection="lambert",
    #         outputDirectory=pathToOutputDir
    #     )
    #     this.get()

    # def test_projections_function5(self):

    #     from breaker.plots import projections
    #     this = projections(
    #         log=log,
    #         gwid="GW151226",
    #         healpixPath=pathToInputDir + "/GW151226-LALInference_skymap_final.fits",
    #         projection="polar",
    #         outputDirectory=pathToOutputDir
    #     )
    #     this.get()

    # def test_projections_function6(self):

    #     from breaker.plots import projections
    #     this = projections(
    #         log=log,
    #         gwid="GW151226",
    #         healpixPath=pathToInputDir + "/GW151226-LALInference_skymap_final.fits",
    #         projection="rectilinear",
    #         outputDirectory=pathToOutputDir
    #     )
    #     this.get()

    def test_projections_function_exception(self):

        from breaker.plots import projections
        try:
            this = projections(
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
