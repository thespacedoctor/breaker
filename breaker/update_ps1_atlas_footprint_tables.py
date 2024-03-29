#!/usr/local/bin/python
# encoding: utf-8
"""
*Update the PS1 footprint table in breaker database and associate with GWs*

:Author:
    David Young

:Date Created:
    October 29, 2015
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
import numpy as np
from fundamentals import tools, times
from fundamentals.mysql import readquery, insert_list_of_dictionaries_into_database_tables, writequery
from astrocalc.coords import unit_conversion
from astrocalc.times import conversions
from HMpTy.mysql import add_htm_ids_to_mysql_database_table
from fundamentals.download import multiobject_download
from fundamentals.renderer import list_of_dictionaries
from fundamentals.mysql import writequery
from datetime import datetime, date, time
import codecs
import csv
import re
import io
import pymysql
from astrocalc.times import now as mjdnow
from fundamentals.renderer import list_of_dictionaries
from fundamentals.mysql import directory_script_runner
from .stats import survey_footprint


class update_ps1_atlas_footprint_tables():
    """
    *Update the PS1 footprint table in breaker database and associate with GWs*

    Metadata for each GW event should be found in the settings file and are used when associating the telescope pointings in the database with a GW event. For example, here are the metadata for the first GW burst:

    .. code-block:: yaml

        gravitational waves:
            G184098:
                human-name: GW150914
                mjd: 57279.90
                plot:
                    raRange: 48.  # CENTRAL WIDTH IN DEGREES
                    decRange: 45.  # CENTRAL HEIGHT IN DEGREES
                    centralCoordinate: [141., 0.0]
                mapPath: "/Users/Dave/Dropbox/notes/astronotes-wiki/projects/GW-G184098/Alternate-Skymap-Stats-For-G184098/maps/LALInference_skymap.fits"

    Both sky-area and the time-range are used in the gw-pointing associations.

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``updateNed`` -- do you want to update NED stream in database (can take a LONG time). Default *False*
        - ``updatePointings`` -- import new pointings from ATLAS and PS1 into the breaker database. Default *True*
        - ``updateAll`` -- update the PS1 database up to 21 days prior to the wave detecion. Default *False*
        - ``updateTransientEventTags`` -- tag transients in the database with their gravity event map probability region locations (tagged for each event). Default *True*

    **Usage:**

        .. code-block:: python

            from breaker import update_ps1_atlas_footprint_tables
            dbUpdater = update_ps1_atlas_footprint_tables(
                log=log,
                settings=settings,
                updateNed=False,
                updatePointings=True,
                updateAll=False,
                updateTransientEventTags=True
            )
            dbUpdater.get()
    """

    def __init__(
            self,
            log,
            settings=False,
            updateNed=False,
            updateAll=False,
            updatePointings=True,
            updateTransientEventTags=True
    ):
        self.log = log
        log.debug("instansiating a new 'update_ps1_atlas_footprint_tables' object")
        self.settings = settings
        self.updateNed = updateNed
        self.updateAll = updateAll
        self.updatePointings = updatePointings
        self.updateTransientEventTags = updateTransientEventTags

        # xt-self-arg-tmpx

        # SETUP THE VARIOUS DATABASE CONNECTIONS REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn, self.atlasDbConn, self.ps13piDbConn = db.get()

        return None

    def get(self):
        """
        *Import the new PS1 pointings and query NED for new data*

        This method:

            * Imports the new PS1 pointings from the PS1 database,
            * Imports the new ATLAS pointings from the ATLAS database,
            * attempts to label these pointings with the ID for an associated GW,
            * queries NED for new data covered by the sky-area of these pointings and
            * adds data to the NED stream database table

        See the ``update_ps1_atlas_footprint_tables`` class of usage info.

        **Return:**
            - None
        """
        self.log.debug('starting the ``get`` method')

        if self.updatePointings:
            if self.updateAll:
                recent = False
            else:
                recent = True
            self.import_new_ps1_pointings(recent=recent)
            self.import_new_atlas_pointings(recent=recent)
            # self.parse_panstarrs_nightlogs(updateAll=self.updateAll)
            self.label_pointings_with_gw_ids()
            # self.populate_ps1_subdisk_table()
        if self.updateNed:
            self.update_ned_database_table()
        self.update_gravity_event_annotations()

        self.log.debug('completed the ``get`` method')
        return None

    def import_new_ps1_pointings(
            self,
            recent=False):
        """
        *Import any new PS1 GW pointings from the ps1gw database into the ``ps1_pointings`` table of the Ligo-Virgo Waves database*

        **Key Arguments:**
            - ``recent`` -- only sync the most recent 1 month of data (speeds things up)

        **Return:**
            - None


         **Usage:**

            .. code-block:: python

                # IMPORT NEW PS1 POINTINGS FROM PS1 GW DATABASE INTO LIGO-VIRGO
                # WAVES DATABASE
                from breaker import update_ps1_atlas_footprint_tables
                dbUpdater = update_ps1_atlas_footprint_tables(
                    log=log,
                    settings=settings
                )
                dbUpdater.import_new_ps1_pointings()
        """
        self.log.debug('starting the ``import_new_ps1_pointings`` method')

        if recent:
            mjd = mjdnow(
                log=self.log
            ).get_mjd()
            recent = mjd - 21
            recent = " and mjd_obs > %(recent)s " % locals()
        else:
            recent = ""

        # SELECT ALL OF THE POINTING INFO REQUIRED FROM THE ps1gw DATABASE
        tables = ["ps1_warp_stack_diff_skycells",
                  "ps1_stack_stack_diff_skycells"]
        filenameMatch = ["ws", "ss"]
        for db in [self.ps1gwDbConn, self.ps13piDbConn]:
            for t, f in zip(tables, filenameMatch):

                sqlQuery = u"""
                    SELECT
                        imageid,
                        ppsub_input,
                        filename,
                        m.exptime exp_time,
                        TRUNCATE(mjd_obs, 8) mjd,
                        LEFT(fpa_filter, 1) AS filter,
                        IF(deteff_counts < 200,
                            m.zero_pt + m.deteff_magref+2.5*log(10,exptime),
                            m.zero_pt + m.deteff_magref + m.deteff_calculated_offset+2.5*log(10,exptime)) AS limiting_mag
                    FROM
                        tcs_cmf_metadata m
                        where filename like "%%.%(f)s.%%" %(recent)s
                """ % locals()

                rows = readquery(
                    log=self.log,
                    sqlQuery=sqlQuery,
                    dbConn=db,
                    quiet=False
                )

                # TIDY RESULTS BEFORE IMPORT
                entries = []

                converter = unit_conversion(
                    log=self.log
                )
                for row in rows:
                    e = {}
                    e["exp_time"] = row["exp_time"]
                    e["mjd"] = row["mjd"]
                    e["filter"] = row["filter"]
                    e["ps1_exp_id"] = row["imageid"]
                    e["limiting_mag"] = row["limiting_mag"]
                    e["filename"] = row["filename"]
                    e["skycell_id"] = (".").join(
                        row["filename"].split(".")[0:5])
                    e["target_image"] = row["ppsub_input"]
                    entries.append(e)

                extras = self.settings["ssh tunnels"][self.settings[
                    "database settings"]["ligo_virgo_waves"]["port"]]
                self.settings[
                    "database settings"]["ligo_virgo_waves"]["tunnel"] = extras

                # ADD THE NEW RESULTS TO THE ps1_pointings TABLE
                insert_list_of_dictionaries_into_database_tables(
                    dbConn=self.ligo_virgo_wavesDbConn,
                    log=self.log,
                    dictList=entries,
                    dbTableName=t,
                    uniqueKeyList=["filename"],
                    dateModified=False,
                    batchSize=2500,
                    replace=True,
                    dbSettings=self.settings[
                        "database settings"]["ligo_virgo_waves"]
                )

        print "PS1 skycells synced between `tcs_cmf_metadata` and `%(t)s` database tables" % locals()

        self.log.debug('completed the ``import_new_ps1_pointings`` method')
        return None

    def import_new_atlas_pointings(
            self,
            recent=False):
        """
        *Import any new ATLAS GW pointings from the atlas3 database into the ``atlas_pointings`` table of the Ligo-Virgo Waves database*

        **Key Arguments:**
            - ``recent`` -- only sync the most recent 1 month of data (speeds things up)

        **Return:**
            - None

         **Usage:**

            .. code-block:: python

                # IMPORT NEW ATLAS POINTINGS FROM ATLAS DATABASE INTO
                # LIGO-VIRGO WAVES DATABASE
                from breaker import update_ps1_atlas_footprint_tables
                dbUpdater = update_ps1_atlas_footprint_tables(
                    log=log,
                    settings=settings
                )
                dbUpdater.import_new_atlas_pointings()
        """
        self.log.debug('starting the ``import_new_atlas_pointings`` method')

        if recent:
            mjd = mjdnow(
                log=self.log
            ).get_mjd()
            recent = mjd - 21
            recent = " mjd_obs > %(recent)s " % locals()
        else:
            recent = "1=1"

        # SELECT ALL OF THE POINTING INFO REQUIRED FROM THE ATLAS DATABASE
        sqlQuery = u"""
            SELECT
                `dec` as `decDeg`,
                `exptime` as `exp_time`,
                `filter`,
                `mjd_obs` as `mjd`,
                `ra` as `raDeg`,
                if(mjd_obs<57855.0,mag5sig-0.75,mag5sig) as `limiting_magnitude`,
                `object` as `atlas_object_id` from atlas_metadata where %(recent)s and object like "TA%%" order by mjd_obs desc;
        """ % locals()
        rows = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.atlasDbConn,
            quiet=False
        )

        # IN ATLAS4 DATABASE atlas_metadata MOVED TO atlas_metadataddc WITH A
        # DIFFERENT SCHEMA - COLLECT THESE ROWS INSTEAD IF atlas_metadata EMPTY
        if not len(rows):
            recent = recent.replace("mjd_obs", "mjd")
            sqlQuery = u"""
            SELECT 
                `dec` AS `decDeg`,
                `texp` AS `exp_time`,
                `filt` AS `filter`,
                `mjd`,
                `ra` AS `raDeg`,
                mag5sig AS `limiting_magnitude`,
                `obj` AS `atlas_object_id`
            FROM
                atlas_metadataddc
            WHERE
                %(recent)s 
            ORDER BY mjd DESC
            """ % locals()
            rows = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.atlasDbConn,
                quiet=False
            )

        # TIDY RESULTS BEFORE IMPORT
        entries = list(rows)

        extras = self.settings["ssh tunnels"][self.settings[
            "database settings"]["ligo_virgo_waves"]["port"]]
        self.settings[
            "database settings"]["ligo_virgo_waves"]["tunnel"] = extras

        # ADD THE NEW RESULTS TO THE ps1_pointings TABLE
        insert_list_of_dictionaries_into_database_tables(
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            dictList=entries,
            dbTableName="atlas_pointings",
            uniqueKeyList=["raDeg", "decDeg", "mjd"],
            dateModified=False,
            batchSize=2500,
            replace=True,
            dbSettings=self.settings[
                "database settings"]["ligo_virgo_waves"]
        )

        # APPEND HTMIDs TO THE ps1_pointings TABLE
        add_htm_ids_to_mysql_database_table(
            raColName="raDeg",
            declColName="decDeg",
            tableName="atlas_pointings",
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            primaryIdColumnName="primaryId"
        )

        print "ATLAS pointings synced between `atlas_metadata` and `altas_pointings` database tables"

        self.log.debug('completed the ``import_new_atlas_pointings`` method')
        return None

    def label_pointings_with_gw_ids(
            self):
        """
        *Attempt to label the PS1 pointing with the GW IDs*

        The GW metadata used to associate PS1 pointings is taken from the settings file

        **Return:**
            - None

         **Usage:**

            .. code-block:: python

                # ATTEMPT TO LABEL PS1 POINTINGS IN DATABASE WITH A GW ID
                from breaker import update_ps1_atlas_footprint_tables
                dbUpdater = update_ps1_atlas_footprint_tables(
                    log=log,
                    settings=settings
                )
                dbUpdater.label_pointings_with_gw_ids()
        """
        self.log.debug('starting the ``label_pointings_with_gw_ids`` method')

        # WAVE METADATA FOUND IN SETTINGS FILE
        for wave in self.settings["gravitational waves"]:

            # UNPACK THE PLOT PARAMETERS FROM THE SETTINGS FILE
            centralCoordinate = self.settings["gravitational waves"][
                wave]["plot"]["centralCoordinate"]
            raRange = float(self.settings["gravitational waves"][
                wave]["plot"]["raRange"])
            decRange = float(self.settings["gravitational waves"][
                wave]["plot"]["decRange"])

            raMax = (centralCoordinate[0] + raRange / 2.) + 5.
            raMin = (centralCoordinate[0] - raRange / 2.) - 5.
            decMax = (centralCoordinate[1] + decRange / 2.) + 5.
            decMin = (centralCoordinate[1] - decRange / 2.) - 5.

            mjdLower = self.settings["gravitational waves"][
                wave]["mjd"] - 10.
            mjdUpper = self.settings["gravitational waves"][
                wave]["mjd"] + 21

            if raMin > 0. and raMax < 360.:
                raWhere = """(raDeg > %(raMin)s and raDeg < %(raMax)s)""" % locals(
                )
            elif raMin < 0.:
                raMin2 = raMin + 360.
                raWhere = """((raDeg > 0. and raDeg < %(raMax)s) or raDeg > %(raMin2)s)""" % locals(
                )
            elif raMax > 360.:
                raMax2 = raMax - 360.
                raWhere = """((raDeg > %(raMin)s and raDeg < 360.) or raDeg < %(raMax2)s)""" % locals(
                )

            decWhere = """(decDeg > %(decMin)s and  decDeg < %(decMax)s)""" % locals(
            )

            mjdWhere = "(mjd>%(mjdLower)s and mjd<%(mjdUpper)s)" % locals()

            sqlQuery = u"""
                update ps1_pointings set gw_id = "%(wave)s" where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is null
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )
            sqlQuery = u"""
                update ps1_pointings set gw_id = CONCAT(gw_id, " %(wave)s") where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is not null and gw_id not like "%%%(wave)s%%";
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )

            sqlQuery = u"""
                update atlas_pointings set gw_id = "%(wave)s" where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is null
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )
            sqlQuery = u"""
                update atlas_pointings set gw_id = CONCAT(gw_id, " %(wave)s") where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is not null and gw_id not like "%%%(wave)s%%";
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )

            mjdWhere = mjdWhere.replace("mjd", "mjd_registered")
            sqlQuery = u"""
                update ps1_nightlogs set gw_id = "%(wave)s" where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is null and type = "OBJECT"
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )
            sqlQuery = u"""
                update ps1_nightlogs set gw_id = CONCAT(gw_id, " %(wave)s") where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is not null and type = "OBJECT" and gw_id not like "%%%(wave)s%%";
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )

        sqlQuery = u"""
            select count(*) as count from ps1_pointings where gw_id is null;
        """ % locals()

        count = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
            quiet=False
        )[0]["count"]

        print "PS1 pointings labelled with their associated GW id"

        if count == 0:
            print "    Note all pointings have been labelled with GW ID"
        else:
            print "    %(count)s pointings remain unlabelled with a GW ID" % locals()

        self.log.debug('completed the ``label_pointings_with_gw_ids`` method')
        return None

    def populate_ps1_subdisk_table(
            self):
        """
        *Calculate 49 subdisks for each of the PS1 pointings (used to query NED in manageable sized batches) and add them to the ``ps1_pointings_subdisks`` table of the database*

        .. image:: http://i.imgur.com/y3G0aax.png
            :width: 600 px

        **Return:**
            - None

         **Usage:**

            .. code-block:: python

                # SPLIT PS1 POINTINGS INTO SUB-DISKS AND ADD TO LV DATABASE
                from breaker import update_ps1_atlas_footprint_tables
                dbUpdater = update_ps1_atlas_footprint_tables(
                    log=log,
                    settings=settings
                )
                dbUpdater.populate_ps1_subdisk_table()
        """
        self.log.debug(
            'completed the ````populate_ps1_subdisk_table`` method')

        # SELECT THE PS1 POINTINGS NEEDING SUBDISKS CALCULATED
        sqlQuery = u"""
            select ps1_exp_id, raDeg, decDeg from ps1_pointings where subdisks_calculated = 0 and raDeg is not null
        """ % locals()

        rows = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
            quiet=False
        )
        ps1PointNum = len(rows)

        # CALCULATE ALL OF THE SUBDISKS
        inserts = []
        expIds = []
        for row in rows:
            subDiskCoordinates = self._get_subdisk_parameters(
                row["raDeg"], row["decDeg"], 1.5)
            ps1_exp_id = row["ps1_exp_id"]
            expIds.append(ps1_exp_id)
            for i, c in enumerate(subDiskCoordinates):
                insert = {
                    "raDeg": c[0],
                    "decDeg": c[1],
                    "ps1_exp_id": ps1_exp_id,
                    "circleId": i + 1
                }
                inserts.append(insert)

        # ADD SUBDISKS TO DATABASE
        if len(inserts):

            extras = self.settings["ssh tunnels"][self.settings[
                "database settings"]["ligo_virgo_waves"]["port"]]
            self.settings[
                "database settings"]["ligo_virgo_waves"]["tunnel"] = extras

            insert_list_of_dictionaries_into_database_tables(
                dbConn=self.ligo_virgo_wavesDbConn,
                log=self.log,
                dictList=inserts,
                dbTableName="ps1_pointings_subdisks",
                uniqueKeyList=["ps1_exp_id", "circleId"],
                dateModified=False,
                batchSize=2500,
                replace=True,
                dbSettings=self.settings[
                    "database settings"]["ligo_virgo_waves"]
            )

            # UPDATE POINTINGS TABLE TO INDICATE SUBDISKS HAVE BEEN CALCULATED
            theseIds = ",".join(expIds)
            sqlQuery = u"""
                update ps1_pointings set subdisks_calculated = 1 where ps1_exp_id in (%(theseIds)s)
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )

        if ps1PointNum == 0:
            print "All PS1 pointings have been split into their 49 sub-disks" % locals()
        else:
            print "%(ps1PointNum)s new PS1 pointings have been split into 49 sub-disks - parameters added to the `ps1_pointings_subdisks` database table" % locals()

        # APPEND HTMIDs TO THE ps1_pointings_subdisks TABLE
        add_htm_ids_to_mysql_database_table(
            raColName="raDeg",
            declColName="decDeg",
            tableName="ps1_pointings_subdisks",
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            primaryIdColumnName="primaryId"
        )

        self.log.debug(
            'completed the ``populate_ps1_subdisk_table`` method')
        return None

    def _get_subdisk_parameters(
            self,
            raDeg,
            decDeg,
            radiusDeg):
        """
        *get subdisk parameters*

        **Key Arguments:**
            - ``raDeg`` -- the central ra of the ps1 pointing
            - ``decDeg`` -- the central dec of the ps1 pointing
            - ``radiusDeg`` -- the radius of the ps1 FOV

        **Return:**
            - ``subDiskCoordinates`` -- the coordinates for 49 subdisks covering the PS1 pointing
        """
        self.log.debug('starting the ``_get_subdisk_parameters`` method')

        import math
        footprintCoords = (raDeg, decDeg)

        shifts = [
            (0, 0),
            (0, math.sqrt(3.) / 2.),
            (3. / 4., math.sqrt(3.) / 4.),
            (3. / 4., -math.sqrt(3.) / 4.),
            (0, -math.sqrt(3.) / 2.),
            (-3. / 4., -math.sqrt(3.) / 4.),
            (-3. / 4., math.sqrt(3.) / 4.)
        ]

        subDiskCoordinates = []
        count = 0
        radius2 = radiusDeg / 2
        for s in shifts:
            x1 = footprintCoords[0] + s[0] * radiusDeg
            y1 = footprintCoords[1] + s[1] * radiusDeg
            for ss in shifts:
                count += 1
                x2 = x1 + ss[0] * radius2
                y2 = y1 + ss[1] * radius2
                self.log.debug("""%(count)s: %(x2)s, %(y2)s""" % locals())
                subDiskCoordinates.append((x2, y2))

        self.log.debug('completed the ``_get_subdisk_parameters`` method')
        return subDiskCoordinates

    def update_ned_database_table(
            self):
        """
        *Use Sherlock & Neddy to query NED and update the catalogues database for previously unseen/stale PS1 footprint areas*

        **Return:**
            - None

        **Usage:**

            .. code-block:: python

                # UPDATE THE NED STREAM FOR NEW PS1 FOOTPRINTS
                from breaker import update_ps1_atlas_footprint_tables
                dbUpdater = update_ps1_atlas_footprint_tables(
                    log=log,
                    settings=settings
                )
                dbUpdater.update_ned_database_table()
        """
        self.log.debug('starting the ``update_ned_database_table`` method')

        from sherlock.update_ned_stream import update_ned_stream

        numDisksToConesearch = 100
        rowCount = 100

        while rowCount > 0:

            sqlQuery = u"""
                select primaryId, raDeg as "ra", decDeg as "dec", htm16ID from ps1_pointings_subdisks where nedQueried = 0 limit %(numDisksToConesearch)s
            """ % locals()
            rows = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
                quiet=False
            )
            rowCount = len(rows)
            ids = []
            ids[:] = [str(row["primaryId"]) for row in rows]
            ids = ",".join(ids)

            if rowCount > 0:
                print "Selecting the next %(rowCount)s subdisks areas to conesearch against NED from the `ps1_pointings_subdisks` table" % locals()
            else:
                print "NED stream is up-to-date, no queries required" % locals()

            update_ned_stream(
                log=self.log,
                cataloguesDbConn=self.cataloguesDbConn,
                settings=self.settings,
                transientsMetadataList=rows
            ).get()

            if len(ids):
                sqlQuery = u"""
                    update ps1_pointings_subdisks set nedQueried = 1 where primaryId in (%(ids)s)
                """ % locals()
                writequery(
                    log=self.log,
                    sqlQuery=sqlQuery,
                    dbConn=self.ligo_virgo_wavesDbConn,
                )

            sqlQuery = u"""
                select count(*) as count from ps1_pointings_subdisks where nedQueried = 0
            """ % locals()
            count = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
                quiet=False
            )
            count = count[0]["count"]

            if rowCount > 0:
                print "NED stream updated for %(rowCount)s PS1 pointing sub-disks (%(count)s to go)" % locals()
                print "-----\n\n"

        self.log.debug('completed the ``update_ned_database_table`` method')
        return None

    def parse_panstarrs_nightlogs(
            self,
            updateAll=False):
        """*download and parse the ps1 night logs from the range of time a wave survey campaign is active*

        The night-log data is added to the ps1_nightlogs table

        **Key Arguments:**
            - ``updateAll`` -- update all of the PS1 nightlogs. This will take sometime, the default is to lift the logs from the last 7 days. Default *False*.

        **Return:**
            - None

        **Usage:**
            ..  todo::

                - add usage info
                - create a sublime snippet for usage
                - update package tutorial if needed

            .. code-block:: python

                usage code

        """
        self.log.debug('starting the ``parse_panstarrs_nightlogs`` method')

        # CONVERTER TO CONVERT MJD TO DATE
        converter = conversions(
            log=self.log
        )

        createStatement = """
CREATE TABLE `ps1_nightlogs` (
  `primaryId` bigint(20) NOT NULL AUTO_INCREMENT COMMENT 'An internal counter',
  `airm` double DEFAULT NULL,
  `comments` varchar(200) DEFAULT NULL,
  `decDeg` double DEFAULT NULL,
  `etime` double DEFAULT NULL,
  `f` varchar(10) DEFAULT NULL,
  `filesetID` varchar(100) DEFAULT NULL,
  `raDeg` double DEFAULT NULL,
  `telescope_pointing` varchar(200) DEFAULT NULL,
  `time_registered` datetime DEFAULT NULL,
  `type` varchar(100) DEFAULT NULL,
  `dateCreated` datetime DEFAULT CURRENT_TIMESTAMP,
  `dateLastModified` datetime DEFAULT CURRENT_TIMESTAMP,
  `updated` varchar(45) DEFAULT '0',
  PRIMARY KEY (`primaryId`),
  UNIQUE KEY `filesetid` (`filesetID`)
) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=latin1;
"""

        from astrocalc.times import now
        mjdNow = now(
            log=self.log
        ).get_mjd()

        # WAVE METADATA FOUND IN SETTINGS FILE
        for wave in self.settings["gravitational waves"]:
            # GIVE A 3 DAY WINDOW EITHER SIDE OF WAVE TIME-RANGE
            mjdLower = int(self.settings["gravitational waves"][
                wave]["mjd"] - 10. - 3.)
            mjdUpper = int(self.settings["gravitational waves"][
                wave]["mjd"] + 21. + 3.)

            if updateAll == False:
                if mjdUpper < mjdNow - 7.:
                    continue
                if mjdUpper > mjdNow:
                    mjdUpper = int(mjdNow)
                if mjdLower < mjdNow - 7.:
                    mjdLower = int(mjdNow - 7.)

            # METRIC NIGHT LOGS FOR EACH NIGHT FOUND AT A URL SIMILAR TO :
            # "http://ipp0022.ifa.hawaii.edu/ps1sc/metrics/2016-12-14/index.html"
            urls = []
            for i in range(mjdUpper - mjdLower + 3):
                mjd = i + mjdLower
                utDate = converter.mjd_to_ut_datetime(
                    mjd=mjd,
                    sqlDate=False,
                    datetimeObject=True
                )
                utDate = utDate.strftime("%Y-%m-%d")
                urls.append("https://ps1wiki.ifa.hawaii.edu/data/metrics/%(utDate)s/index.html" % locals(
                ))

            localUrls = multiobject_download(
                urlList=urls,
                downloadDirectory="/tmp",
                log=self.log,
                timeStamp=True,
                timeout=180,
                concurrentDownloads=2,
                resetFilename=False,
                credentials=False,  # { 'username' : "...", "password", "..." }
                longTime=True,
                indexFilenames=False
            )

            for url in localUrls:
                if not url:
                    continue
                pathToReadFile = url
                try:
                    self.log.debug("attempting to open the file %s" %
                                   (pathToReadFile,))
                    readFile = codecs.open(
                        pathToReadFile, encoding='utf-8', mode='r')
                    thisData = readFile.read()
                    readFile.close()
                except IOError, e:
                    message = 'could not open the file %s' % (pathToReadFile,)
                    self.log.critical(message)
                    raise IOError(message)
                readFile.close()

                regex = re.compile(r'<pre>\s*# (filesetID.*?)</pre>', re.S)
                matchObject = re.finditer(
                    regex,
                    thisData
                )

                for match in matchObject:
                    csvReader = csv.DictReader(
                        io.StringIO(match.group(1)), delimiter='|')
                    nightLog = []
                    for row in csvReader:
                        cleanDict = {}
                        for k, v in row.iteritems():
                            cleanDict[k.strip().replace(" ", "_")] = v.strip()
                        if "telescope_pointing" in cleanDict:
                            cleanDict["raDeg"] = cleanDict["telescope_pointing"].split()[
                                0]
                            cleanDict["decDeg"] = cleanDict["telescope_pointing"].split()[
                                1]
                        if "time_registered" in cleanDict:
                            cleanDict["time_registered"] = cleanDict[
                                "time_registered"].replace("Z", "")
                        nightLog.append(cleanDict)

                dataSet = list_of_dictionaries(
                    log=self.log,
                    listOfDictionaries=nightLog
                )
                # Recursively create missing directories
                if not os.path.exists("/tmp/ps1_nightlogs"):
                    os.makedirs("/tmp/ps1_nightlogs")
                mysqlData = dataSet.mysql(
                    tableName="ps1_nightlogs", filepath="/tmp/ps1_nightlogs/ps1_nightlog_%(utDate)s.sql" % locals(), createStatement=createStatement)

                directory_script_runner(
                    log=self.log,
                    pathToScriptDirectory="/tmp/ps1_nightlogs",
                    databaseName=self.settings["database settings"][
                        "ligo_virgo_waves"]["db"],
                    loginPath=self.settings["database settings"][
                        "ligo_virgo_waves"]["loginPath"],
                    successRule="delete",
                    failureRule="failed"
                )

        # APPEND HTMIDs TO THE ps1_pointings TABLE
        add_htm_ids_to_mysql_database_table(
            raColName="raDeg",
            declColName="decDeg",
            tableName="ps1_nightlogs",
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            primaryIdColumnName="primaryId"
        )

        sqlQuery = """update ps1_nightlogs set mjd_registered = ((UNIX_TIMESTAMP(time_registered)/ 86400.0)+ 2440587.5-2400000.5) where mjd_registered is null;""" % locals()
        writequery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        self.log.debug('completed the ``parse_panstarrs_nightlogs`` method')
        return None

    def update_gravity_event_annotations(
            self):
        """*update gravity event annotations*

        **Key Arguments:**
            # -

        **Return:**
            - None

        **Usage:**
            ..  todo::

                - add usage info
                - create a sublime snippet for usage
                - write a command-line tool for this method
                - update package tutorial with command-line tool info if needed

            .. code-block:: python

                usage code

        """
        self.log.debug(
            'starting the ``update_gravity_event_annotations`` method')

        from breaker.transients import annotator

        # CREATE THE ANNOTATION HELPER TABLES IF THEY DON"T EXIST
        moduleDirectory = os.path.dirname(__file__)
        mysql_scripts = moduleDirectory + "/resources/mysql"

        self.log.info(
            'creating tables for transient and exposure annotations')
        for db in ["ps1gw", "ps13pi", "atlas"]:
            directory_script_runner(
                log=self.log,
                pathToScriptDirectory=mysql_scripts + "/create_tables",
                databaseName=self.settings["database settings"][db]["db"],
                loginPath=self.settings["database settings"][db]["loginPath"],
                waitForResult=True,
                successRule=None,
                failureRule=None
            )
        for db in ["ligo_virgo_waves"]:
            directory_script_runner(
                log=self.log,
                pathToScriptDirectory=mysql_scripts + "/ligo_virgo_waves/create_tables",
                databaseName=self.settings["database settings"][db]["db"],
                loginPath=self.settings["database settings"][db]["loginPath"],
                waitForResult=True,
                successRule=None,
                failureRule=None
            )

        # UPDATE THE TABLE WITH THE METADATA OF EACH GRAVITY EVENT
        self.log.info(
            'updating GW metadata table')
        sqlQuery = ""
        for g in self.settings["gravitational waves"]:
            h = self.settings["gravitational waves"][g]["human-name"]
            m = self.settings["gravitational waves"][g]["mjd"]
            cmd = """insert ignore into tcs_gravity_events (`gracedb_id`, `gravity_event_id`, `mjd`) VALUES ("%(g)s", "%(h)s", %(m)s) on duplicate key update mjd=%(m)s;\n""" % locals(
            )
            sqlQuery += cmd
        for db in [self.atlasDbConn, self.ps1gwDbConn, self.ps13piDbConn]:
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=db
            )
        sqlQuery = sqlQuery.replace("tcs_gravity_events", "gravity_events")
        writequery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
        )
        self.log.info(
            'populating transient tables with placeholders before transient annotations')
        for db in ["ps1gw", "ps13pi", "atlas"]:
            directory_script_runner(
                log=self.log,
                pathToScriptDirectory=mysql_scripts + "/annotate",
                databaseName=self.settings["database settings"][db]["db"],
                loginPath=self.settings["database settings"][db]["loginPath"],
                waitForResult=True,
                successRule=None,
                failureRule=None
            )
        for db in ["ligo_virgo_waves"]:
            directory_script_runner(
                log=self.log,
                pathToScriptDirectory=mysql_scripts + "/ligo_virgo_waves/annotate",
                databaseName=self.settings["database settings"][db]["db"],
                loginPath=self.settings["database settings"][db]["loginPath"],
                waitForResult=True,
                successRule=None,
                failureRule=None
            )

        dbDict = {
            "ps1gw": self.ps1gwDbConn,
            "atlas": self.atlasDbConn,
            "ps13pi": self.ps13piDbConn,
            "ligo_virgo_waves": self.ligo_virgo_wavesDbConn
        }

        mjd = mjdnow(
            log=self.log
        ).get_mjd()

        for db in dbDict.keys():

            print "Adding transient/exposure annotations to %(db)s database" % locals()

            for g in self.settings["gravitational waves"]:

                m = self.settings["gravitational waves"][g]["mjd"]
                if m < mjd - 22:
                    continue

                h = self.settings["gravitational waves"][g]["human-name"]
                print "Annotating new transient/exposure associated with gravity event %(h)s" % locals()
                m = self.settings["gravitational waves"][g]["mjd"]
                mapPath = self.settings["gravitational waves"][g]["mapPath"]
                mapName = os.path.basename(mapPath)

                thisDbConn = dbDict[db]

                if thisDbConn in [self.ps1gwDbConn, self.ps13piDbConn]:

                    self.log.info(
                        """selecting transients from `%(db)s` database""" % locals())

                    # DELETE STALE MAP ANNOTATIONS
                    sqlQuery = u"""update tcs_gravity_event_annotations set map_name = null, enclosing_contour = null where gracedb_id = "%(g)s" and map_name != "%(mapName)s";""" % locals(
                    )
                    writequery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                    )

                    sqlQuery = u"""
                        SELECT
                            a.transient_object_id, a.gracedb_id, t.ra_psf, t.dec_psf
                        FROM
                            tcs_transient_objects t,
                            tcs_gravity_event_annotations a
                        WHERE
                            a.transient_object_id = t.id
                                AND t.detection_list_id != 0
                                AND (a.map_name !=
                                     "%(mapName)s"  or a.map_name is null)
                                AND a.gracedb_id="%(g)s";
                    """ % locals()

                    rows = readquery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                        quiet=False
                    )

                    transients = {}
                    for r in rows:
                        transients[r["transient_object_id"]] = (
                            r["ra_psf"], r["dec_psf"])

                    self.log.info(
                        """annotating selected transients from `%(db)s` database""" % locals())
                    an = annotator(
                        log=self.log,
                        settings=self.settings,
                        gwid=g
                    )
                    transientNames, probs = an.annotate(transients)

                if thisDbConn in [self.atlasDbConn]:
                    # DELETE STALE MAP ANNOTATIONS
                    sqlQuery = u"""update tcs_gravity_event_annotations set map_name = null, enclosing_contour = null where gracedb_id = "%(g)s" and map_name != "%(mapName)s";""" % locals(
                    )
                    writequery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                    )

                    self.log.info(
                        """selecting transients from `%(db)s` database""" % locals())
                    sqlQuery = u"""
                        SELECT
                            a.transient_object_id, a.gracedb_id, t.ra, t.dec
                        FROM
                            atlas_diff_objects t,
                            tcs_gravity_event_annotations a
                        WHERE
                            a.transient_object_id = t.id
                                AND t.detection_list_id != 0
                                AND (a.map_name !=
                                     "%(mapName)s"  or a.map_name is null)
                                AND a.gracedb_id="%(g)s";
                    """ % locals()

                    rows = readquery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                        quiet=False
                    )

                    transients = {}
                    for r in rows:
                        transients[r["transient_object_id"]] = (
                            r["ra"], r["dec"])

                    self.log.info(
                        """annotating selected transients from `%(db)s` database""" % locals())
                    an = annotator(
                        log=self.log,
                        settings=self.settings,
                        gwid=g
                    )
                    transientNames, probs = an.annotate(transients)

                if thisDbConn in [self.ligo_virgo_wavesDbConn]:

                    print "Annotating PS1 skycells associated with gravity event %(h)s" % locals()

                    # DELETE STALE MAP ANNOTATIONS
                    sqlQuery = u"""update ps1_skycell_gravity_event_annotations set prob_coverage = null, map_name = null where gracedb_id = "%(g)s" and map_name != "%(mapName)s";""" % locals(
                    )
                    writequery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                    )

                    # PANSTARRS SKYCELLS
                    self.log.info(
                        """selecting PS1 skycells from `%(db)s` database""" % locals())
                    sqlQuery = u"""
                        SELECT 
                                a.skycell_id, a.gracedb_id, t.raDeg, t.decDeg
                            FROM
                                ps1_skycell_map t,
                                ps1_skycell_gravity_event_annotations a
                            WHERE
                                a.skycell_id = t.skycell_id
                                AND (a.map_name != "%(mapName)s"  or a.map_name is null)
                                AND a.gracedb_id="%(g)s"; 
                    """ % locals()
                    rows = readquery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                        quiet=False
                    )

                    if len(rows):
                        exposures = {}
                        for r in rows:
                            exposures[r["skycell_id"]] = (
                                r["raDeg"], r["decDeg"])

                        stats = survey_footprint(
                            log=self.log,
                            settings=self.settings,
                            gwid=g
                        )
                        self.log.info(
                            """annotating PS1 skycells from `%(db)s` database""" % locals())
                        exposureIDs, probs = stats.annotate_exposures(
                            exposures=exposures,
                            pointingSide=0.4
                        )

                        dataList = []
                        for p, t in zip(probs, exposureIDs):
                            dataList.append({
                                "skycell_id": t,
                                "prob_coverage": p,
                                "gracedb_id": g,
                                "map_name": mapName
                            })
                        tableName = "ps1_skycell_gravity_event_annotations"

                        dataSet = list_of_dictionaries(
                            log=self.log,
                            listOfDictionaries=dataList,
                            reDatetime=re.compile(
                                '^[0-9]{4}-[0-9]{2}-[0-9]{2}T')
                        )
                        # RECURSIVELY CREATE MISSING DIRECTORIES
                        if not os.path.exists("/tmp/mysqlinsert/%(db)s" % locals()):
                            os.makedirs("/tmp/mysqlinsert/%(db)s" % locals())
                        now = datetime.now()
                        now = now.strftime("%Y%m%dt%H%M%S%f")
                        mysqlData = dataSet.mysql(
                            tableName=tableName, filepath="/tmp/mysqlinsert/%(db)s/%(now)s.sql" % locals(), createStatement=False)

                    print "Annotating ATLAS exposures associated with gravity event %(h)s" % locals()

                    # DELETE STALE MAP ANNOTATIONS
                    sqlQuery = u"""update atlas_exposure_gravity_event_annotations set prob_coverage = null, map_name = null  where gracedb_id = "%(g)s" and map_name != "%(mapName)s";""" % locals(
                    )
                    writequery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                    )

                    # ATLAS EXPOSURES
                    self.log.info(
                        """selecting atlas exposures from `%(db)s` database""" % locals())
                    sqlQuery = u"""
                        SELECT 
                                atlas_object_id, gracedb_id, raDeg, decDeg
                            FROM
                                atlas_exposure_gravity_event_annotations
                            WHERE
                                (map_name != "%(mapName)s"  or map_name is null)
                                AND gracedb_id="%(g)s"; 
                    """ % locals()
                    rows = readquery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=thisDbConn,
                        quiet=False
                    )

                    if len(rows):
                        total = len(rows)
                        print "%(total)s ATLAS exposures still need annotating, annotating these now" % locals()
                        exposures = {}
                        for r in rows:
                            exposures[r["atlas_object_id"]] = (
                                r["raDeg"], r["decDeg"])

                        stats = survey_footprint(
                            log=self.log,
                            settings=self.settings,
                            gwid=g
                        )
                        self.log.info(
                            """annotating atlas exposures from `%(db)s` database""" % locals())
                        exposureIDs, probs = stats.annotate_exposures(
                            exposures=exposures,
                            pointingSide=5.46
                        )

                        dataList = []
                        for p, t in zip(probs, exposureIDs):
                            dataList.append({
                                "atlas_object_id": t,
                                "prob_coverage": p,
                                "gracedb_id": g,
                                "map_name": mapName
                            })
                        tableName = "atlas_exposure_gravity_event_annotations"

                        dataSet = list_of_dictionaries(
                            log=self.log,
                            listOfDictionaries=dataList,
                            reDatetime=re.compile(
                                '^[0-9]{4}-[0-9]{2}-[0-9]{2}T')
                        )
                        # RECURSIVELY CREATE MISSING DIRECTORIES
                        if not os.path.exists("/tmp/mysqlinsert/%(db)s" % locals()):
                            os.makedirs("/tmp/mysqlinsert/%(db)s" % locals())
                        now = datetime.now()
                        now = now.strftime("%Y%m%dt%H%M%S%f")
                        mysqlData = dataSet.mysql(
                            tableName=tableName, filepath="/tmp/mysqlinsert/%(db)s/%(now)s.sql" % locals(), createStatement=False)

                if thisDbConn not in [self.ligo_virgo_wavesDbConn]:
                    dataList = []

                    for p, t in zip(probs, transientNames):
                        # if p < 91:
                        if p < 101:
                            dataList.append({
                                "transient_object_id": t,
                                "enclosing_contour": p,
                                "gracedb_id": g,
                                "map_name": mapName
                            })
                    tableName = "tcs_gravity_event_annotations"

                    dataSet = list_of_dictionaries(
                        log=self.log,
                        listOfDictionaries=dataList,
                        reDatetime=re.compile('^[0-9]{4}-[0-9]{2}-[0-9]{2}T')
                    )
                    # RECURSIVELY CREATE MISSING DIRECTORIES
                    if not os.path.exists("/tmp/mysqlinsert/%(db)s" % locals()):
                        os.makedirs("/tmp/mysqlinsert/%(db)s" % locals())
                    now = datetime.now()
                    now = now.strftime("%Y%m%dt%H%M%S%f")
                    mysqlData = dataSet.mysql(
                        tableName=tableName, filepath="/tmp/mysqlinsert/%(db)s/%(now)s.sql" % locals(), createStatement=False)

                exists = os.path.exists("/tmp/mysqlinsert/%(db)s" % locals())
                if exists:

                    # IMPORT TO DATABASE
                    directory_script_runner(
                        log=self.log,
                        pathToScriptDirectory="/tmp/mysqlinsert/%(db)s" % locals(
                        ),
                        databaseName=self.settings[
                            "database settings"][db]["db"],
                        loginPath=self.settings[
                            "database settings"][db]["loginPath"],
                        waitForResult=True,
                        successRule="delete",
                        failureRule="failed"
                    )

        self.log.debug(
            'completed the ``update_gravity_event_annotations`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
